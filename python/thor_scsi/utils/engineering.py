"""Handle engineering tolerances

The task is split up in two different entites:

    1. create commands that describe a random distribution

    2. create commands that will apply a variation to a given
       property of the element.

The first a sequence of commands are used to describe the varations
(see :class:`EngineeringDistributionCommand`).

Using a :class:`ǹumpy.random.Generator` these commands are used to
create a sequence of commands that describe how each element property
should be changed precisely. These commands are then executed by
:func:`àpply_factors`.

This functionallity is intended to support engineering tolerance
studies. The user is expected to create commands describing how
the user wants to variy the different parameters. If the identical
variation can be used for different elements,
:func:`create_distribution_commands` will be handy.

The specific commands that describe what specific changes will be
applied to the different elements (and its properties) can then be
created by :func:`create_commands_for_randomised_state`.

the function :func:`apply_factors` can then be used to apply these
commands to all the different elements. This function has been split
from the functionality of creating the specific command so that the
user can inspect or store this specific set of commands.


The pseudo random number generator to be used is under users control.
It has to conform to the :class:`ǹumpy.random.Generator` API. Thus
the user can control the seed and similar aspects by the pseudo random
number generator's API.

If replay of random number generation is of interest, please note
that the sequence of engineering distribution commands supplied to
:func:`create_commands_for_randomised_state` is important. It will
make calls to the random number generator in the order as the
commands are applied. Please also note that changing the vector size
has an impact. Thus it is currently discussed if the element properties,
that are vectors should benignly ignore elements of the scale or
add vector that are beyond the size of the vector to work on.
This would allow e.g. to apply the same engineering values to the
multipoles if e.g. 4 or 10 multipole are used, as long as the
same vector length is used in the instances of
:class:`EngineeringDistributionCommand`.
"""
import thor_scsi.lib as tslib
import numpy as np
import numpy.random

# import xarray as xr
from dataclasses import dataclass
from typing import Sequence
import copy


@dataclass
class EngineeringCommand:
    """Command for applying a specific set or value for exeucting engineering studies

    Todo:
       * consider json import / export: thus long term storage could be feasible
    """

    # Method name on how to retrieve the object and to set it
    # calls for a better abstraction
    set_method_name: str
    get_method_name: str

    # I assume that I can scale all factors
    scale_factors: np.ndarray
    add_factors: np.ndarray

    # on what to works
    element_index: int


@dataclass
class EngineeringDistributionCommand:
    """Command for creating engineering command abstract class

    Todo:
        How to handle vector length:
         Let's compute too many factors and ignore if there
         are too many?

         Shall that be done by the appropriate element?

       consider json import / export: thus long term storage could be feasible
    """

    # Method name on how to set the string
    set_method_name: str
    get_method_name: str

    # I assume that a numpy random generator distribution will be
    # used. many of these just use these properties
    loc: float
    size: float
    distribution_name: str

    # on what to work
    element_index: int

    # find out how to use size properly
    vector_length: int

    def _create_factors(self, method):
        raise AssertionError("Requires to be implemented in a derived class")

    def create_command(self, rng: numpy.random.Generator):
        rng_method = getattr(rng, self.distribution_name)
        add_factors, scale_factors = self._create_factors(rng_method)
        cmd = EngineeringCommand(
            element_index=self.element_index,
            scale_factors=scale_factors,
            add_factors=add_factors,
            set_method_name=self.set_method_name,
            get_method_name=self.get_method_name,
        )
        return cmd


class ScaleEngineeringDistributionCommand(EngineeringDistributionCommand):
    """Create a command for scaling a property
    """

    def _create_factors(self, method):
        add_factors = 0.0
        scale_factors = method(size=self.vector_length)
        return add_factors, scale_factors


class AddEngineeringDistributionCommand(EngineeringDistributionCommand):
    """Create a command for scaling a property
    """

    def _create_factors(self, method):
        add_factors = 0.0
        scale_factors = method(size=self.vector_length)
        return add_factors, scale_factors


def apply_factors(
    acc: tslib.Accelerator, commands: Sequence[EngineeringCommand], copy: bool = True
) -> tslib.Accelerator:
    """Applies factors found in the dataset to the different objects of the accelerator


    Please note the factors here are static ... user is expected to
    generate a different factors dataset for any single setting: i.e if random
    studies are made look to function XXX that will support you in
    generating these factors.

    Furthermore user can apply it at again or chain calls with different files

    Todo:
       * How to apply it to each individual multipole?
       * How to apply it to the different offsets and angle ?
       *  Or how to apply it to scalars and vectors ?

    """

    if copy:
        acc = acc.copy()

    # Iterate over the dataset and apply
    for cmd in commands:
        element = cmd.element_index
        # XXX check get / set method name against a positive list?
        get_method = getattr(element, cmd.get_method_name)
        set_method = getattr(element, cmd.set_method_name)

        a_property = get_method()
        new_property = cmd.scale_factors * a_property + cmd.add_factors
        set_method(new_property)


def create_commands(
    *,
    element_indices: Sequence[int],
    add_factors: np.ndarray,
    scale_factors: np.ndarray,
    get_method_name: str,
    set_method_name: str,
) -> Sequence[EngineeringCommand]:
    """create commands given the available factors

    typically the user will not call it directly but it will be
    rather called by other functions

    Warning:
       just prototype of a function not checked
    """

    n_elem = len(element_indices)

    def check(obj, object_name: str):
        l2 = len(obj)
        if l2 != n_elem:
            txt = (
                f"For object {object_name},"
                f" length of object: len({obj}) == {l2} does not "
                f" match number of element indices {n_elem}"
            )
            raise AssertionError(txt)

    check(add_factors)
    check(scale_factors)

    r = [
        EngineeringCommand(
            element_index=idx,
            add_factors=af,
            scale_factors=sf,
            set_method_name=set_method_name,
            get_method_name=get_method_name,
        )
        for idx, af, sf in zip(element_indices, add_factors, scale_factors)
    ]
    return r


def create_commands_for_randomised_state(
    eng_dist_cmds: Sequence[EngineeringDistributionCommand], rng: numpy.random.Generator
) -> Sequence[EngineeringCommand]:
    """Create commands for current randomised state

    Gets a set of distribution commands: asks each to provide a command of a specific state

    Just iterates over the commands and applies the random number generator to it
    """

    r = [ed.create_command(rng) for ed in eng_dist_cmds]
    return r


def create_distribution_commands(
    element_indices: Sequence[int], reference_command: EngineeringDistributionCommand
) -> Sequence[EngineeringDistributionCommand]:
    """
    """

    def new_command(element_index: int) -> EngineeringDistributionCommand:
        cmd = copy.copy(reference_command)
        cmd.element_index = element_index
        return cmd

    r = [new_command(idx) for idx in element_indices]
    return r


__all__ = [
    "create_distribution_commands",
    "create_commands_for_randomised_state",
    "create_commands",
    "EngineeringDistributionCommand",
]
# def c
#    for idx, af, sf, gm, sm in zip
# def create_random_
