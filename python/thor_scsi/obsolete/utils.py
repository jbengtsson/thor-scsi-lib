from typing import Sequence
import numpy as np
import numpy.ma as ma


__all__ = ["ss_vect_to_masked_matrix"]


def ss_vect_to_masked_matrix(t_lists: Sequence[Sequence[float]])-> ma.MaskedArray:
    """Convert ss_vect lists to an masked array

    The lists are not of equal length. Thus a matrix is returned
    with the elements that are not set masked out.

    Todo:
      Check: is that function not already obsolete?
      See :func:`thor_scsi.`
    """
    # Contains the data of the phase space vector
    ps_tmp = t_lists[0]

    ps = np.array(ps_tmp[:-1], dtype=np.float)

    # First order Taylor expansion: jacobian derivatives
    mat_lists = t_lists[0:-1]
    # this extra element is currently not added at the export function
    # remove it here as soon as the export code is fixed.
    padding = t_lists[-1] + [0]

    tmp = mat_lists + [padding]
    l_cols = [len(l) for l in tmp]
    n_cols = np.max(l_cols)
    # print(f"mat: Number of columns {l_cols}, {n_cols}")
    mat = np.array(tmp, dtype=np.float)
    return ps, mat
