from thor_scsi.flame import GLPSParser, Config
from thor_scsi.lib import Accelerator
import os


def parse_config_file(file_name: str, t_dir=None) -> Config:
    """parse an accelerator configuration file

    Args:
        filename: filename of the configuration file
        t_dir:    directory where further configuration is to be found?

    Returns:
        a :class:`thor_scsi.flame.Config` instance representing
        the information of the configuration file in an abstract
        synatx tree

    Todo:
        Find out how t_dir works precicesly

    """
    if t_dir is None:
        # where should be the default directory to look into?
        t_dir = os.environ["HOME"]

    with open(file_name) as fp:
        text = fp.read()

    C = GLPSParser().parse_byte(text, t_dir)
    return C


def accelerator_from_config(filename: str) -> Accelerator:
    '''Create an accelerator object from configuration file

    Reads configuration and returns accelerator object
    '''
    C = parse_config_file(filename)
    acc = Accelerator(C)
    return acc


__all__ = ["parse_config_file", "accelerator_from_config"]
