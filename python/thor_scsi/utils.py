from typing import Sequence
import numpy as np
import numpy.ma as ma


__all__ = ["ss_vect_to_masked_matrix"]


def ss_vect_to_masked_matrix(t_lists: Sequence[Sequence[float]])-> ma.MaskedArray :
    """Convert ss_vect lists to an masked array

    The lists are not of equal length. Thus a matrix is returned
    with the elements that are not set masked out.
    """
    n_rows = len(t_lists)
    n_cols = np.max([len(l) for l in t_lists])

    data = np.zeros([n_rows, n_cols], dtype=np.float)
    data[:, :] = np.nan
    mask = np.ones([n_rows, n_cols], np.bool)
    for i in range(len(t_lists)):
        t_list = t_lists[i]
        n = len(t_list)
        data[i, :n] = t_list
        mask[i, :n] = False

    res = ma.array(data=data, mask=mask)
    return res
