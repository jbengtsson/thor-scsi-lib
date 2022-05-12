"""Some math tools

Not found yet somewhere else. If found please replace
"""
import numpy as np


def distances_for_sequence(a: np.ndarray) -> np.ndarray:
    """Compute matrix of distances, spoofed trace

    Trace is set to the largest distance. The returned
    distance matrix is intended to be used for tests
    """
    a = np.atleast_1d(a)
    # Ensure only one dimension
    n_elem, = a.shape

    d = a[:, np.newaxis] - a[np.newaxis, :]
    d = np.absolute(d)

    # Ignore trace ... lions (not only in Cairo)
    d_max = d.max()
    idx = np.arange(n_elem)
    d[idx, idx] = d_max

    return d


def minimum_distance_above_threshold(a: np.ndarray, threshold: float):
    """Minimum distance between elements above threshold

    Returns:
         flag, index1, index2

    flag is True if above threshold (then index1 and index2 are None)

    Otherwise the indices of the corresponding elements are returned
    (typically the closest ones)
    """
    d = distances_for_sequence(a)
    t_min = d.min()
    if t_min > threshold:
        # As expected
        return True, None, None

    idx = d.argmin()
    n_elem = len(a)
    n_col = idx % n_elem
    n_row = idx // n_elem

    return False, n_col, n_row
