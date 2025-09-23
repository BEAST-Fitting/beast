import numpy as np

__all__ = ["compute_grid_weights", "compute_bin_boundaries"]


def compute_grid_weights(in_x, log=False):
    """
    Compute the grid weights.  Needed for marginalization (aka integration).  The
    weights are the relative widths of of each x bin.

    Parameters
    ----------
    x : numpy array
       centers of each bin

    log : boolean
       set if values are in log units

    Returns
    -------
    weights : numpy array
        weights as bin widths divided by the average width
    """
    # ensure x values are monotonically increasing
    sindxs = np.argsort(in_x)
    x = in_x[sindxs]

    n_x = len(x)
    bin_hdiffs = np.diff(x) / 2.0

    # define the bin min and max boundaries
    # handling the two edge cases
    bin_mins = np.zeros(n_x)
    bin_mins[1:] = x[1:] - bin_hdiffs
    bin_mins[0] = x[0] - bin_hdiffs[0]

    bin_maxs = np.zeros(n_x)
    bin_maxs[0:-1] = x[0:-1] + bin_hdiffs
    bin_maxs[-1] = x[-1] + bin_hdiffs[-1]

    if log:
        weights = (10**bin_maxs) - (10**bin_mins)
    else:
        weights = bin_maxs - bin_mins

    # put the weights in the same order as in_x
    out_weights = np.zeros(n_x)
    out_weights[sindxs] = weights

    # return normalized weights to avoid numerical issues
    return out_weights / np.average(out_weights)


def compute_bin_boundaries(tab, noneg=False):
    """
    Computes the boundaries of bins

    The bin boundaries are defined as the midpoint between each value in tab.
    At the two edges, 1/2 of the bin width is subtracted/added to the
    min/max of tab.

    Parameters
    ----------
    tab : numpy array
       centers of each bin

    Returns
    -------
    tab2 : numpy array
       boundaries of the bins
    """
    temp = tab[1:] - np.diff(tab) / 2.0
    tab2 = np.zeros(len(tab) + 1)
    tab2[0] = tab[0] - np.diff(tab)[0] / 2.0
    if noneg & (tab2[0] < 0.0):
        tab2[0] = 0.5 * tab[0]
    tab2[-1] = tab[-1] + np.diff(tab)[-1] / 2.0
    tab2[1:-1] = temp
    return tab2
