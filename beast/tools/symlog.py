import numpy as np

__all__ = ["symlog", "inverse_symlog"]


symlog_linthreshold = 1e-20


def symlog(y, C=symlog_linthreshold):
    """
    Compute the symlog value that transitions between linear and log scales
    based on the theshold value C.  Symlog works for both positive and
    negative input values.

    Parameters
    ----------
    y : float
        intput value

    C : float
        value of the threshold where the symlog changes from linear to log

    Returns
    -------
    y_symlog : float
        symlog of y
    """
    return np.sign(y) * np.log10(1.0 + np.abs(y / C))


def inverse_symlog(y, C=symlog_linthreshold):
    """
    Compute the linear value based on the input symlog value.

    Parameters
    ----------
    y : float
        intput symlog value

    C : float
        value of the threshold where the symlog changes from linear to log

    Returns
    -------
    y_linear : float
        linear version of y
    """
    return np.sign(y) * C * (-1.0 + 10.0 ** np.abs(y))
