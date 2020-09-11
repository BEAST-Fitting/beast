import numpy as np

__all__ = ["symlog", "inverse_symlog"]


symlog_linthreshold = 1e-20


def symlog(y, C=symlog_linthreshold):
    return np.sign(y) * np.log10(1.0 + np.abs(y / C))


def inverse_symlog(y, C=symlog_linthreshold):
    return np.sign(y) * C * (-1.0 + 10.0 ** np.abs(y))
