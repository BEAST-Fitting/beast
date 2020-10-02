import numpy as np

from beast.tools.symlog import symlog, inverse_symlog


def test_symlog():

    y = np.array([-1e-10, -1e-19, -1e-20, -1e-24, 0.0, 1e-24, 1e-20, 1e-19, 1e-10])

    exp_symlog = [
        -1.00000000e01,
        -1.04139269e00,
        -3.01029996e-01,
        -4.34272769e-05,
        0.00000000e00,
        4.34272769e-05,
        3.01029996e-01,
        1.04139269e00,
        1.00000000e01,
    ]

    # test going to symlog
    np.allclose(symlog(y), exp_symlog)

    # test converting back to linear
    np.allclose(inverse_symlog(exp_symlog), y)
