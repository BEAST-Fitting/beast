import numpy as np
from astropy.tests.helper import remote_data

from beast.tools.read_beast_data import (
    read_lnp_data,
    read_noise_data,
    read_sed_data,
    get_lnp_grid_vals,
)
from beast.tests.helpers import download_rename


@remote_data
def test_read_noise_data():
    noise_trim_fname = download_rename("beast_example_phat_noisemodel_trim.grid.hd5")
    ndata = read_noise_data(noise_trim_fname)

    exp_keys = ["bias", "completeness", "error"]
    for ckey in ndata.keys():
        assert ckey in exp_keys, f"{ckey} not in noise data expected keys"

    # check an entry for a single model (caching current values 18 Apr 2020)
    # fmt: off
    exp_bias = [-6.77602149e-20, -1.36353610e-20, 2.87448605e-20,
                -2.38253474e-21, -1.70330281e-20, -2.70390708e-20]
    exp_error = [1.63128160e-19, 7.50503350e-20, 7.65873857e-20,
                 2.48842055e-20, 9.41313147e-20, 2.79650823e-20]
    exp_compl = [1.0, 0.95552407, 1.0, 0.74733078, 0.77777778, 0.42857143]
    # fmt: on
    np.testing.assert_allclose(
        ndata["bias"][10], exp_bias, err_msg="Expected bias values not correct",
    )
    np.testing.assert_allclose(
        ndata["error"][10], exp_error, err_msg="Expected error values not correct",
    )
    np.testing.assert_allclose(
        ndata["completeness"][10],
        exp_compl,
        err_msg="Expected completeness values not correct",
    )


if __name__ == "__main__":
    test_read_noise_data()
