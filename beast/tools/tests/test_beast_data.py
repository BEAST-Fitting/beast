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
def test_read_lnp_data():
    lnp_fname = download_rename("beast_example_phat_lnp.hd5")
    ldata = read_lnp_data(lnp_fname)

    exp_keys = ["vals", "indxs"]
    for ckey in ldata.keys():
        assert ckey in exp_keys, f"{ckey} not in lnp data expected keys"

    # check an entry for a single model (caching current values 20 Apr 2020)
    # fmt: off
    exp_vals = [-56.83604431, -76.34762573, -17.55770874, -18.23323059, -10.53744507]
    exp_indxs = [14639., 15015., 296., 12636., 1336.]
    # fmt: on
    np.testing.assert_allclose(
        ldata["vals"][0][0:5],
        exp_vals,
        err_msg="Expected posterior (vals) values not correct",
    )
    np.testing.assert_allclose(
        ldata["indxs"][0][0:5], exp_indxs, err_msg="Expected index values not correct",
    )


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


@remote_data
def test_read_sed_data():
    seds_trim_fname = download_rename("beast_example_phat_seds_trim.grid.hd5")
    requested_params = ["Av", "Rv", "f_A", "M_ini", "logA", "Z", "distance"]

    # check that when return_params=True, then just a list of parameters is returned
    sparams = read_sed_data(seds_trim_fname, return_params=True)
    assert isinstance(sparams, list), "Returned params are not a list"
    checknames = requested_params + ["seds", "lamb"]
    for cname in checknames:
        assert cname in sparams, f"{cname} not in sed parameter list"

    # check that otherwise, the requested sed data is returned
    sdata = read_sed_data(seds_trim_fname, param_list=requested_params)
    expected_values = {
        "Av": 0.0,
        "Rv": 2.0,
        "f_A": 1.0,
        "M_ini": 4.0073261261,
        "logA": 6.0,
        "Z": 0.008,
        "distance": 783429.642766212,
    }
    for cname in requested_params:
        assert cname in sdata.keys(), f"requsted parameter {cname} not in sed data"
        np.testing.assert_allclose(
            sdata[cname][10],
            expected_values[cname],
            err_msg=f"expected value of {cname} is not found",
        )


@remote_data
def test_get_lnp_grid_vals():
    lnp_fname = download_rename("beast_example_phat_lnp.hd5")
    ldata = read_lnp_data(lnp_fname)

    seds_trim_fname = download_rename("beast_example_phat_seds_trim.grid.hd5")
    requested_params = ["Av", "Rv", "f_A", "M_ini", "logA", "Z", "distance"]
    sdata = read_sed_data(seds_trim_fname, param_list=requested_params)

    lgvals_data = get_lnp_grid_vals(sdata, ldata)

    # check that otherwise, the requested lgvals data is returned
    expected_values = {
        "Av": [0.0, 0.0, 0.0, 0.0, 0.0],
        "Rv": [2.0, 2.0, 2.0, 2.0, 2.0],
        "f_A": [1.0, 1.0, 1.0, 1.0, 1.0],
        "M_ini": [3.89416909, 3.92726111, 3.95603228, 2.04966068, 2.04999995],
        "logA": [6.0, 6.0, 6.0, 9.0, 9.0],
        "Z": [0.03, 0.03, 0.03, 0.004, 0.004],
        "distance": [
            783429.64276621,
            783429.64276621,
            783429.64276621,
            783429.64276621,
            783429.64276621,
        ],
    }
    for cname in requested_params:
        assert (
            cname in lgvals_data.keys()
        ), f"requsted parameter {cname} not in sed data"
        np.testing.assert_allclose(
            lgvals_data[cname][0:5, 10],
            expected_values[cname],
            err_msg=f"expected value of {cname} is not found",
        )


if __name__ == "__main__":
    test_get_lnp_grid_vals()
