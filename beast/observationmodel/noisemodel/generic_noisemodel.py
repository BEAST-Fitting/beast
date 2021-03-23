"""
Generates a generic noise model from artifical star tests (ASTs) results
using the toothpick method.  Using ASTs results in a noise model that
includes contributions from measurement (photon) noise *and* crowding
noise.

Toothpick assumes that all bands are independent - no covariance.
This is a conservative assumption.  If there is true covariance
more accurate results with smaller uncertainties on fit parameters
can be achieved using the trunchen method.  The trunchen method
requires significantly more complicated ASTs and many more of them.
"""
import numpy as np
import h5py
import tables

from beast.observationmodel.noisemodel import toothpick

__all__ = [
    "make_toothpick_noise_model",
    "get_noisemodelcat",
]


def make_toothpick_noise_model(
    outname,
    astfile,
    sedgrid,
    vega_fname=None,
    absflux_a_matrix=None,
    nfluxbins=50,
    **kwargs,
):
    """ toothpick noise model assumes that every filter is independent with
    any other.

    Parameters
    ----------
    outname : str
        path and filename into which save the noise model

    astfile : str
        path to the file into which are ASTs results

    sedgrid : SEDGrid instance
        sed model grid for everyone of which we will evaluate the model

    absflux_a_matrix : ndarray
        absolute calibration a matrix giving the fractional uncertainties
        including correlated terms (off diagonals)

    nfluxbins : int (default=50)
        number of flux bins

    returns
    -------
    noisefile: str
        noisemodel file name
    """

    # read in AST results
    model = toothpick.MultiFilterASTs(astfile, sedgrid.filters, vega_fname=vega_fname)

    # set the column mappings as the external file is BAND_VEGA or BAND_IN
    model.set_data_mappings(in_pair=("in", "in"), out_pair=("out", "rate"), upcase=True)

    # compute binned biases and uncertainties as a function of flux
    model.fit_bins(nbins=nfluxbins)

    # evaluate the noise model for all the models in sedgrid
    bias, sigma, compl = model(sedgrid)

    # absolute flux calibration uncertainties
    #  currently we are ignoring the off-diagnonal terms
    if absflux_a_matrix is not None:
        if absflux_a_matrix.ndim == 1:
            abs_calib_2 = absflux_a_matrix[:] ** 2
        else:  # assumes a cov matrix
            abs_calib_2 = np.diag(absflux_a_matrix)

        noise = np.sqrt(abs_calib_2 * sedgrid.seds[:] ** 2 + sigma ** 2)

        # check if the noise model has been extrapolated at the faint or bright flux levels
        # if so, then set the noise to a negative value (later may be used to
        # trim the model of "invalid" models)
        # if the noise model has been extrapolated, the completeness is set to zeros
        for k in range(len(model.filters)):
            (indxs,) = np.where(compl[:, k] <= 0.0)
            if len(indxs) > 0:
                noise[indxs, k] *= -1.0

    else:
        noise = sigma

    print("Writing to disk into {0:s}".format(outname))
    with tables.open_file(outname, "w") as outfile:
        outfile.create_array(outfile.root, "bias", bias)
        outfile.create_array(outfile.root, "error", noise)
        outfile.create_array(outfile.root, "completeness", compl)

    return outname


def get_noisemodelcat(filename):
    """
    returns the noise model

    Parameters
    ----------
    filename: str
        file containing the outputs from OneD_ASTs_ModelGenerator

    Returns
    -------
    ntable : dict
        dictonary containing the elements of the noise model
    """
    nfile = h5py.File(filename, "r")

    # create a dictonary of the elements
    ntable = {}
    for ckey in nfile.keys():
        ntable[ckey] = np.array(nfile[ckey])

    nfile.close()

    # check that at least the 3 basic elements are included
    expected_elements = ["error", "bias", "completeness"]
    for cexp in expected_elements:
        if cexp not in ntable.keys():
            raise ValueError(f"{cexp} values not found in noisemodel")

    return ntable
