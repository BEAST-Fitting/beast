"""
splinter noise model assumes that every photometric band is independent
from the others and has a fractional flux uncertainty and no bias.

Method
------
Create a noise model that has sigmas that are frac_unc times sed_flux and
zeros for the bias terms.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

import tables


__all__ = ['make_splinter_noise_model']


def make_splinter_noise_model(outname, sedgrid, frac_unc=0.10,
                              absflux_a_matrix=None, **kwargs):
    """
    Splinter noise model assumes that every filter is independent with
    any other.  And assumes a fractional uncertainty at all fluxes.
    No ASTs are used.

    Parameters
    ----------
    outname: str
        path and filename into which save the noise model

    sedgrid: SEDGrid instance
        sed model grid for everyone of which we will evaluate the model

    frac_unc: float [default = 0.10 (10%)]
        fractional flux uncertainy

    absflux_a_matrix: ndarray
        absolute calibration a matrix giving the fractional uncertainties
        including correlated terms (off diagonals)
        for the splinter model, only the diagonal terms are used

    returns
    -------

    noisefile: str
        noisemodel file name
    """

    n_models, n_filters = sedgrid.seds.shape

    # fill the bias vector with zeros
    bias = np.full((n_models, n_filters), 0.0)

    # fill the sigma vector with uncertainties based on the
    #   input fraction uncertainty
    sigma = sedgrid.seds[:]*frac_unc

    # fill the completeness vector with ones
    compl = np.full((n_models, n_filters), 1.0)

    # absolute flux calibration uncertainties
    #  off-diagnonal terms are ignored for the splinter
    if absflux_a_matrix is not None:
        if absflux_a_matrix.ndim == 1:
            abs_calib_2 = absflux_a_matrix[:] ** 2
        else:   # assumes a cov matrix
            abs_calib_2 = np.diag(absflux_a_matrix)

        noise = np.sqrt(abs_calib_2 * sedgrid.seds[:] ** 2 + sigma ** 2)
    else:
        noise = sigma

    # save to disk
    print('Writting to disk into {0:s}'.format(outname))
    with tables.open_file(outname, 'w') as outfile:
        outfile.create_array(outfile.root, 'bias', bias)
        outfile.create_array(outfile.root, 'error', noise)
        outfile.create_array(outfile.root, 'completeness', compl)

    return outname
