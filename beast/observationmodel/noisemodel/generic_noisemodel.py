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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import tables

from . import toothpick

__all__ = ['Generic_ToothPick_Noisemodel','make_toothpick_noise_model',
           'get_noisemodelcat']

class Generic_ToothPick_Noisemodel(toothpick.MultiFilterASTs):

    def set_data_mappings(self):
        """
        hard code mapping directly with the interface to ASTs format
        """
        for k in self.filters:
            try:
                self.data.set_alias(k + '_out',
                                    k.split('_')[-1].upper() + '_VEGA')
                #self.data.set_alias(k + '_rate',
                #                    k.split('_')[-1].upper() + '_RATE')
                self.data.set_alias(k + '_in',
                                    k.split('_')[-1].upper() + '_IN')
            except Exception as e:
                print(e)
                print('Warning: Mapping failed. This could lead to ' +
                      'wrong results')


def make_toothpick_noise_model(outname, astfile, sedgrid,
                               use_rate=False, vega_fname=None,
                               absflux_a_matrix=None, **kwargs):
    """ toothpick noise model assumes that every filter is independent with
    any other.

    Parameters
    ----------
    outname: str
        path and filename into which save the noise model

    astfile: str
        path to the file into which are ASTs results

    sedgrid: SEDGrid instance
        sed model grid for everyone of which we will evaluate the model

    use_rate: boolean
        set to use the rate column (normalized vega flux) 
        instead of out column (mags)

    absflux_a_matrix: ndarray
        absolute calibration a matrix giving the fractional uncertainties
        including correlated terms (off diagonals)

    returns
    -------
    noisefile: str
        noisemodel file name
    """
    
    # read in AST results
    model = Generic_ToothPick_Noisemodel(astfile, sedgrid.filters,
                                         vega_fname=vega_fname)

    # compute binned biases and uncertainties as a function of flux
    if use_rate:
        # change the mappings for the out column to the rate column
        for cfilt in sedgrid.filters:
            model.data.set_alias(cfilt + '_out',
                                 cfilt.split('_')[-1].upper() + '_RATE')
        model.fit_bins(nbins=30, completeness_mag_cut=-10)
    else:
        model.fit_bins(nbins=30, completeness_mag_cut=80)

    #for k in range(len(model.filters)):
    #    print(model.filters[k])
    #    print(model._fluxes[:,k])
    #    print(model._sigmas[:,k]/model._fluxes[:,k])
    #    print(model._biases[:,k]/model._fluxes[:,k])
    #    print(model._compls[:,k])

    # evaluate the noise model for all the models in sedgrid
    bias, sigma, compl = model(sedgrid)

    # absolute flux calibration uncertainties
    #  currently we are ignoring the off-diagnonal terms
    if absflux_a_matrix is not None:
        if absflux_a_matrix.ndim == 1:
            abs_calib_2 = absflux_a_matrix[:] ** 2
        else:   # assumes a cov matrix
            abs_calib_2 = np.diag(absflux_a_matrix)

        noise = np.sqrt(abs_calib_2 * sedgrid.seds[:] ** 2 + sigma ** 2)
    else:
        noise = sigma

    # check if the noise model has been extrapolated at the faint flux levels
    # if so, then set the noise to a negative value (later may be used to
    # trim the model of "invalid" models)
    # we are assuming that extrapolation at high fluxes is ok as the noise
    # will be very small there
    for k in range(len(model.filters)):
        indxs, = np.where(sedgrid.seds[:,k] <= model._minmax_asts[0,k])
        if len(indxs) > 0:
            noise[indxs,k] *= -1.0

    print('Writting to disk into {0:s}'.format(outname))
    with tables.open_file(outname, 'w') as outfile:
        outfile.create_array(outfile.root,'bias', bias)
        outfile.create_array(outfile.root,'error', noise)
        outfile.create_array(outfile.root,'completeness', compl)

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
    table: pytables.Table
        table containing the elements of the noise model
    """
    return tables.open_file(filename)


if __name__ == '__main__':

    pass
