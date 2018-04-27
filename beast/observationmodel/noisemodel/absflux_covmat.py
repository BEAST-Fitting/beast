"""
Function to generate fractional absolute flux covariance matrix for
HST and (potentially other) photometric filters
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from astropy.io.fits import getdata

from .. import phot
from ...config import __ROOT__

from ...tools.pbar import Pbar

def hst_frac_matrix(filters, spectrum=None, progress=True,
                    hst_fname=None, filterLib=None):
    """ Uses the Bohlin et al. (2013) provided spectroscopic
    absolute flux covariance matrix to generate the covariance matrix
    for the input set of HST filters.

    Parameters
    ----------
    filters : filter names

    Keywords
    --------
    progress: bool, optional
            if set, display a progress bar
    spectrum : 2 element tuple
               (wave, sed)
               wave = 1D numpy array with wavelengths in XX units
               spectrum = 1D numpy array with flux in ergs...
    hst_fname : str
                file with hst absflux covariance matrix
    filterLib:  str
        full filename to the filter library hd5 file


    Returns
    -------
    2D numpy array giving the fractional covariance matrix
      (must be multiplied by the SED flux (x2) to get the
       true covariance matrix)

    ToDos:
    ------
    - Probably better to do a proper integration than just a weighted
      sum (minor).
    """

    # get the HST fractional covariance matrix at spectroscopic resolution
    if hst_fname is None:
        hst_fname = __ROOT__+'/hst_whitedwarf_frac_covar.fits'

    hst_data = getdata(hst_fname,1)

    waves = hst_data['WAVE'][0]
    frac_spec_covar = hst_data['COVAR'][0]
    n_waves = len(waves)

    # define a flat spectrum if it does not exist
    if spectrum is None:
        spectrum = (waves, np.full((n_waves),1.0))

    # read in the filter response functions
    flist = phot.load_filters(filters, filterLib=filterLib,
                              interp=True, lamb=waves)

    # setup multiplication images to make it easy to compute the results
    n_filters = len(filters)
    mult_image = np.empty((n_waves,n_waves,n_filters))
    mult_image_spec = np.empty((n_waves,n_waves,n_filters))
    image_ones = np.full((n_waves,n_waves),1.0)
    band_ones = np.full((n_filters,n_filters),1.0)

    for i in range(n_filters):
        mult_image[:,:,i] = image_ones*flist[i].transmit

    # handle single spectrum or many spectra
    if len(spectrum[1].shape) > 1:
        n_models = spectrum[1].shape[0]
        results = np.empty((n_models, n_filters, n_filters))
    else:
        n_models = 1
        progress = False

    # setup the progress bar
    if progress is True:
        it = Pbar(desc='Calculating AbsFlux Covariance ' + \
                  'Matrices').iterover(list(range(n_models)))
    else:
        it = list(range(n_models))

    frac_covar_bands = np.empty((n_filters,n_filters))
    for k in it:
        if n_models == 1:
            interp_spectrum = np.interp(waves, spectrum[0], spectrum[1])
        else:
            interp_spectrum = np.interp(waves, spectrum[0], spectrum[1][k,:])


        for i in range(n_filters):
            mult_image_spec[:,:,i] = mult_image[:,:,i]*interp_spectrum

        for i in range(n_filters):
            for j in range(i,n_filters):
                frac_covar_bands[i,j] = np.sum(frac_spec_covar*
                                               mult_image_spec[:,:,i]*
                                               mult_image_spec[:,:,j].T)
                frac_covar_bands[i,j] /= np.sum(mult_image_spec[:,:,i]*
                                                mult_image_spec[:,:,j].T)

        # fill in the symmetric terms
        for i in range(n_filters):
            for j in range(0,i):
                frac_covar_bands[i,j] = frac_covar_bands[j,i]

        # add the term accounting for the uncertainty in the overall
        #  zero point of the flux scale
        #  (e.g., uncertainty in Vega at 5555 A)
        frac_covar_bands += 4.9e-5

        if n_models > 1:
            results[k,:,:] = frac_covar_bands

    if n_models == 1:
        return frac_covar_bands
    else:
        return results
