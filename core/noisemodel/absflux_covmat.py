"""
Function to generate fractional absolute flux covariance matrix for
HST and (potentially other) photometric filters
"""

import numpy as np
from astropy.io.fits import getdata    

from .. import phot
from ...config import __ROOT__

def hst_frac_matrix(filters, spectrum=None):
    """ Uses the Bohlin et al. (2013) provided spectroscopic
    absolute flux covariance matrix to generate the covariance matrix
    for the input set of HST filters.

    Keywords
    ----------
    filters : filter names
    spectrum : 2 element tuple
               (wave, sed)
               wave = 1D numpy array with wavelengths in XX units
               spectrum = 1D numpy array with flux in ergs...

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
    hst_data = getdata(__ROOT__+'/libs/hst_whitedwarf_frac_covar.fits',1)

    waves = hst_data['WAVE'][0]
    frac_spec_covar = hst_data['COVAR'][0]
    n_waves = len(waves)

    # define a flat spectrum if it does not exist
    if spectrum is None:
        spectrum = (waves, np.full((n_waves),1.0))

    # read in the filter response functions
    flist = phot.load_filters(filters, interp=True, lamb=waves)

    # setup multiplication images to make it easy to compute the results
    n_filters = len(filters)
    mult_image = np.empty((n_waves,n_waves,n_filters))
    image_ones = np.full((n_waves,n_waves),1.0)

    # handle single spectrum or many spectra
    if len(spectrum[1].shape) > 1:
        n_models = spectrum[1].shape[0]
        results = np.empty((n_models, n_filters, n_filters))
    else:
        n_models = 1

    print(n_models)

    for k in range(0,n_models,100):
        print(k)
        if n_models == 1:
            interp_spectrum = np.interp(waves, spectrum[0], spectrum[1])
        else:
            interp_spectrum = np.interp(waves, spectrum[0], spectrum[1][k,:])

        frac_covar_bands = np.empty((n_filters,n_filters))
    
        for i in range(n_filters):
            mult_image[:,:,i] = image_ones*flist[i].transmit*interp_spectrum

        for i in range(n_filters):
            for j in range(n_filters):
                frac_covar_bands[i,j] = np.sum(frac_spec_covar*
                                               mult_image[:,:,i]*
                                               mult_image[:,:,j].T)
                frac_covar_bands[i,j] /= np.sum(mult_image[:,:,i]*
                                                mult_image[:,:,j].T)

        frac_covar_bands += 4.9e-5
        print(np.sqrt(frac_covar_bands))

        if n_models > 1:
            results[k,:,:] = frac_covar_bands


    if n_models == 1:
        return frac_covar_bands
    else:
        return results
