"""
Function to generate fractional absolute flux covariance matrix for
HST and (potentially other) photometric filters
"""

import numpy as np
from astropy.io.fits import getdata    

from .. import phot
from ...config import __ROOT__

def hst_frac_matrix(filters):
    """ Uses the Bohlin et al. (2013) provided spectroscopic
    absolute flux covariance matrix to generate the covariance matrix
    for the input set of HST filters.

    Keywords
    ----------
    filters : filter names

    Returns
    -------
    2D numpy array giving the fractional covariance matrix
      (must be multiplied by the SED flux (x2) to get the
       true covariance matrix)

    ToDos:
    ------
    - Probably better to do a proper integration than just a weighted sum (minor).
    - allow for the input of a spectrum (or many spectra) to compute the covariance matrix
      for a specific spectrum (currently assumes a flat spectrum)
      this would allow for the covariance matrix to be computed for each model SED
    """

    # get the HST fractional covariance matrix at spectroscopic resolution
    hst_data = getdata(__ROOT__+'/libs/hst_whitedwarf_frac_covar.fits',1)

    waves = hst_data['WAVE'][0]
    frac_spec_covar = hst_data['COVAR'][0]

    # read in the filter response functions
    flist = phot.load_filters(filters, interp=True, lamb=waves)

    # setup multiplication images to make it easy to compute the results
    n_filters = len(filters)
    n_waves = len(waves)

    mult_image_x = np.empty((n_waves,n_waves,n_filters))
    mult_image_y = np.empty((n_waves,n_waves,n_filters))
    for k in range(n_filters):
        for i in range(n_waves):
            mult_image_x[i,:,k] = flist[k].transmit
            mult_image_y[:,i,k] = flist[k].transmit

    frac_covar_bands = np.empty((n_filters,n_filters))
    
    for i in range(n_filters):
        for j in range(n_filters):
            frac_covar_bands[i,j] = np.sum(frac_spec_covar*mult_image_x[:,:,i]*mult_image_y[:,:,j])/np.sum(mult_image_x[:,:,i]*mult_image_y[:,:,j])

    frac_covar_bands += 4.9e-5

    return frac_covar_bands
