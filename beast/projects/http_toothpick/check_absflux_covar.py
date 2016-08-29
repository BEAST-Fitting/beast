
import numpy as np
from beast.core.noisemodel import absflux_covmat

filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F555W',
           'HST_ACS_WFC_F775W', 'HST_WFC3_F775W',
           'HST_WFC3_F110W', 'HST_WFC3_F160W']

filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
           'HST_ACS_WFC_F814W',
           'HST_WFC3_F110W', 'HST_WFC3_F160W']

a = absflux_covmat.hst_frac_matrix(filters)

            
print(100.*np.sqrt(a))
