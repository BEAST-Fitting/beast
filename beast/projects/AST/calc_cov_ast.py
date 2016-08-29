"""
Script to extract all the information from Artificial Star Tests
Computes the covariance matrices, their determinants (lnQ) and the inverse of their Cholesky decomposition.
Write the results in a hd5 file. 
Author: HA 
Date: 11/27/13
"""

import numpy as np    
import pyfits
import tables
from beast.core.hdfstore import *
from beast.proba.helpers import calc_covar_mat

# Reading ASTs file
dir_data = '/home/arab/beast_data/' # path to change
a = pyfits.getdata(dir_data+'12059_M31-B17-F03.gst.fake.fits')

max_asts = len(a)
N_bands = 6

bands = ['F275W','F336W','F475W','F814W','F110W','F160W']

# Extract fluxes and pack into variables
imags = np.zeros((N_bands, max_asts)) # input magnitudes
rflux = np.zeros((N_bands,max_asts))  # recovered fluxes
rmags = np.zeros((N_bands,max_asts))  # recovered magnitudes
ave_rflux=np.zeros((N_bands))         # average recovered fluxes
unc_rflux=np.zeros((N_bands))         # uncertainty recovered flluxes

for k in np.arange(6):
    imags[k,:] = a[bands[k]+'_IN']
    rflux[k,:] = a[bands[k]+'_RATE']
    rmags[k,:] = a[bands[k]+'_VEGA']

iflux = np.power(10, -0.4*imags) # Conversion to flux unit/vega_flux

# Find the stars
N_stars = len(np.unique(imags[5,:],return_index=True)[0]) # Number of artificial stars
uindxs = np.unique(imags[5,:],return_index=True)[1]       # Indices of first realization of a given star

# Output variable definitions
cov_mat = np.zeros((N_stars,N_bands, N_bands))
inv_cholesky = np.zeros((N_stars,N_bands, N_bands))
lnQ =  np.zeros((N_stars))
outname="AST_results2.hd5" # Output file name

with tables.openFile(outname, 'w') as outfile:
    for i in np.arange(N_stars):
        tuple_indxs = np.where(imags[5,:] == imags[5,uindxs[i]])
        indxs=tuple_indxs[0]
        N_indxs = len(indxs) # Number of realizations
        
        # Find non-detections (At least 1 detection)
        gtindxs = np.ones(N_indxs)
        for j in np.arange(N_indxs):
            l = indxs[j]
            tuple_tindxs = np.where(rmags[:,l] > 90)
            N_tindxs = len(tuple_tindxs[0])
            if (N_tindxs > 5): 
                gtindxs[j]=0
                
        tuple_nindxs = np.where(gtindxs == 1)
        nindxs = tuple_nindxs[0]
        N_indxs = len(nindxs)
        # At least 20 realizations
        if (N_indxs > 0): 
            indxs = indxs[nindxs]
        if (N_indxs > 20):
            # Compute the covariance matrices, the inverse cholesky decomp. and lnQ
            (cov_mat[i,:,:], inv_cholesky[i,:,:], lnQ[i]) = calc_covar_mat(iflux[:,indxs].T,rflux[:,indxs].T)
            # Compute the mean of the recovered fluxes + uncertainties
            for k in np.arange(N_bands):
                ave_rflux[k] = np.mean(rflux[k,indxs])
                unc_rflux[k] =np.std(rflux[k,indxs])
        
        # Write the outputs
        star_group = outfile.createGroup('/', 'star_%d' % i, title="star %d" % i)
        outfile.createArray(star_group, 'input_flux', iflux[:,indxs[0]])
        outfile.createArray(star_group, 'recov_flux', ave_rflux[:])
        outfile.createArray(star_group, 'unc_flux', unc_rflux[:])
        outfile.createArray(star_group, 'cov_mat', cov_mat[i,:,:])
        outfile.createArray(star_group, 'inv_cholesky', inv_cholesky[i,:,:])
        outfile.createArray(star_group, 'lnQ', [lnQ[i]])
        outfile.flush()
outfile.close()

