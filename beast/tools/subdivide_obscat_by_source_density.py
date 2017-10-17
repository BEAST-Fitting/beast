#!/usr/bin/env python
#
# Split the observed catalog by source density.
# Needed as the noise model is dependent on source density (e.g. crowding).
#
# Based on code by Heddy Arab (2015)
#
# Modified by Karl Gordon Sep 2015
#    to be more general (e.g. can be used for other surveys like HTTP)

import argparse

from astropy.table import Table
import astropy.io.fits as pyfits
import numpy as np
import os

def split_obs_by_source_density(catfile, bin_width=1, sort_col='F475W_RATE', 
                                Ns_file = 6250):
    """
    Splits the observation in different source density files

    !!! Note: This function uses the updated catalog file containing the
    source density for each star in the brick 
    (check for create_source_density_map.py) !!!

    INPUTS:
    -------
    catfile: filename of the observed catalog
    bin_width: float
        Width of source density bin in star/arcsec
    sort_col: column with which to sort the split files
    Ns_file: integer
        Number of sources per subfile
    
    OUTPUT: None
    -------
    """

    obs = Table.read(catfile)
     
    SD_cut = obs['SourceDensity']
    SD_max = int(np.max(SD_cut)) + 1
    SD_min = int(np.min(SD_cut))
    print('Source density min/max = ', SD_min, SD_max)

    bins = np.arange(SD_min, SD_max, bin_width)  # source density binning
    inds = np.digitize(SD_cut,bins)              # finding in which source density bin each star fall
    inds = inds-1                                # np.digitize output indexes start at 1 
    
    uinds = np.unique(inds)           # Selecting the populated bins
    for ek, ok in enumerate(uinds):
        ii, = np.where(inds == ok)
        cc = SD_cut[ii]
        print(ok,ok+1, cc.min(), cc.max())
        assert((cc.min() >= ok) & (cc.max() < ok+1))   # checking we are in the right bin
        aa = np.str(bins[ok])                          # lower value in bin (string naming purpose)
        bb = np.str(bins[ok] + bin_width)              # higher value in bin (string naming purpose)
        print('SD ' + aa + '-' + bb + ' # sources = ' + str(len(ii)))
        sdobs = obs[ii]
        sdfile = catfile.replace('.fits','_SD_' + aa + '-' + bb + '.fits')
        sdobs.write(sdfile,overwrite=True)

        sindxs = np.argsort(sdobs[sort_col])# Sort the stars 

        N = len(sindxs)
        Nb_files = int(N/Ns_file) + 1   # Computing the number of subfiles
        print('dividing into ' + str(Nb_files) + ' subfiles for later fitting speed')

        # Writing the files
        for i in range(Nb_files):
            min_k = i*Ns_file
            if i < Nb_files:
                max_k = (i+1)*Ns_file
            else:
                max_k = N

            sdobs[sindxs[min_k:max_k]].write(sdfile.replace('.fits','_sub' + str(i) + '.fits'),overwrite=True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("catfile", type=str,
                        help='catalog FITS file') 
    parser.add_argument("--n_per_file", type=int, default=6250,
                        help='Number of sources per subfile')
    parser.add_argument("--sort_col", type=str, default='F475W_RATE',
                        help='Column in catalog file for sorting')
    args = parser.parse_args()

    split_obs_by_source_density(args.catfile, 
                                sort_col=args.sort_col,
                                Ns_file=args.n_per_file)
