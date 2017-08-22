#!/usr/bin/env python
"""
Creates a source density map based on an input catalog
Used to split the observed catalog and ASTs into source density bins 
   for the BEAST

Based on code by Heddy Arab (2015)

Modified by Karl Gordon Sep 2015
   to be more general (e.g. can be used for other surveys like HTTP)
Updated to be even more general (KDG - Aug 2017)
"""
from __future__ import print_function

import argparse

import numpy as np
import matplotlib.pyplot as plt
import math
from astropy import wcs
from astropy.io import fits
from astropy.table import Table

def make_source_dens_map(catfile, 
                         pix_size = 10., 
                         mag_name = 'F475W_VEGA',
                         mag_cut = [24.5,27]):
    """
    Computes the source density map and store it in a pyfits HDU
    Also writes a text file storing the source density for each source

    INPUTS:
    -------
    catfile: filename of observed catalog
    pix_size: float
         size of pixels in which the source density is computed
    mag_name: string
         name of magnitude column in table
    mag_cut: 2-element list
         magnitude range on which the source density is computed

    OUTPUT:
    -------
    hdu: astropy header data unit
         source density map

    written by HA 01/21/15
    modified by KDG 09/10/15
    """

    cat = Table.read(catfile)

    # Setting map fame
    min_ra = cat['RA'].min()
    max_ra = cat['RA'].max()
    min_dec = cat['DEC'].min()
    max_dec = cat['DEC'].max()
    
    print(min_ra, max_ra, min_dec, max_dec)

    #Compute number of pixel alog each axis pix_size in arcsec
    dec_delt = pix_size/3600.
    n_y = np.fix(np.round((max_dec - min_dec)/dec_delt))
    ra_delt = dec_delt
    n_x = np.fix(np.round(math.cos(0.5*(max_dec+min_dec)*math.pi/180.)
                          *(max_ra-min_ra)/ra_delt))
    ra_delt *= -1.

    n_x = int(np.max([n_x,1]))
    n_y = int(np.max([n_y,1]))

    print('# of x & y pixels = ', n_x, n_y)

    ra_limits = min_ra + ra_delt*np.arange(0,n_x+1, dtype=float)
    dec_limits = min_dec + dec_delt*np.arange(0,n_y+1, dtype=float)
    
    cdelt = [ra_delt, dec_delt]
    crpix = np.asarray([n_x, n_y], dtype = float) / 2.
    crval = np.asarray([(min_ra + max_ra), (min_dec+max_dec)]) / 2.

    w = wcs.WCS(naxis = 2)
    w.wcs.crpix = crpix
    w.wcs.cdelt = cdelt
    w.wcs.crval = crval
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    N_stars = len(cat)

    world = np.zeros((N_stars,2),float)
    world[:,0] = cat['RA']
    world[:,1] = cat['DEC']
    print('working on converting ra, dec to pix x,y')
    pixcrd = w.wcs_world2pix(world, 1)
    pix_x = pixcrd[:,0]
    pix_y = pixcrd[:,1]

    npts_map = np.empty([n_x,n_y])
    source_dens = np.empty(N_stars)

    for i in range(n_x):
        print('x = %s out of %s' % (str(i+1), str(n_x)))
        for j in range(n_y):
            indxs,= np.where((pix_x > i) & (pix_x <= i+1) 
                             & (pix_y > j) & (pix_y <= j+1))
            n_indxs = len(indxs)
            indxs_for_SD,= np.where((cat[mag_name][indxs] >= mag_cut[0]) 
                                    & (cat[mag_name][indxs] <= mag_cut[1]))
            n_indxs = len(indxs_for_SD)
            if n_indxs > 0:
                npts_map[i,j] = n_indxs
            else:
                npts_map[i,j] = 0
            source_dens[indxs] = npts_map[i,j]/(pix_size**2)

    header = w.to_header()
    hdu = fits.PrimaryHDU(npts_map.T, header=header)

    # Save to FITS file
    hdu.writeto(catfile.replace('.fits','_source_den_image.fits'),
                overwrite=True)

    # Save the source density for individual stars in a new catalog file
    cat['SourceDensity'] = source_dens
    cat.write(catfile.replace('.fits','_with_source_density.fits'),
              overwrite=True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("catfile", type=str,
                        help='catalog FITS file')
    parser.add_argument("--pixsize", type=float, default=10.,
                        help='pixel size')
    args = parser.parse_args()
    catfile = args.catfile
    print(catfile)

    make_source_dens_map(args.catfile, pix_size=args.pixsize)
