#!/usr/bin/env python
"""
Creates a source density map based on an input catalog
Used to split the observed catalog and ASTs into source density bins
   for the BEAST
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
    FITS files written to disk.
    """

    cat = Table.read(catfile)
    # force catalog column names to be upper case
    for name in cat.colnames:
         cat.rename_column(name, name.upper())
    # force filter magnitude name to be upper case to match column names
    mag_name = mag_name.upper()

    # get the columns with fluxes
    rate_cols = [s for s in cat.colnames if s[-4:] == 'RATE']
    n_filters = len(rate_cols)

    # create the indexs where any of the rates are zero and non-zero
    #   zero = missing data, etc. -> bad for fitting
    #   non-zero = good data, etc. -> great for fitting
    initialize_zero = False
    band_zero_indxs = {}
    print('band, good, zero')
    for cur_rate in rate_cols:
        cur_good_indxs, = np.where(cat[cur_rate] != 0.0)
        cur_indxs, = np.where(cat[cur_rate] == 0.0)
        print(cur_rate, len(cur_good_indxs), len(cur_indxs))
        if not initialize_zero:
            initialize_zero = True
            zero_indxs = cur_indxs
            nonzero_indxs = cur_good_indxs
        else:
            zero_indxs = np.union1d(zero_indxs, cur_indxs)
            nonzero_indxs = np.intersect1d(nonzero_indxs, cur_good_indxs)

        # save the zero indexs for each band
        band_zero_indxs[cur_rate] = zero_indxs

    print('all bands', len(nonzero_indxs), len(zero_indxs))

    # Setting map fame
    min_ra = cat['RA'].min()
    max_ra = cat['RA'].max()
    min_dec = cat['DEC'].min()
    max_dec = cat['DEC'].max()

    #Compute number of pixel alog each axis pix_size in arcsec
    dec_delt = pix_size/3600.
    n_y = np.fix(np.round((max_dec - min_dec)/dec_delt))
    ra_delt = dec_delt
    n_x = np.fix(np.round(math.cos(0.5*(max_dec+min_dec)*math.pi/180.)
                          *(max_ra-min_ra)/ra_delt))
    #ra_delt *= -1. #Not sure why the ra delta would want to be negative...

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

    npts_map = np.zeros([n_x,n_y], dtype=float)
    npts_zero_map = np.zeros([n_x,n_y], dtype=float)
    npts_band_zero_map = np.zeros([n_x,n_y,n_filters], dtype=float)
    source_dens = np.zeros(N_stars, dtype=float)

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
                npts_map[i,j] = n_indxs/(pix_size**2)

                # now make a map of the sources with zero fluxes in
                #   at least one band
                zindxs, = np.where((pix_x[zero_indxs] > i)
                                   & (pix_x[zero_indxs] <= i+1)
                                   & (pix_y[zero_indxs] > j)
                                   & (pix_y[zero_indxs] <= j+1))
                if len(zindxs) > 0:
                    npts_zero_map[i,j] = len(zindxs)

                # do the same for each band
                for k, cur_rate in enumerate(rate_cols):
                    tindxs = band_zero_indxs[cur_rate]
                    zindxs, = np.where((pix_x[tindxs] > i)
                                       & (pix_x[tindxs] <= i+1)
                                       & (pix_y[tindxs] > j)
                                       & (pix_y[tindxs] <= j+1))
                    if len(zindxs) > 0:
                        npts_band_zero_map[i,j,k] = len(zindxs)

            # save the source density as an entry for each source
            source_dens[indxs] = npts_map[i,j]

    # Save to FITS file
    header = w.to_header()
    hdu = fits.PrimaryHDU(npts_map.T, header=header)
    hdu.writeto(catfile.replace('.fits',
                                '_source_den_image.fits'),
                overwrite=True)

    # Save to FITS file (zero flux sources)
    header = w.to_header()
    hdu = fits.PrimaryHDU(npts_zero_map.T, header=header)
    hdu.writeto(catfile.replace('.fits','_npts_zero_fluxes_image.fits'),
                overwrite=True)

    for k, cur_rate in enumerate(rate_cols):
        # Save to FITS file (zero flux sources)
        header = w.to_header()
        hdu = fits.PrimaryHDU(npts_band_zero_map[:,:,k].T, header=header)
        hdu.writeto(catfile.replace('.fits','_npts_zero_fluxes_'
                                    + cur_rate + '_image.fits'),
                    overwrite=True)

    # Save the source density for individual stars in a new catalog file
    cat['SourceDensity'] = source_dens
    cat.write(catfile.replace('.fits','_with_sourceden_inc_zerofluxes.fits'),
              overwrite=True)

    # Save the source density for individual stars in a new catalog file
    #   only those that have non-zero fluxes in all bands
    cat[nonzero_indxs].write(catfile.replace('.fits',
                                             '_with_sourceden.fits'),
                             overwrite=True)
    
    bin_details = Table(names=['i_ra', 'i_dec', 'sourcedens', 'min_ra', 'max_ra', 'min_dec', 'max_dec'])

    for i in range(n_x):

        for j in range(n_y):

            bin_details.add_row([i, j, npts_map[i, j], ra_limits[i], ra_limits[i + 1], dec_limits[j], dec_limits[j + 1]])

    bin_details.write(catfile.replace('.fits', '_source_den_bins_table.fits'))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("catfile", type=str,
                        help='catalog FITS file')
    parser.add_argument("--pixsize", type=float, default=10.,
                        help='pixel size')
    args = parser.parse_args()

    make_source_dens_map(args.catfile, pix_size=args.pixsize)
