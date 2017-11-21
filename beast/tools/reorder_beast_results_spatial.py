#!/usr/bin/env python
#
# reorder the BEAST results to be in ra/dec bins
#
# optimized BEAST runs are done by sets of stars in source density bins
#   sorted by flux and subdivided into smaller files
#   this allow for the BEAST grid to be cut and speeds up the fitting
# 
# but this is non-ideal for almost any analysis of the results
#   especially the MegaBEAST!
#
# History: based on work by Heddy Arab prior to 2016
#          reordered code and flow to optimize for speed by Karl Gordon (09/16)

import os
import glob
import math

import h5py
from tqdm import trange, tqdm

import argparse
import numpy as np

from astropy import wcs
from astropy.io import fits
from astropy.table import Table

def setup_spatial_regions(cat_filename,
                          pix_size=10.0):
    """
    The spatial regions are setup via a WCS object
        
    Parameters
    ----------
    cat_filename : string
       filename of catalog

    pix_size : float
       size of pixels/regions in arcsec

    Returns
    -------
    wcs_info: astropy WCS object
    """
    
    # read in the catalog file
    cat = Table.read(cat_filename)

    # min/max ra
    min_ra = cat['RA'].min()
    max_ra = cat['RA'].max()
    min_dec = cat['DEC'].min()
    max_dec = cat['DEC'].max()

    # ra/dec delta values
    dec_delt = pix_size/3600.
    ra_delt = dec_delt

    # compute the number of pixels and 
    n_y = int(np.rint((max_dec - min_dec)/dec_delt) + 1)
    n_x = int(np.rint(math.cos(0.5*(max_dec+min_dec)*math.pi/180.)*
                      (max_ra-min_ra)/ra_delt) + 1)

    # ra delta should be negative
    ra_delt *= -1.

    print('# of x & y pixels = ', n_x, n_y)

    w = wcs.WCS(naxis = 2)
    w.wcs.crpix = np.asarray([n_x, n_y], dtype = float) / 2.
    w.wcs.cdelt = [ra_delt, dec_delt]
    w.wcs.crval = np.asarray([(min_ra+max_ra), (min_dec+max_dec)]) / 2.
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    return (w, n_x, n_y)

def regions_for_objects(ra,
                        dec,
                        wcs_info):
    """
    Generate the x,y coordinates for each object based on the input
    ra/dec and already created WCS information.
        
    Parameters
    ----------
    ra : array of float
       right ascension of the objects

    dec : array of float
       declination of the objects

    wcs_info: astropy WCS object
       previously generated WCS object based on the full catalog

    Returns
    -------
    dictonary of:

    x : int array
      x values of regions

    y : int array
      y values of regions

    name : str array
      string array composed of x_y
    """

    # generate the array needed for fast conversion
    world = np.empty((len(ra),2),float)
    world[:,0] = ra
    world[:,1] = dec

    # convert
    pixcrd = wcs_info.wcs_world2pix(world, 1)

    # get the arrays to return
    x = pixcrd[:,0].astype(int)
    y = pixcrd[:,1].astype(int)
    xy_name = [None]*len(ra)

    x_str = x.astype(np.string_)
    y_str = y.astype(np.string_)
                  
    for k in range(len(x)):
        xy_name[k] = str(x[k]) + '_' + str(y[k])
    
    # return the results as a dictonary
    #   values are truncated to provide the ids for the subregions
    return {'x': x, 'y': y, 'name': xy_name}

if __name__ == '__main__':

    # command line params to specify the run directory
    #   and any other needed parameters

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--bricknum", 
                        help="PHAT brick num shortcut" + \
                        " (superceeds other inputs)")
    parser.add_argument("-s","--stats_filename", 
                        help="Filename of the full run stats")
    parser.add_argument("-r","--region_filebase", 
                        help="Filebase of the run regions")
    parser.add_argument("-o","--output_filebase", 
                        help="Filebase to use for output")
    parser.add_argument("-p","--reg_size", default=10., type=float,
                        help="spatial region size [arcsec]")
    args = parser.parse_args()

    if args.bricknum:
        brick = str(args.bricknum)
        cat_filename = '/astro/dust_kg2/harab/toothpick_results/v1_1/b' + \
            brick + '_stats_v1_1.fits'
        reg_filebase = '/astro/dust_kg2/kgordon/BEAST_production/b' + \
            brick + '/b' + brick 
        out_dir = '/astro/dust_kg2/kgordon/BEAST_production/b' + \
            brick + '/spatial'
        out_filebase = out_dir + '/b' + brick
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    elif (args.stats_filename and args.region_filebase):
        cat_filename = args.stats_filename    
        reg_filebase = args.region_filebase
        out_filebase = args.output_filebase
    else:
        parser.print_help()
        exit()

    # size of the regions (square, units in arcsec)
    reg_size = args.reg_size

    # read in the full brick catalog and setup the spatial subdivided regions
    wcs_info, n_x, n_y = setup_spatial_regions(cat_filename,
                                               pix_size=reg_size)

    # setup array to store number of stars per pixel
    wcs_nstars = np.zeros((n_y, n_x), dtype=int)

    # find all the subdivided BEAST files for this brick
    sub_files = glob.glob(reg_filebase + '*_stats.fits')

    # loop over the files and output to the appropriate spatial region files
    # in loop:
    #      read in the locations of each star and calculated spatial region
    #      append the 1D pdfs
    #      append the sparse nD pdfs
    #      append the completeness function (??)

    #sub_files = sub_files[0:2]
    for cur_file in tqdm(sub_files, desc="orig sub files"):
        # read in the stats file
        cur_cat = Table.read(cur_file)
        n_objs = len(cur_cat)

        # read in the pdf1d file
        #   a number of extensions, one for each entry in the stats file
        hdulist = fits.open(cur_file.replace('_stats.fits','_pdf1d.fits'))
        cur_pdf1d_vals = []
        cur_pdf1d_name = []
        n_qnames = len(hdulist) - 1
        for k in range(n_qnames):
            cur_pdf1d_name.append(hdulist[k+1].header['EXTNAME'])
            cur_pdf1d_vals.append(hdulist[k+1].data)
        hdulist.close()

        # open the lnp file for reading
        cur_lnpfile = h5py.File(cur_file.replace('_stats.fits',
                                                 '_lnp.hd5'), 'r')
        
        # get the source density and subregion tag
        # allows for unique filenames for the spatial regions for output
        orig_reg_tag = cur_file[cur_file.find('_sd'):cur_file.find('_stats')]
        
        # determine the subregions for all the objects
        xy_vals = regions_for_objects(cur_cat['RA'],
                                      cur_cat['DEC'],
                                      wcs_info)

        # get the unique xy regions
        xy_names = np.squeeze(xy_vals['name'])
        uniq_xy_names, rindxs = np.unique(xy_names, 
                                          return_inverse=True)

        #uniq_xy_names = uniq_xy_names[0:5]

        # loop over the unique xy regions
        for k in trange(len(uniq_xy_names), desc="outputing " + orig_reg_tag, 
                        leave=False):
            uxy_name = uniq_xy_names[k]

            # get the indexes for the objects in this region
            indxs, = np.where(rindxs == k)

            # add the number of stars found to summary array
            wcs_nstars[xy_vals['y'][indxs[0]], 
                       xy_vals['x'][indxs[0]]] += len(indxs)

            # create region directory if it does not exist
            reg_dir = out_filebase + '_' + uxy_name
            if not os.path.exists(reg_dir):
                os.makedirs(reg_dir)

            # base filename for output
            reg_filebase = out_filebase + '_' + uxy_name + '/' + \
                uxy_name + orig_reg_tag

            # write the stats info
            reg_stats_file = reg_filebase + '_stats.fits'
            cur_cat[indxs].write(reg_stats_file, overwrite=True)

            # write the pdf1d info
            reg_pdf1d_file = reg_filebase + '_pdf1d.fits'

            # setup the primary header and hdulist
            hdulist = fits.HDUList([fits.PrimaryHDU()])

            # generate the extensions
            for k, qname in enumerate(cur_pdf1d_name):
                # get the 1D PDFs for the cur objects 
                #   plus the last column giving the values of the bins
                cur_reg_pdf1d = cur_pdf1d_vals[k][np.append(indxs,n_objs),:]

                chdu = fits.PrimaryHDU(cur_reg_pdf1d)
                chdu.header.set('XTENSION','IMAGE') 
                chdu.header.set('EXTNAME',qname) 

                hdulist.append(chdu)

            # write the 1D PDFs
            hdulist.writeto(reg_pdf1d_file, clobber=True)

            # write the nD sparse likelihood info
            reg_lnp_file = reg_filebase + '_lnp.hd5'

            # open the file (overwrites an existing file)
            reg_lnpfile = h5py.File(reg_lnp_file, 'w')
            
            for i, k in enumerate(indxs):
                star_group = reg_lnpfile.create_group('star_%d' % i)

                # transfer the information
                for cp_name, cp_value in cur_lnpfile['star_%d' % k].items():
                    star_group.create_dataset(cp_name,data=cp_value.value)

            reg_lnpfile.close()

        cur_lnpfile.close()

        # Now, write out the WCS info and number of stars per pixel to file
        #   do every subregion file to have an on-the-fly check
        header = wcs_info.to_header()
        hdu = fits.PrimaryHDU(wcs_nstars, header=header)

        # Save to FITS file
        hdu.writeto(out_filebase+'_nstars.fits', clobber=True)
