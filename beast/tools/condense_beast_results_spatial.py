#!/usr/bin/env python
#
# condense the spatially reordered BEAST results to single
#   stats, pdf1d, and lnp files per spatial pixel
#
# History: 
#  created Oct 2016 by Karl Gordon

import os
import glob
import math

import h5py
from tqdm import trange, tqdm

import argparse
import numpy as np

from astropy.io import fits
from astropy.table import Table, Column, vstack

def condense_stats_files(bname,
                         cur_dir,
                         out_dir):

    # get all the stats files
    stats_files = glob.glob(cur_dir + '*_stats.fits')

    # loop through the stats files, building up the output table
    cats_list = []
    for cur_stat in stats_files:
        cur_cat = Table.read(cur_stat)

        # get the source density and subregion name
        #  in other words, the reordering tag
        #  bpos is the location after the 2nd underscore of pix coords
        #  epos is before the _stats.fits ending string
        bpos = cur_stat.find(bname+'/')+len(bname)+1
        bpos = cur_stat.find('_',bpos)
        bpos = cur_stat.find('_',bpos+1) + 1
        epos = cur_stat.find('_stats')
        reorder_tag = cur_stat[bpos:epos]
            
        # add the reorder tag to each entry in the current catalog
        n_entries = len(cur_cat)
        cur_cat.add_column(Column([reorder_tag] * n_entries,
                                  name='reorder_tag'))

        # append to list
        cats_list.append(cur_cat)
    
    # concatenate all the small catalogs together
    full_cat = vstack(cats_list)
        
    # output the full pixel catalog
    full_cat.write(out_dir+'/' + bname + '_stats.fits', overwrite=True)

    # return the number of sources in the catalog for later use
    return len(full_cat)

def condense_pdf1d_files(bname,
                         cur_dir,
                         out_dir,
                         n_sources):

    # useful default for the "failure" case of all negative values and log spacing requested
    n_bins_default = 50 

    # get all the files
    pdf1d_files = glob.glob(cur_dir + '*_pdf1d.fits')

    # loop over the pdf1d files, accumulating the condensed 2D format 1d pdfs
    init_condensed = False
    cond_pdf1d_vals = []
    cond_pdf1d_name = []
    for cur_pdf1d in pdf1d_files:

        hdulist = fits.open(cur_pdf1d)
        n_qnames = len(hdulist) - 1
        for k in range(n_qnames):
            pdf1d_histo = hdulist[k+1].data
            n_cur_source, n_bins = pdf1d_histo.shape
            n_cur_source -= 1

            # initialize condensed versions
            if not init_condensed:
                pos_source = 0

                if n_bins == 1:
                    n_bins = n_bins_default

                cond_pdf1d_vals.append(np.empty((n_sources+1, n_bins),
                                                dtype=float))
                cond_pdf1d_name.append(hdulist[k+1].header['EXTNAME'])

                # copy in the bin values to the last column
                cond_pdf1d_vals[k][n_sources,:] = pdf1d_histo[-1,:]

            # insert the 1d pdfs into the 2D format structure
            #print(k, len(cond_pdf1d_vals))
            cond_pdf1d_vals[k][pos_source:pos_source+n_cur_source,:] = \
                pdf1d_histo[0:-1,:]
        
        pos_source += n_cur_source
        init_condensed = True

    # output the condensed 2D format 1D pdfs
    hdulist = fits.HDUList([fits.PrimaryHDU()])

    # generate the extensions
    for k, qname in enumerate(cond_pdf1d_name):
        chdu = fits.PrimaryHDU(cond_pdf1d_vals[k])
        chdu.header.set('XTENSION','IMAGE') 
        chdu.header.set('EXTNAME',qname) 
        
        hdulist.append(chdu)

    # write the 1D PDFs
    hdulist.writeto(out_dir+'/'+bname+'_pdf1d.fits', overwrite=True)

def condense_lnp_files(bname,
                       cur_dir,
                       out_dir):

    # get all the files
    lnp_files = glob.glob(cur_dir + '*_lnp.hd5')

    # open the condensed hd5 file for writing
    #   remove it if it already exisits (h5py does not overwrite)
    clfile = out_dir+'/'+bname+'_lnp.hd5'
    try:
        os.remove(clfile)
    except OSError:
        pass

    cond_lnp_file = h5py.File(clfile)

    # loop over the small lnp files and copy to main lnp file
    k = 0
    for cur_lnp in lnp_files:
        cur_lnpfile = h5py.File(cur_lnp, 'r')
        
        # loop over all the stars (groups)
        for sname in cur_lnpfile.keys():
            star_group = cond_lnp_file.create_group('star_%d' % k)
            for cp_name, cp_value in cur_lnpfile[sname].items():
                star_group.create_dataset(cp_name,data=cp_value.value)

            k += 1

        cur_lnpfile.close()

    cond_lnp_file.close()

if __name__ == '__main__':

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--bricknum", 
                        help="PHAT brick num shortcut" + \
                        " (superceeds other input)")
    parser.add_argument("-d","--filedir", 
                        help="Directory to condense results")
    args = parser.parse_args()

    if args.bricknum:
        brick = str(args.bricknum)
        out_dir = '/astro/dust_kg2/kgordon/BEAST_production/b' + \
            brick + '/spatial'
    elif args.filedir:
        out_dir = args.filedir
    else:
        parser.print_help()
        exit()

    if not os.path.exists(out_dir):
        print(out_dir + ' directory does not exist')
        exit()

    # get the list of directories
    #    each directory is a different pixel
    pix_dirs = sorted(glob.glob(out_dir + '/*/'))

    # loop over each subdirectory and condense the files as appropriate
    for cur_dir in tqdm(pix_dirs, desc='spatial regions'):

        # get the base name
        spos = cur_dir.rfind('/',0,len(cur_dir)-1)
        bname = cur_dir[spos+1:-1]

        # process that catalog (stats) files
        n_sources = condense_stats_files(bname, cur_dir, out_dir)

        # process the pdf1d files
        condense_pdf1d_files(bname, cur_dir, out_dir, n_sources)

        # process the nD lnp files
        condense_lnp_files(bname, cur_dir, out_dir)
