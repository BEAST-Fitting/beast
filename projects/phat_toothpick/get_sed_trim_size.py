#!/usr/bin/env python

"""
Code to save the sed trim file sizes for later use

.. history::
    Written 3jun15 by KDG.

"""

from __future__ import print_function
import os
import glob
import string
import argparse
#####
import numpy as np
import datamodel_production as datamodel
from astropy.table import Table
from astropy.io import fits
#####

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--faint", help="Faint set of 4 band detections (instead of bright)",
                        action="store_true")
    parser.add_argument("-t", "--tag", help="tag for files")
    parser.add_argument("bricks", metavar='N', type=str, nargs='+',
                        help='bricks to use')
    args = parser.parse_args()
    bricks = args.bricks

    print('doing brisk = ', bricks)
    
    if args.faint:
        ext_brick = 'f'
    else:
        ext_brick = ''

    if args.tag:
        ext_tag = '_' + args.tag
    else:
        ext_tag = ''

    basename = 'obscat/'

    basepath = 'BEAST_production/'
    for k, brick in enumerate(bricks):
        icat_files = glob.glob(basepath + 'b'+brick+ext_brick+'/'+basename + 'b*.fits')
        if k == 0:
            cat_files = icat_files
        else:
            cat_files = np.concatenate((cat_files,icat_files))
    cat_files = np.array(cat_files)

    n_cat_files = len(cat_files)
    n_pernode_files = 250

    # get the size of the sed_trim files
    sed_size = np.zeros(n_cat_files)
    for i, cat_file in enumerate(cat_files):
        bpos = string.find(cat_file,'obscat/')
        dpos = string.find(cat_file,'SD-')
        spos = string.find(cat_file,'sub')
        ppos = string.rfind(cat_file,'.')
        brick_num = cat_file[bpos+8:bpos+10]
        sd_num = cat_file[dpos+3:spos-1]
        sub_num = cat_file[spos+3:ppos]

        sed_file = basepath+'b'+brick_num+ext_brick+'/b'+brick_num+'_sd'+sd_num+'_sub'+sub_num+'_sed_trim.grid.hd5'
        if not os.path.isfile(sed_file):
            print('no sed_trim file for ', cat_file)
            sed_size[i] = 0.0
        else:
            sed_size[i] = os.path.getsize(sed_file)

    # now remove all the missing sed_trim subregions
    gindxs, = np.where(sed_size > 0.0)
    cat_files = cat_files[gindxs]
    sed_size = sed_size[gindxs]

    # save the sed trim sizes to an ASCII file
    job_path = basepath+'/refit_xsede_jobs/'
    sedtrim_file = job_path+'beast_xsede_refit_sed_trim_size.dat'
    stf = open(sedtrim_file,'w')

    stf.write('cat_file     size\n')
    for i, cat_file in enumerate(cat_files):
        stf.write(cat_file + "  " + "%.2f" % sed_size[i] + '\n')

    stf.close()
