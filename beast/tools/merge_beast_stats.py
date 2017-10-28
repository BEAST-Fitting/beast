#!/usr/bin/env python
#
# condense separate BEAST stats files into a single stats file

import os
import glob
import math

import argparse
import numpy as np

from astropy.io import fits
from astropy.table import Table, Column, vstack


def merge_stats_files(stats_files,
                      out_stats_filebase):

    # loop through the stats files, building up the output table
    cats_list = []
    for cur_stat in stats_files:
        cur_cat = Table.read(cur_stat)

        # get the source density and subregion name
        #  in other words, the reordering tag
        #  bpos is the location after the 2nd underscore of pix coords
        #  epos is before the _stats.fits ending string
        bpos = cur_stat.find('_sd') + 1
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
    full_cat.write(out_stats_filebase + '_stats.fits', overwrite=True)

    # return the number of sources in the catalog for later use
    return len(full_cat)

if __name__ == '__main__':

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filebase", 
                        help="filebase to use (e.g., xxx_*_stats.fits)")
    args = parser.parse_args()

    # get the files to merge
    stats_files = glob.glob(args.filebase + '*_stats.fits')

    # do the merge
    n_objs = merge_stats_files(stats_files, args.filebase)

