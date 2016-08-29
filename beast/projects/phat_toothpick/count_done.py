#!/usr/bin/env python

"""
Code to setup count the number of stars in each brick that are done

.. history::
   Written 15 Apr 2015

"""

from __future__ import print_function
import os
import glob
import string
import argparse
#####
import numpy as np
#import datamodel_production as datamodel
from astropy.table import Table
#####


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--faint", help="Faint set of 4 band detections (instead of bright)",
                        action="store_true")
    parser.add_argument("bricks", metavar='N', type=str, nargs='+',
                        help='bricks to use')
    #parser.add_argument("brick", help="brick number")
    args = parser.parse_args()
    bricks = args.bricks

    basename = 'obscat/'
    basepath = 'BEAST_production/'

    if args.faint:
        ext_brick = 'f'
    else:
        ext_brick = ''
        
    tot_all_poss = 0
    tot_all_done = 0

    for k, brick in enumerate(bricks):
        cat_files = glob.glob(basepath + 'b'+brick+ext_brick+'/'+basename + 'b*.fits')

        tot_poss = 0
        tot_done = 0
        for i, cat_file in enumerate(cat_files):
            bpos = string.find(cat_file,'obscat/')
            dpos = string.find(cat_file,'SD-')
            spos = string.find(cat_file,'sub')
            ppos = string.rfind(cat_file,'.')
            brick_num = cat_file[bpos+8:bpos+10]
            sd_num = cat_file[dpos+3:spos-1]
            sub_num = cat_file[spos+3:ppos]

            # get the number of stars in the obscat
            obs = Table.read(cat_file)
            tot_poss += len(obs)

            # read the stats file and see if this subregion is done yet
            results_path = basepath+'b'+brick_num+ext_brick+'/'
            stats_file = results_path+'b'+brick_num+'_sd'+sd_num+'_sub'+sub_num+'_stats.fits'

            if os.path.isfile(stats_file):
                t = Table.read(stats_file)
                indxs, = np.where(t['Pmax'] != 0.0)
                #print(brick_num, sd_num, sub_num, len(indxs), len(t['Pmax']))
                tot_done += len(indxs)
                
        tot_all_poss += tot_poss
        tot_all_done += tot_done
        print('brick ',brick, tot_poss, tot_done, float(tot_done)/tot_poss)
        
    print('totals: ', tot_all_poss, tot_all_done, float(tot_all_done)/tot_all_poss)
