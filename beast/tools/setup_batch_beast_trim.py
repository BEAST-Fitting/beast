#!/usr/bin/env python

"""
Code to setup the batch files for BEAST trim grid runs 
"""

from __future__ import print_function
import os
import glob

import argparse

import numpy as np


def setup_batch_beast_trim(project,
                               datafile,
                               astfile,
                               num_subtrim=5,
                               nice=None):
    """
    Sets up batch files for submission to the 'at' queue on linux (or similar) systems


    Parameters
    ----------
    project : string
        project name to use (basename for files)

    datafile : string
        file with the observed data (FITS file) - the observed photometry, not the sub-files

    astfile : string
        file with ASTs

    num_subtrim : int (default = 5)
        number of trim batch jobs

    nice : int (default = None)
        set this to an integer (-20 to 20) to prepend a "nice" level to the trimming command
  
    """

    ast_file = astfile
    n_subtrim_files = num_subtrim

    
    full_model_filename = "%s/%s_seds.grid.hd5"%(project,project)

    cat_files = np.array(glob.glob(datafile.replace('.fits','*_sub*.fits')))

    n_cat_files = len(cat_files)

    # make sure n_subtrim_files isn't larger than the number of catalog sub-files
    if n_subtrim_files > n_cat_files:
        n_subtrim_files = n_cat_files

    # max number of files per process
    n_per_subtrim = int(n_cat_files/n_subtrim_files)
    if n_cat_files % n_subtrim_files != 0:
        n_per_subtrim += 1


    print('# trim files per process = ',n_per_subtrim)

    # setup the subdirectory for the batch and log files
    job_path = project+'/trim_batch_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path+'logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    sd_nums = np.empty(n_cat_files, dtype=int)
    for i, cat_file in enumerate(cat_files):
        # get the sd number
        dpos = cat_file.find('SD_')
        ddpos = cat_file.find('-',dpos+4)
        sd_nums[i] = int(cat_file[dpos+3:ddpos])

    # now sort on sd num
    sindxs = np.argsort(sd_nums)
    cat_files = cat_files[sindxs]

    # prepend nice
    nice_str = ''
    if nice is not None:
        nice_str = 'nice -n' + str(int(nice)) + ' '


    joblist_file = job_path+'beast_batch_trim.joblist'
    pf = open(joblist_file,'w')

    bt_f = []
    for i in range(n_subtrim_files):
        trimfile = job_path+'BEAST_' + str(i+1)
        bt_f.append(open(trimfile,'w'))
        bt_f[-1].write(project + '\n')
        #bt_f[-1].write(full_model_filename + '\n')
        pf.write(nice_str + 'python -m beast.tools.trim_many_via_obsdata '+trimfile+' > '
                 +log_path+'beast_trim_tr'+str(i+1)+'.log\n')
    pf.close()
    
    k = 0
    n_cur = 0
    for i, cat_file in enumerate(cat_files):
        # get the sd number
        dpos = cat_file.find('SD_')
        spos = cat_file.find('sub')
        ppos = cat_file.rfind('.')
        sd_num = cat_file[dpos+3:spos-1]
        sub_num = cat_file[spos+3:ppos]

        bt_f[k].write(sd_num + ' ' + sub_num)
        bt_f[k].write(' ' + cat_file)
        bt_f[k].write(' ' + ast_file)
        bt_f[k].write('\n')
        n_cur += 1
        if n_cur >= n_per_subtrim:
            n_cur = 0
            k += 1

    for i in range(n_subtrim_files):
        bt_f[i].close()



if __name__ == '__main__':

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("projectname",
                        help="project name to use (basename for files)")
    parser.add_argument("datafile",
                        help="file with the observed data (FITS file)")
    parser.add_argument("astfile",
                        help="file with ASTs")
    parser.add_argument("--num_subtrim", default=5, type=int,
                        help="number of trim batch jobs")
    args = parser.parse_args()

    project = args.projectname
    datafile = args.datafile
    ast_file = args.astfile
    n_subtrim_files = args.num_subtrim

    setup_batch_beast_trim(project, datafile, ast_file, n_subtrim_files=5)
