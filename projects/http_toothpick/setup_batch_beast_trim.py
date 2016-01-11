#!/usr/bin/env python

"""
Code to setup the batch files for BEAST trim grid runs on HTTP data

.. history::
    Written 23sep15 by KDG.

"""

from __future__ import print_function
import os
import glob
import string
import argparse
#####
import numpy as np
import datamodel_production as datamodel
#####

if __name__ == '__main__':

    basename = 'obscat/'

    basepath = 'HTTP_production/'
    cat_files = np.array(glob.glob(basepath+'/'+basename + 'http*.fits'))

    n_cat_files = len(cat_files)
    n_subtrim_files = 8
    n_per_subtrim = int(n_cat_files/n_subtrim_files) + 1

    print('# trim files per process = ',n_per_subtrim)

    # setup the subdirectory for the xsede slurm and log files
    job_path = basepath+'/trim_batch_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path+'/logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    sd_nums = np.empty(n_cat_files, dtype=int)
    for i, cat_file in enumerate(cat_files):
        # get the sd number
        dpos = string.find(cat_file,'SD_')
        ddpos = string.find(cat_file,'-',dpos+4)
        sd_nums[i] = int(cat_file[dpos+3:ddpos])

    # now sort on sd num
    sindxs = np.argsort(sd_nums)
    cat_files = cat_files[sindxs]

    joblist_file = job_path+'beast_xsede_trim.joblist'
    pf = open(joblist_file,'w')

    bt_f = []
    for i in range(n_subtrim_files):
        trimfile = job_path+'HTTPTR_' + str(i+1)
        bt_f.append(open(trimfile,'w'))
        bt_f[-1].write('HTTP_production/HTTP_production_seds.grid.hd5\n')
        pf.write('./trim_many_via_obsdata.py '+trimfile+' > '+log_path+'http_trim_tr'+str(i+1)+'.log\n')
    pf.close()
    
    k = 0
    n_cur = 0
    for i, cat_file in enumerate(cat_files):
        # get the sd number
        dpos = string.find(cat_file,'SD_')
        spos = string.find(cat_file,'sub')
        ppos = string.rfind(cat_file,'.')
        sd_num = cat_file[dpos+3:spos-1]
        sub_num = cat_file[spos+3:ppos]
        
        bt_f[k].write(sd_num + ' ' + sub_num + '\n')
        n_cur += 1
        if n_cur >= n_per_subtrim:
            n_cur = 0
            k += 1

    for i in range(n_subtrim_files):
        bt_f[i].close()
