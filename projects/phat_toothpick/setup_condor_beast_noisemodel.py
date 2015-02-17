#!/usr/bin/env python

"""
Code to setup the condor job files for BEAST runs on PHAT data

.. history::
    Written 6Feb15 by KDG.

"""

from __future__ import print_function
import os
import glob
import string
#####
import numpy as np
import datamodel_production as datamodel
#####

if __name__ == '__main__':
    
    basename = 'merged_asts/'

    basepath = 'BEAST_production/'
    nm_files = np.array(glob.glob(basepath + basename + '*.fits'))

    # setup the subdirectory for the beast job files
    job_path = basepath+'nm_condor_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = basepath+'nm_condor_jobs/logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    # open the shell file to population with condor_submit commands for the noisemodel creation
    nm_filename = job_path+'condor_submit_nm'
    nm = open(nm_filename, 'w')
    nm.write('#!/bin/bash\n')
    nm.write('#\n')
    nm.write('\n')

    # get the 2nd value for the sd region
    sd_num = np.zeros(len(nm_files), dtype=np.int8)
    for k, nm_file in enumerate(nm_files):
        # get the sd number
        dpos = string.rfind(nm_file,'_')
        ppos = string.rfind(nm_file,'.')
        sd_num[k] = int(nm_file[dpos+1:ppos])
        #print(sd_num[k])

    # now sort
    sindxs = np.argsort(sd_num)
    nm_files = nm_files[sindxs]

    for nm_file in nm_files:
        # get the sd number
        dpos = string.find(nm_file,'SD_')
        sd_num = nm_file[dpos+3:len(nm_file)-5]
        print(sd_num)

        # write the noise model condor job files
        jobbase = 'beast_phat_toothpick_sd'+sd_num+'_nm'
        f = open(job_path+jobbase+'.job', 'w')
        f.write('Universe = Vanilla\n')
        f.write('Executable = ./run_production.py\n')
        f.write('Arguments = -n 0 ' + sd_num + ' 0\n')
        f.write('GetEnv = True\n')
        #f.write('Requirements  = Memory >= 22528\n')
        f.write('InitialDir = /user/kgordon/BEAST/beast/projects/phat_toothpick/\n')
        f.write('Log = '+log_path+jobbase+'.log\n')
        f.write('Error = '+log_path+jobbase+'_error.log\n')
        f.write('Output = '+log_path+jobbase+'.out\n')
        f.write('Queue\n')
            
        f.close()
        
        nm.write('condor_submit ' + job_path+jobbase+'.job\n')

    nm.write('\n')
    nm.write('exit 0\n')
    nm.close()
    
    os.chmod(nm_filename, 0744)
