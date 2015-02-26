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
    
    brick = '19'

    basename = 'obscat/'

    basepath = 'BEAST_production/'
    cat_files = glob.glob(basepath + 'b'+brick+'/'+basename + 'b'+brick+'*.fits')

    # setup the subdirectory for the beast job files
    job_path = basepath+'b'+brick+'/condor_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path+'logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    # setup the subdirectory for the beast batc files
    bjob_path = basepath+'b'+brick+'/batch_jobs/'
    if not os.path.isdir(bjob_path):
        os.mkdir(bjob_path)

    blog_path = bjob_path+'logs/'
    if not os.path.isdir(blog_path):
        os.mkdir(blog_path)

    # open the shell file to population with condor_submit commands
    cs_filename = job_path+'condor_submit_b'+brick
    cs = open(cs_filename, 'w')
    cs.write('#!/bin/bash\n')
    cs.write('#\n')
    cs.write('\n')

    # open the shell file to population with condor_submit commands
    tcs_filename = job_path+'condor_submit_b'+brick+'_trim'
    tcs = open(tcs_filename, 'w')
    tcs.write('#!/bin/bash\n')
    tcs.write('#\n')
    tcs.write('\n')

    # open a txt fle for batch trimming of models
    bt_filename = basepath+'b'+brick + '/b' + brick + '_many_trim.txt' 
    bt = open(bt_filename,'w')
    bt.write('BEAST_production/BEAST_production_seds.grid.hd5\n')
   
    for cat_file in cat_files:
        # get the sd number
        dpos = string.find(cat_file,'SD-')
        spos = string.find(cat_file,'sub')
        ppos = string.rfind(cat_file,'.')
        sd_num = cat_file[dpos+3:spos-1]
        sub_num = cat_file[spos+3:ppos]

        ##########
        # print the min/max observed flux in the file

        # read in the observed data
        #obsdata = datamodel.get_obscat(cat_file, 10.0, datamodel.filters)
        #filtername = 'f475w_rate'
        #print(np.amin(obsdata.data[:][filtername]),np.amax(obsdata.data[:][filtername]))
        ##########

        # write the beast job files
        jobbase = 'beast_phat_toothpick_b'+brick+'_sd'+sd_num+'_sub'+sub_num
        print(jobbase)

        # trim job
        f = open(job_path+jobbase+'_trim.job', 'w')
        f.write('Universe = Vanilla\n')
        f.write('Executable = ./run_production.py\n')
        f.write('Arguments = -t ' + brick + ' ' + sd_num + ' ' + sub_num+'\n')
        f.write('GetEnv = True\n')
        f.write('InitialDir = /user/kgordon/BEAST/beast/projects/phat_toothpick/\n')
        f.write('Log = '+log_path+jobbase+'_trim.log\n')
        f.write('Error = '+log_path+jobbase+'_trim_error.log\n')
        f.write('Output = '+log_path+jobbase+'_trim.out\n')
        f.write('Queue\n')

        f.close()

        # fit job
        f = open(job_path+jobbase+'.job', 'w')
        f.write('Universe = Vanilla\n')
        f.write('Executable = ./run_production.py\n')
        f.write('Arguments = -f ' + brick + ' ' + sd_num + ' ' + sub_num+'\n')
        f.write('GetEnv = True\n')
        f.write('InitialDir = /user/kgordon/BEAST/beast/projects/phat_toothpick/\n')
        f.write('Log = '+log_path+jobbase+'.log\n')
        f.write('Error = '+log_path+jobbase+'_error.log\n')
        f.write('Output = '+log_path+jobbase+'.out\n')
        f.write('Queue\n')

        f.close()

        # fit batch job
        f = open(bjob_path+jobbase+'.batch', 'w')
        f.write('nice -n 10 ./run_production.py -f ' + brick + ' ' + sd_num + ' ' + sub_num+' &> ' + blog_path+jobbase + '.out\n')
        f.close()

        tcs.write('condor_submit ' + job_path+jobbase+'_trim.job\n')
        cs.write('condor_submit ' + job_path+jobbase+'.job\n')

        bt.write(brick + ' ' + sd_num + ' ' + sub_num + '\n')
        
    tcs.write('\n')
    tcs.write('exit 0\n')
    tcs.close()

    cs.write('\n')
    cs.write('exit 0\n')
    cs.close()

    bt.close()
    
    os.chmod(tcs_filename, 0744)
    os.chmod(cs_filename, 0744)

