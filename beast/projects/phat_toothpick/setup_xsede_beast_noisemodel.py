#!/usr/bin/env python

"""
Code to setup the slurm job files for BEAST noisemodel runs on PHAT data

.. history::
    Written 18mar15 by KDG.

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

    # setup the subdirectory for the xsede slurm and log files
    job_path = basepath+'nm_xsede_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = basepath+'nm_xsede_jobs/logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    # open the noise model xsede slurm and param files
    sf = open(job_path+'beast_xsede_noisemodel.slurm','w')
    joblist_file = job_path+'beast_xsede_noisemodel.joblist'
    pf = open(joblist_file,'w')

    # fill the slurm file with all the needed lines (lots of them)
    sf.write('#!/bin/csh\n')
    sf.write('#\n')
    sf.write('# Simple SLURM script for submitting multiple serial\n')
    sf.write('# jobs (e.g. parametric studies) using a script wrapper\n')
    sf.write('# to launch the jobs.\n')
    sf.write('#\n')
    sf.write('# To use, build the launcher executable and your\n')
    sf.write('# serial application(s) and place them in your WORKDIR\n')
    sf.write('# directory.  Then, edit the CONTROL_FILE to specify \n')
    sf.write('# each executable per process.\n')
    sf.write('#-------------------------------------------------------\n')
    sf.write('#-------------------------------------------------------\n')
    sf.write('# \n')
    sf.write('#         <------ Setup Parameters ------>\n')
    sf.write('#\n')
    sf.write('#SBATCH -J BeastNM             # Job name\n')
    sf.write('#SBATCH -N 1                   # Total number of nodes (16 cores/node)\n')
    sf.write('#SBATCH -n 32                  # Total number of tasks\n')
    sf.write('#SBATCH -p largemem            # Queue name\n')
    sf.write('#SBATCH -o '+log_path+'BeastNM.o%j         # Name of stdout output file (%j expands to jobid)\n')
    sf.write('#SBATCH -t 01:00:00            # Run time (hh:mm:ss)\n')
    sf.write('#      <------------ Account String ------------>\n')
    sf.write('# <--- (Use this ONLY if you have MULTIPLE accounts) --->\n')
    sf.write('##SBATCH -A \n')
    sf.write('#------------------------------------------------------\n')
    sf.write('\n')
    sf.write('module load launcher\n')
    sf.write('setenv EXECUTABLE     $TACC_LAUNCHER_DIR/init_launcher \n')
    sf.write('setenv CONTROL_FILE   '+joblist_file+'\n')
    sf.write('setenv WORKDIR        .\n')
    sf.write('# \n')
    sf.write('# Variable description:\n')
    sf.write('#\n')
    sf.write('#  EXECUTABLE     = full path to the job launcher executable\n')
    sf.write('#  CONTROL_FILE   = text input file which specifies\n')
    sf.write('#                   executable for each process\n')
    sf.write('#                   (should be located in WORKDIR)\n')
    sf.write('#  WORKDIR        = location of working directory\n')
    sf.write('#\n')
    sf.write('#      <------ End Setup Parameters ------>\n')
    sf.write('#--------------------------------------------------------\n')
    sf.write('#--------------------------------------------------------\n')
    sf.write('\n')
    sf.write('#----------------\n')
    sf.write('# Error Checking\n')
    sf.write('#----------------\n')
    sf.write('\n')
    sf.write('if ( ! -e $WORKDIR ) then\n')
    sf.write('echo " "\n')
    sf.write('echo "Error: unable to change to working directory."\n')
    sf.write('echo "       $WORKDIR"\n')
    sf.write('echo " "\n')
    sf.write('echo "Job not submitted."\n')
    sf.write('exit\n')
    sf.write('endif\n')
    sf.write('\n')
    sf.write('if ( ! -f $EXECUTABLE ) then\n')
    sf.write('echo " "\n')
    sf.write('echo "Error: unable to find launcher executable $EXECUTABLE."\n')
    sf.write('echo " "\n')
    sf.write('echo "Job not submitted."\n')
    sf.write('exit\n')
    sf.write('endif\n')
    sf.write('\n')
    sf.write('if ( ! -f $WORKDIR/$CONTROL_FILE ) then\n')
    sf.write('echo " "\n')
    sf.write('echo "Error: unable to find input control file $CONTROL_FILE."\n')
    sf.write('echo " "\n')
    sf.write('echo "Job not submitted."\n')
    sf.write('exit\n')
    sf.write('endif\n')
    sf.write('\n')
    sf.write('\n')
    sf.write('#----------------\n')
    sf.write('# Job Submission\n')
    sf.write('#----------------\n')
    sf.write('\n')
    sf.write('cd $WORKDIR/\n')
    sf.write('echo " WORKING DIR:   $WORKDIR/"\n')
    sf.write('\n')
    sf.write('$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE\n')
    sf.write('\n')
    sf.write('echo " "\n')
    sf.write('echo " Parameteric Job Complete"\n')
    sf.write('echo " "\n')
    sf.close()

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

        pf.write('./run_production.py -n 0 ' + sd_num + ' 0 > '+log_path+'beast_nm_sd'+sd_num+'.log\n')

    pf.close()
