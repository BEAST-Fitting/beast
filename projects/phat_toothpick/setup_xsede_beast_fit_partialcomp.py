#!/usr/bin/env python

"""
Code to setup the slurm job files for BEAST trim grid runs on PHAT data

.. history::
    Written 18mar15 by KDG.

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
#####

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("brick", help="brick number")
    args = parser.parse_args()
    brick = args.brick
    
    basename = 'obscat/'

    basepath = 'BEAST_production/'
    cat_files = glob.glob(basepath + 'b'+brick+'/'+basename + 'b*.fits')
    cat_files = np.array(cat_files)

    n_cat_files = len(cat_files)
    n_pernode_files = 500

    # setup the subdirectory for the xsede slurm and log files
    job_path = basepath+'b'+brick+'/fit_xsede_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path+'logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    # get the size of the sed_trim files
    sed_size = np.zeros(n_cat_files)
    for i, cat_file in enumerate(cat_files):
        dpos = string.find(cat_file,'SD-')
        spos = string.find(cat_file,'sub')
        ppos = string.rfind(cat_file,'.')
        sd_num = cat_file[dpos+3:spos-1]
        sub_num = cat_file[spos+3:ppos]
        sed_file = basepath+'b'+brick+'/b'+brick+'_sd'+sd_num+'_sub'+sub_num+'_sed_trim.grid.hd5'
        sed_size[i] = os.path.getsize(sed_file)

    # now sort on size to put similar sized files together for efficiency
    sindxs = np.argsort(sed_size)
    sindxs = sindxs[::-1]
    cat_files = cat_files[sindxs]
    sed_size = sed_size[sindxs]

    cur_f = 0
    cur_total_size = 0.0
    j = -1
    for i, cat_file in enumerate(cat_files):
        # get the sd number
        dpos = string.find(cat_file,'SD-')
        spos = string.find(cat_file,'sub')
        ppos = string.rfind(cat_file,'.')
        sd_num = cat_file[dpos+3:spos-1]
        sub_num = cat_file[spos+3:ppos]

        # read the stats file and see if this subregion is done yet
        results_path = basepath+'b'+brick+'/'
        stats_file = results_path+'b'+brick+'_sd'+sd_num+'_sub'+sub_num+'_stats.fits'
        t = Table.read(stats_file)
        indxs, = np.where(t['Pmax'] != 0.0)

        if len(indxs) == len(t['Pmax']):
            print(stats_file + ' done')
        else:
            print(stats_file + ' not done - adding to continue fitting list (' + str(len(indxs)) + '/' + str(len(t['Pmax'])) + ')')

            j += 1
            if j%n_pernode_files == 0:
                cur_f += 1

                # close previous files
                if j != 0:
                    pf.close()
                    print('total sed_trim size [Gb] = ', cur_total_size/(1024.*1024.*1024.))
                    cur_total_size = 0.0

                # open the slurm and param files
                sf = open(job_path+'beast_xsede_refit_'+str(cur_f)+'.slurm','w')
                joblist_file = job_path+'beast_xsede_refit_'+str(cur_f)+'.joblist'
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
                sf.write('#SBATCH -J BeastF'+str(cur_f)+'             # Job name\n')
                sf.write('#SBATCH -N 1                   # Total number of nodes (16 cores/node)\n')
                sf.write('#SBATCH -n 32                  # Total number of tasks\n')
                sf.write('#SBATCH -p largemem            # Queue name\n')
                sf.write('#SBATCH -o '+log_path+'BeastF'+str(cur_f)+'.o%j         # Name of stdout output file (%j expands to jobid)\n')
                sf.write('#SBATCH -t 48:00:00            # Run time (hh:mm:ss)\n')
                sf.write('#SBATCH --mail-user=kgordon@stsci.edu\n')
                sf.write('#SBATCH --mail-type=begin  # email me when the job starts\n')
                sf.write('#SBATCH --mail-type=end    # email me when the job finishes\n')
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

            pf.write('./run_production_memory.py -f -r ' + brick + ' ' + sd_num + ' '+sub_num+' > '+log_path+'beast_fit_resume_b'+brick+'_sd'+sd_num+'_sub'+sub_num+'.log\n')
            cur_total_size += sed_size[i]

    print('total sed_trim size [Gb] = ', cur_total_size/(1024.*1024.*1024.))
    pf.close()
