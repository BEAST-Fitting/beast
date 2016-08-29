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
#####

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--faint", help="Faint set of 4 band detections (instead of bright)",
                        action="store_true")
    #parser.add_argument("bricks", help="comma list of bricks (e.g, 2,3,4) [no spaces]")
    parser.add_argument("bricks", metavar='N', type=str, nargs='+',
                        help='bricks to use')
    args = parser.parse_args()
    bricks = args.bricks

    basename = 'obscat/'

    if args.faint:
        ext_brick = 'f'
    else:
        ext_brick = ''

    n_bricks = len(bricks)

    basepath = 'BEAST_production/'
    for k, brick in enumerate(bricks):
        icat_files = glob.glob(basepath + 'b'+brick+ext_brick+'/'+basename + 'b*.fits')
        if k == 0:
            cat_files = np.array(icat_files)
        else:
            cat_files = np.concatenate((cat_files,icat_files))

    # find the files that need new trim files
    need_trim = np.empty(len(cat_files),dtype=int)
    for i, cat_file in enumerate(cat_files):
        # get the sd number
        bpos = string.find(cat_file,'obscat/')
        dpos = string.find(cat_file,'SD-')
        spos = string.find(cat_file,'sub')
        ppos = string.rfind(cat_file,'.')
        brick_num = cat_file[bpos+8:bpos+10]
        sd_num = cat_file[dpos+3:spos-1]
        sub_num = cat_file[spos+3:ppos]
        sed_file = basepath+'b'+brick_num+ext_brick+'/b'+brick_num+'_sd'+sd_num+'_sub'+sub_num+'_sed_trim.grid.hd5'
        noise_file = basepath+'b'+brick_num+ext_brick+'/b'+brick_num+'_sd'+sd_num+'_sub'+sub_num+'_noisemodel_trim.hd5'
        if (not os.path.isfile(sed_file)) or (not os.path.isfile(noise_file)):
            need_trim[i] = 1
            if os.path.isfile(sed_file):
                os.remove(sed_file)
            if os.path.isfile(noise_file):
                os.remove(noise_file)
        else:
            need_trim[i] = 0

    gindxs, = np.where(need_trim == 1)
    cat_files = cat_files[gindxs]

    n_cat_files = len(cat_files)
    n_subtrim_files = 32
    n_per_subtrim = int(n_cat_files/32) + 1

    tot_hours = int(n_per_subtrim/4.)

    print('# trim files per process = ',n_per_subtrim)
    print('total hours = ',tot_hours)

    # setup the subdirectory for the xsede slurm and log files
    job_path = basepath+'/trim_xsede_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path+'/logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    # open the noise model xsede slurm and param files
    basefile = 'beast_xsede_trim'
    if args.faint:
        basefile += '_faint'
    sf = open(job_path+basefile+'.slurm','w')
    joblist_file = job_path+basefile+'.joblist'
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
    sf.write('#SBATCH -J BeastTR             # Job name\n')
    sf.write('#SBATCH -N 1                   # Total number of nodes (16 cores/node)\n')
    sf.write('#SBATCH -n 32                  # Total number of tasks\n')
    sf.write('#SBATCH -p largemem            # Queue name\n')
    sf.write('#SBATCH -o '+log_path+'BeastTR.o%j         # Name of stdout output file (%j expands to jobid)\n')
    sf.write('#SBATCH -t '+str(tot_hours)+':00:00            # Run time (hh:mm:ss)\n')
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

    bt_f = []
    for i in range(n_subtrim_files):
        trimfile = job_path+'BeastTR_' + str(i+1)
        bt_f.append(open(trimfile,'w'))
        bt_f[-1].write('BEAST_production/BEAST_production_seds.grid.hd5\n')
        if args.faint:
            pf.write('./trim_many_via_obsdata.py -a '+trimfile+' > '+log_path+'beast_trim_tr'+str(i+1)+'.log\n')
        else:
            pf.write('./trim_many_via_obsdata.py '+trimfile+' > '+log_path+'beast_trim_tr'+str(i+1)+'.log\n')
    pf.close()
    
    sd_nums = np.empty(n_cat_files, dtype=int)
    for i, cat_file in enumerate(cat_files):
        # get the sd number
        dpos = string.find(cat_file,'SD-')
        ddpos = string.find(cat_file,'-',dpos+4)
        sd_nums[i] = int(cat_file[dpos+3:ddpos])

    # now sort on sd num
    sindxs = np.argsort(sd_nums)
    cat_files = cat_files[sindxs]

    k = 0
    n_cur = 0
    for i, cat_file in enumerate(cat_files):
        # get the sd number
        bpos = string.find(cat_file,'obscat/')
        dpos = string.find(cat_file,'SD-')
        spos = string.find(cat_file,'sub')
        ppos = string.rfind(cat_file,'.')
        brick_num = cat_file[bpos+8:bpos+10]
        sd_num = cat_file[dpos+3:spos-1]
        sub_num = cat_file[spos+3:ppos]
        
        bt_f[k].write(brick_num + ' ' + sd_num + ' ' + sub_num + '\n')
        n_cur += 1
        if n_cur >= n_per_subtrim:
            n_cur = 0
            k += 1

    for i in range(n_subtrim_files):
        bt_f[i].close()
