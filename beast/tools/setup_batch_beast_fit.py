#!/usr/bin/env python

"""
Code to create the batch files for fitting
"""

from __future__ import print_function
import os
import glob
import argparse
import tables
#####
import numpy as np
from astropy.table import Table
from astropy.io import fits
#####


def setup_batch_beast_fit(projectname,
                              datafile,
                              num_percore=5,
                              nice=None,
                              overwrite_logfile=True):
    """
    Sets up batch files for submission to the 'at' queue on linux (or similar) systems


    Parameters
    ----------
    project : string
        project name to use (basename for files)

    datafile : string
        file with the observed data (FITS file) - the observed photometry, not the sub-files

    num_percore : int (default = 5)
        number of fitting runs per core

    nice : int (default = None)
        set this to an integer (-20 to 20) to prepend a "nice" level to the fitting command

    overwrite_logfile : boolean (default = True)
        if True, will overwrite the log file; if False, will append to existing log file
  
    """
    

    project = projectname

    cat_files = np.array(glob.glob(datafile.replace('.fits','*_sub*.fits')))

    n_cat_files = len(cat_files)
    n_pernode_files = num_percore

    # setup the subdirectory for the batch and log files
    job_path = project+'/fit_batch_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path+'logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    pf_open = False
    cur_f = 0
    cur_total_size = 0.0
    j = -1

    #cat_files = cat_files[0:2]

    for i, cat_file in enumerate(cat_files):
        # get the sd number
        dpos = cat_file.find('SD_')
        spos = cat_file.find('sub')
        ppos = cat_file.rfind('.')
        sd_num = cat_file[dpos+3:spos-1]
        sub_num = cat_file[spos+3:ppos]

        # read the stats file and see if this subregion is done yet
        results_path = project + '/'
        basename = "%s/%s_sd%s_sub%s"%(project,project,sd_num,sub_num)

        stats_file = basename + '_stats.fits'
        pdf1d_file = basename + '_pdf1d.fits'
        lnp_file = basename + '_lnp.hd5'

        reg_run = False
        run_done = False
        if not os.path.isfile(stats_file):
            reg_run = True
            print('no stats file')
        if not os.path.isfile(pdf1d_file):
            reg_run = True
            print('no pdf1d file')
        if not os.path.isfile(lnp_file):
            reg_run = True
            print('no lnp file')

        # first check if the pdf1d mass spacing is correct
        if not reg_run:
            hdulist = fits.open(pdf1d_file)
            delta1 = hdulist['M_ini'].data[-1,1] - hdulist['M_ini'].data[-1,0]
            if delta1 > 1.0:  # old linear spacing
                print('pdf1d lin mass spacing - full refitting needed')
                old_mass_spacing = True
            else:
                old_mass_spacing = False
                print('pdf1d log mass spacing - ok')

            if old_mass_spacing:
                run_done = False
                reg_run = True

        # now check if the number of results is the same as 
        #    the number of observations
        if not reg_run:
            # get the observed catalog
            obs = Table.read(cat_file)

            # get the fit results catalog
            t = Table.read(stats_file)
            # get the number of stars that have been fit
            indxs, = np.where(t['Pmax'] != 0.0)

            # get the number of entries in the lnp file
            f = tables.open_file(lnp_file, 'r')
            nlnp = f.root._v_nchildren - 2
            f.close()

            print('# obs, stats, lnp = ',len(obs), len(indxs), nlnp)
            if (len(indxs) == len(obs)) & (nlnp == len(obs)):

                # final check, is the pdf1d file correctly populated
                tot_prob = np.sum(hdulist['M_ini'].data, axis=1)
                tindxs, = np.where(tot_prob > 0.0)
                print('# good pdf1d = ', len(tindxs) - 1)
                if len(tindxs) == (len(obs) + 1):
                    run_done = True

        if run_done:
            print(stats_file + ' done')
        else:

            j += 1
            if j%n_pernode_files == 0:
                cur_f += 1

                # close previous files
                if j != 0:
                    pf.close()
                    print('total sed_trim size [Gb] = ', 
                          cur_total_size/(1024.*1024.*1024.))
                    cur_total_size = 0.0

                # open the slurm and param files
                pf_open = True
                joblist_file = job_path+'beast_batch_fit_'+str(cur_f) \
                               +'.joblist'
                pf = open(joblist_file,'w')
                

            ext_str = ''
            if reg_run:
                print(stats_file 
                      + ' does not exist ' + 
                      '- adding job as a regular fit job (not resume job)')
            else:
                print(stats_file 
                      + ' not done - adding to continue fitting list (' + \
                      str(len(indxs)) + '/' + str(len(t['Pmax'])) + ')')
                ext_str = '-r'

            nice_str = ''
            if nice is not None:
                nice_str = 'nice -n' + str(int(nice)) + ' '

            pipe_str = ' > '
            if not overwrite_logfile:
                pipe_str = ' >> '

            job_command = nice_str + 'python run_beast_production.py -f ' + ext_str + ' ' + \
                          sd_num + ' '+sub_num + pipe_str \
                          + log_path+'beast_fit' + \
                          '_sd'+sd_num+'_sub'+sub_num+'.log'

            pf.write(job_command+'\n')

    if pf_open:
        pf.close()




if __name__ == '__main__':

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("projectname",
                        help="project name to use (basename for files)")
    parser.add_argument("datafile",
                        help="file with the observed data (FITS file)")
    parser.add_argument("--num_percore", default=5, type=int,
                        help="number of fitting runs per core")
    args = parser.parse_args()

    project = args.projectname
    datafile = args.datafile
    n_pernode_files = args.num_percore

    setup_batch_beast_fit(project, datafile, n_pernode_files=5)
