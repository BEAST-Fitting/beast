#!/usr/bin/env python

"""
Code to create the batch files for fitting
"""

from __future__ import print_function
import os
import argparse
import tables
#####
import numpy as np
from astropy.table import Table
from astropy.io import fits
#####

from beast.tools import verify_params
from beast.run_beast import create_filenames
import datamodel
import importlib


def setup_batch_beast_fit(num_percore=5,
                          nice=None,
                          overwrite_logfile=True,
                          prefix=None,
                          use_sd=True,
                          nsubs=1, nprocs=1):
    """
    Sets up batch files for submission to the 'at' queue on
    linux (or similar) systems

    Parameters
    ----------
    num_percore : int (default = 5)
        number of fitting runs per core

    nice : int (default = None)
        set this to an integer (-20 to 20) to prepend a "nice" level
        to the fitting command

    overwrite_logfile : boolean (default = True)
        if True, will overwrite the log file; if False, will append to
        existing log file

    prefix : string (default=None)
        Set this to a string (such as 'source activate astroconda') to prepend
        to each batch file (use '\n's to make multiple lines)

    use_sd : boolean (default=True)
        If True, split runs based on source density (determined by finding
        matches to datamodel.astfile with SD info)

    nsubs : int (default=1)
        number of subgrids used for the physics model

    nprocs : int (default=1)
        Number of parallel processes to use when doing the fitting
        (currently only implemented for subgrids)


    Returns
    -------
    run_info_dict : dict
        Dictionary indicating which catalog files have complete modeling, and
        which job files need to be run

    """

    # before doing ANYTHING, force datamodel to re-import (otherwise, any
    # changes within this python session will not be loaded!)
    importlib.reload(datamodel)
    # check input parameters
    verify_params.verify_input_format(datamodel)


    # setup the subdirectory for the batch and log files
    job_path = datamodel.project+'/fit_batch_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path+'logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    
    # get file name lists (to check if they exist and/or need to be resumed)
    file_dict = create_filenames.create_filenames(use_sd=use_sd, nsubs=nsubs)    

    # - input files
    photometry_files = file_dict['photometry_files']
    #modelsedgrid_files = file_dict['modelsedgrid_files']
    #noise_files = file_dict['noise_files']

    # - output files
    stats_files = file_dict['stats_files']
    pdf_files = file_dict['pdf_files']
    lnp_files = file_dict['lnp_files']

    # - total number of files
    n_files = len(photometry_files)

    # - other useful info
    sd_sub_info = file_dict['sd_sub_info']
    gridsub_info = file_dict['gridsub_info']

    
    # names of output log files
    log_files = []
    
    for i in range(n_files):
        
        sd_piece = ''
        if use_sd == True:
            sd_piece = '_SD'+sd_sub_info[i][0]+'_sub'+sd_sub_info[i][1]

        gridsub_piece = ''
        if nsubs > 1:
           gridsub_piece = '_gridsub'+str(gridsub_info[i])

        log_files.append('beast_fit'+sd_piece+gridsub_piece+'.log')



    # start making the job files!

    pf_open = False
    cur_f = 0
    cur_total_size = 0.0
    j = -1

    # keep track of which files are done running
    run_info_dict = {'phot_file': photometry_files,
                     'done': np.full(n_files, False),
                     'files_to_run': []}


    for i, phot_file in enumerate(photometry_files):

        print('')

        # check if this is a full run
        reg_run = False
        run_done = False
        if not os.path.isfile(stats_files[i]):
            reg_run = True
            print('no stats file')
        if not os.path.isfile(pdf_files[i]):
            reg_run = True
            print('no pdf1d file')
        if not os.path.isfile(lnp_files[i]):
            reg_run = True
            print('no lnp file')

        # first check if the pdf1d mass spacing is correct
        if not reg_run:
            hdulist = fits.open(pdf_files[i])
            delta1 = (hdulist['M_ini'].data[-1, 1]
                      - hdulist['M_ini'].data[-1, 0])
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
            obs = Table.read(photometry_files[i])

            # get the fit results catalog
            t = Table.read(stats_files[i])
            # get the number of stars that have been fit
            indxs, = np.where(t['Pmax'] != 0.0)

            # get the number of entries in the lnp file
            f = tables.open_file(lnp_files[i], 'r')
            nlnp = f.root._v_nchildren - 2
            f.close()

            print('# obs, stats, lnp = ', len(obs), len(indxs), nlnp)
            if (len(indxs) == len(obs)) & (nlnp == len(obs)):

                # final check, is the pdf1d file correctly populated
                tot_prob = np.sum(hdulist['M_ini'].data, axis=1)
                tindxs, = np.where(tot_prob > 0.0)
                print('# good pdf1d = ', len(tindxs) - 1)
                if len(tindxs) == (len(obs) + 1):
                    run_done = True

        if run_done:
            print(stats_files[i] + ' done')
            run_info_dict['done'][i] = True
        else:

            j += 1
            if j % num_percore == 0:
                cur_f += 1

                # close previous files
                if j != 0:
                    pf.close()
                    print('total sed_trim size [Gb] = ',
                          cur_total_size/(1024.*1024.*1024.))
                    cur_total_size = 0.0

                # open the slurm and param files
                pf_open = True
                joblist_file = (job_path + 'beast_batch_fit_' + str(cur_f)
                                + '.joblist')
                pf = open(joblist_file, 'w')
                run_info_dict['files_to_run'].append(joblist_file)

                # write out anything at the beginning of the file
                if prefix is not None:
                    pf.write(prefix+'\n')


            # flag for resuming
            resume_str = ''
            if reg_run:
                print(stats_files[i]
                      + ' does not exist ' +
                      '- adding job as a regular fit job (not resume job)')
            else:
                print(stats_files[i]
                      + ' not done - adding to continue fitting list ('
                      + str(len(indxs)) + '/' + str(len(t['Pmax'])) + ')')
                resume_str = '-r'

            # prepend a `nice` value
            nice_str = ''
            if nice is not None:
                nice_str = 'nice -n' + str(int(nice)) + ' '

            # choose whether to append or overwrite log file
            pipe_str = ' > '
            if not overwrite_logfile:
                pipe_str = ' >> '

            # set SD+sub option
            sd_str = ''
            if use_sd == True:
                sd_str = ' --choose_sd_sub "{0}" "{1}" '.format(sd_sub_info[i][0], sd_sub_info[i][1])

            # set gridsub option
            gs_str = ''
            if nsubs > 1:
                gs_str = ' --choose_subgrid {0} '.format(gridsub_info[i])

            job_command = (nice_str + 'python -m beast.run_beast.run_fitting '
                           + resume_str + sd_str + gs_str
                           + ' --nsubs '+str(nsubs)
                           + ' --nprocs '+str(nprocs)
                           + pipe_str + log_path+log_files[i])

            pf.write(job_command+'\n')

    if pf_open:
        pf.close()

    # return the info about completed modeling
    return run_info_dict


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
    num_percore = args.num_percore

    setup_batch_beast_fit(project, datafile, num_percore=num_percore)
