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
                           nice=None,
                           seds_fname=None,
                           prefix=None):
    """
    Sets up batch files for submission to the 'at' queue on
    linux (or similar) systems


    Parameters
    ----------
    project : string
        project name to use (basename for files)

    datafile : string
        file with the observed data (FITS file) - the observed photometry,
        not the sub-files

    astfile : string
        file with ASTs

    num_subtrim : int (default = 5)
        number of trim batch jobs

    nice : int (default = None)
        set this to an integer (-20 to 20) to prepend a "nice" level
        to the trimming command

    seds_fname : string (default = None)
        full filename to the SED grid

    prefix : string (default=None)
        Set this to a string (such as 'source activate astroconda') to prepend
        to each batch file (use '\n's to make multiple lines)

    """
    ast_file = astfile

    if seds_fname is None:
        full_model_filename = "%s/%s_seds.grid.hd5" % (project, project)
    else:
        full_model_filename = seds_fname

    # photometry files
    cat_files = sorted(glob.glob(datafile.replace('.fits', '*_sub*.fits')))
    n_cat_files = len(cat_files)


    # setup the subdirectory for the batch and log files
    job_path = project+'/trim_batch_jobs/'

    # construct file bases
    filebase_list = []
    for cat_file in cat_files:
        # get the sd/sub number
        dpos = cat_file.find('SD_')
        spos = cat_file.find('sub')
        ppos = cat_file.rfind('.')
        curr_sd = cat_file[dpos+3:spos-1]
        curr_sub = cat_file[spos+3:ppos]

        filebase_list.append("%s/%s_sd%s_sub%s"%(project,
                                                 project,
                                                 curr_sd,
                                                 curr_sub))

        
    # call the generic batch trim code
    generic_batch_trim(full_model_filename,
                       ["%s/%s_noisemodel.hd5"%(project,project)] * n_cat_files,
                       cat_files,
                       [ast_file] * n_cat_files,
                       filebase_list,
                       job_path=job_path,
                       file_prefix='BEAST',
                       num_subtrim=num_subtrim,
                       nice=nice,
                       prefix=prefix)
                       
                           
 

def generic_batch_trim(model_grid_file,
                           noise_model_files,
                           data_files,
                           ast_files,
                           trimmed_file_bases,
                           job_path=None,
                           file_prefix='BEAST',
                           num_subtrim=5,
                           nice=None,
                           prefix=None):
    """
    Sets up batch files for submission to the 'at' queue on
    linux (or similar) systems

    This is a more generic version of the above code, in which all
    files are specified at the input, rather than constructed.

    The parameters that are lists of strings need to be the same length.

    Parameters
    ----------
    model_grid_file : string
        full filename to the SED grid (note that having the same model grid
        file is what makes the batch trimming efficient)

    noise_model_files : list of strings
        full filename to the noise models

    data_files : list of strings
        file with the observed data (FITS file)

    ast_files : list of strings
        file with ASTs

    trimmed_file_bases : string
        prefix (path+name) for the eventual trimmed model grid and noise model

    job_path : string (default=None)
        path in which to store info about the jobs
    
    file_prefix : string (default='BEAST')
        prefix for the name of the file in which to store the trimming commands

    num_subtrim : int (default = 5)
        number of trim batch jobs

    nice : int (default = None)
        set this to an integer (-20 to 20) to prepend a "nice" level
        to the trimming command

    prefix : string (default=None)
        Set this to a string (such as 'source activate astroconda') to prepend
        to each batch file (use '\n's to make multiple lines)

    """


    # check that everything is the same length that needs to be
    keyword_list = [noise_model_files, data_files, ast_files, trimmed_file_bases]
    length_array = np.array( [len(key) for key in keyword_list ] )
    if len(np.unique(length_array)) != 1:
        print('input lists are not the same length!')
        return


    # number of photometry catalog files
    n_cat_files = len(data_files)

    # make sure n_subtrim_files isn't larger than the number of
    # catalog sub-files
    if num_subtrim > n_cat_files:
        num_subtrim = n_cat_files

    # max number of files per process
    n_per_subtrim = int(n_cat_files/num_subtrim)
    if n_cat_files % num_subtrim != 0:
        n_per_subtrim += 1

    print('# trim files per process = ', n_per_subtrim)

    # setup the subdirectory for the batch and log files
    if job_path is None:
        job_path = os.path.dirname(trimmed_file_bases[0])+'/trim_batch_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path+'logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)


    # prepend nice
    nice_str = ''
    if nice is not None:
        nice_str = 'nice -n' + str(int(nice)) + ' '

    joblist_file = job_path+file_prefix+'_batch_trim.joblist'
    pf = open(joblist_file, 'w')

    # write out anything at the beginning of the file
    if prefix is not None:
        pf.write(prefix+'\n')

    
    bt_f = []
    for i in range(num_subtrim):
        trimfile = job_path+file_prefix+'_' + str(i+1)
        bt_f.append(open(trimfile, 'w'))
        bt_f[-1].write(model_grid_file + '\n')
        pf.write(nice_str + 'python -m beast.tools.trim_many_via_obsdata '
                 + trimfile + ' > '
                 + log_path + file_prefix+'_trim_tr'+str(i+1)+'.log\n')
    pf.close()


    k = 0
    n_cur = 0
    for i in range(n_cat_files):
        bt_f[k].write(noise_model_files[i])
        bt_f[k].write(' ' + data_files[i])
        bt_f[k].write(' ' + ast_files[i])
        bt_f[k].write(' ' + trimmed_file_bases[i])
        bt_f[k].write('\n')
        n_cur += 1
        if n_cur >= n_per_subtrim:
            n_cur = 0
            k += 1

    for i in range(num_subtrim):
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

    setup_batch_beast_trim(project, datafile, ast_file,
                           num_subtrim=n_subtrim_files)
