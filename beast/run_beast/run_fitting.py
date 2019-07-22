# system imports
from __future__ import (absolute_import, division, print_function)
import os
import argparse
import time
import numpy as np
import glob
import pickle


# BEAST imports
from beast.fitting import fit
from beast.physicsmodel.grid import FileSEDGrid
from beast.tools import verify_params
from beast.run_beast.helper_functions import (subcatalog_fname,
                                                 parallel_wrapper,
                                                 get_modelsubgridfiles)


import datamodel
import importlib

import pdb



def run_fitting(use_sd=True, nsubs=1, nprocs=1,
                choose_sd_sub=None, choose_subgrid=None,
                resume=False):
    """
    Run the fitting.  If nsubs > 1, this will find existing subgrids.
    If use_sd is True, will also incorporate source density info.

    The additional choose_* options are to make queue scripts usable,
    by specifying a given SD+sub and/or subgrid for the fitting run.


    Parameters
    ----------
    use_sd : boolean (default=True)
        If True, create source density dependent noise models (determined by
        finding matches to datamodel.astfile with SD info)

    nsubs : int (default=1)
        number of subgrids used for the physics model

    nprocs : int (default=1)
        Number of parallel processes to use
        (currently only implemented for subgrids)

    choose_sd_sub : list of two strings (default=None)
        If this is set, the fitting will just be for this combo of SD+sub,
        rather than all of them.  Overrides use_sd.
        format of the list: ['#-#','#']

    choose_subgrid : int (default=None)
        If this is set, the fitting with just be for this subgrid index.
        If nsubs=1, this is ignored.

    resume : boolean (default=False)
        choose whether to resume existing run or start over

    """

    # before doing ANYTHING, force datamodel to re-import (otherwise, any
    # changes within this python session will not be loaded!)
    importlib.reload(datamodel)
    # check input parameters
    verify_params.verify_input_format(datamodel)


    # keep track of time
    start_time = time.clock()



    # --------------------
    # make lists of file names
    # --------------------

    file_dict = create_filenames(use_sd=use_sd, nsubs=nsubs,
                                     choose_sd_sub=choose_sd_sub,
                                     choose_subgrid=choose_subgrid)    
    # input files
    photometry_files = file_dict['photometry_files']
    modelsedgrid_files = file_dict['modelsedgrid_files']
    noise_files = file_dict['noise_files']

    # output files
    stats_files = file_dict['stats_files']
    pdf_files = file_dict['pdf_files']
    lnp_files = file_dict['lnp_files']

    # total number of files
    n_files = len(photometry_files)

    # other potentially useful info
    sd_sub_info = file_dict['sd_sub_info']
    gridsub_info = file_dict['gridsub_info']


    if nsubs > 1:

        gridpickle_files = file_dict['gridpickle_files']

        # make the grid dictionary file:
        # File where the ranges and number of unique values for the grid
        # will be stored (this can take a while to calculate)
        for i in range(len(gridpickle_files)):
            if not os.path.isfile(gridpickle_files[i]):
                # list of corresponding SED grids and noise models
                modelsedgrid_list = sorted(glob.glob(
                    os.path.join(datamodel.project, '**/*gridsub'+str(gridsub_info[i])+'_sed_trim*' ) ))
                noise_list = sorted(glob.glob(
                    os.path.join(datamodel.project, '**/*gridsub'+str(gridsub_info[i])+'_noisemodel_trim*' ) ))
                print('creating grid_info_dict for ' + gridpickle_files[i])
                grid_info_dict = subgridding_tools.reduce_grid_info(
                    modelsedgrid_list,
                    noise_list, nprocs=nprocs)

                with open(gridpickle_files[i], 'wb') as p:
                    pickle.dump(grid_info_dict, p)
                print('wrote grid_info_dict to ' + gridpickle_files[i])




    # --------------------
    # do the fitting!
    # --------------------

    # set up function inputs
    
    if nsubs == 1:

        input_list = [(photometry_files[i], modelsedgrid_files[i], noise_files[i],
                        stats_files[i], pdf_files[i], lnp_files[i], None, resume)
                            for i in range(n_files)]

    if nsubs > 1:

        input_list = [(photometry_files[i], modelsedgrid_files[i], noise_files[i],
                        stats_files[i], pdf_files[i], lnp_files[i], gridpickle_files[i], resume)
                            for i in range(n_files)]


    # run the fitting (via parallel wrapper)

    parallel_wrapper(fit_submodel, input_list, nprocs=nprocs)


    # see how long it took!
    new_time = time.clock()
    print('time to fit: ',(new_time - start_time)/60., ' min')







def create_filenames(use_sd=True, nsubs=1,
                     choose_sd_sub=None, choose_subgrid=None):
    """
    Helper function to make all of the filenames.  Arguments are same as run_fitting.

    Returns dictionary with the lists of filenames.
    """

    # input files
    photometry_files = []
    modelsedgrid_files = []
    noise_files = []

    # output files
    stats_files = []
    pdf_files = []
    lnp_files = []

    # other potentially useful things
    sd_sub_info = []
    gridsub_info = []
    

    

    # ** no subgrids **

    if nsubs == 1:

        # -- SD+sub specified
        if choose_sd_sub is not None:

            photometry_files.append( datamodel.obsfile.replace('.fits',
                '_SD{0}_sub{1}.fits'.format(choose_sd_sub[0], choose_sd_sub[1])) )
            modelsedgrid_files.append( '{0}/{0}_SD{1}_sub{2}_seds_trim.grid.hd5'.format(
                datamodel.project, choose_sd_sub[0], choose_sd_sub[1]) )
            noise_files.append( '{0}/{0}_SD{1}_sub{2}_noisemodel_trim.hd5'.format(
                datamodel.project, choose_sd_sub[0], choose_sd_sub[1]) )

            stats_files.append( '{0}/{0}_SD{1}_sub{2}_stats.fits'.format(
                datamodel.project, choose_sd_sub[0], choose_sd_sub[1]) )
            pdf_files.append( '{0}/{0}_SD{1}_sub{2}_pdf1d.fits'.format(
                datamodel.project, choose_sd_sub[0], choose_sd_sub[1]) )
            lnp_files.append( '{0}/{0}_SD{1}_sub{2}_lnp.hd5'.format(
                datamodel.project, choose_sd_sub[0], choose_sd_sub[1]) )

            sd_sub_info.append([choose_sd_sub[0], choose_sd_sub[1]])
              
        # -- using source density info
        elif use_sd == True:

            photometry_files = sorted( glob.glob(data.obsfile.replace('.fits','_SD*_sub*.fits')) )

            for phot_file in photometry_files:
                # get the sd/sub number
                dpos = phot_file.find('SD')
                spos = phot_file.find('sub')
                ppos = phot_file.rfind('.')
                curr_sd = phot_file[dpos+2:spos-1]
                curr_sub = phot_file[spos+3:ppos]
                # construct other file names
                modelsedgrid_files.append( '{0}/{0}_SD{1}_sub{2}_seds_trim.grid.hd5'.format(
                    datamodel.project, curr_sd, curr_sub) )
                noise_files.append( '{0}/{0}_SD{1}_sub{2}_noisemodel_trim.hd5'.format(
                    datamodel.project, curr_sd, curr_sub) )

                stats_files.append( '{0}/{0}_SD{1}_sub{2}_stats.fits'.format(
                    datamodel.project, curr_sd, curr_sub) )
                pdf_files.append( '{0}/{0}_SD{1}_sub{2}_pdf1d.fits'.format(
                    datamodel.project, curr_sd, curr_sub) )
                lnp_files.append( '{0}/{0}_SD{1}_sub{2}_lnp.hd5'.format(
                    datamodel.project, curr_sd, curr_sub) )

                sd_sub_info.append([curr_sd, curr_sub])
 

        # -- no source density splitting
        else:
    
            photometry_files.append( datamodel.obsfile )
            modelsedgrid_files.append( '{0}/{0}_seds_trim.grid.hd5'.format(datamodel.project) )
            noise_files.append( '{0}/{0}_noisemodel_trim.hd5'.format(datamodel.project) )

            stats_files.append( '{0}/{0}_stats.fits'.format(datamodel.project) )
            pdf_files.append( '{0}/{0}_pdf1d.fits'.format(datamodel.project) )
            lnp_files.append( '{0}/{0}_lnp.hd5'.format(datamodel.project) )


    # ** with subgrids **

    # subgrids require a pickle file with grid info
    gridpickle_files = []

    
    if nsubs > 1:

        # start with getting the model grid files (note these aren't trimmed ones)
        outdir = os.path.join('.', datamodel.project)
        subgrid_names_file = os.path.join(outdir, 'subgrid_fnames.txt')
        temp = get_modelsubgridfiles(subgrid_names_file)
        gridsub_list = np.arange(len(temp))
        # or a subset if set
        if choose_subgrid is not None:
            gridsub_list = [choose_subgrid]

        # -- SD+sub specified
        if choose_sd_sub is not None:

            photometry_files.append( datamodel.obsfile.replace('.fits',
                '_SD{0}_sub{1}.fits'.format(choose_sd_sub[0], choose_sd_sub[1])) )

            for gridsub in gridsub_list:
                modelsedgrid_files.append(
                    '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}_gridsub{3}_sed_trim.grid.hd5'.format(
                    datamodel.project, choose_sd_sub[0], choose_sd_sub[1], gridsub) )
                noise_files.append(
                    '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}_gridsub{3}_noisemodel_trim.hd5'.format(
                    datamodel.project, choose_sd_sub[0], choose_sd_sub[1], gridsub) )

                stats_files.append(
                    '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}_gridsub{3}_stats.fits'.format(
                    datamodel.project, choose_sd_sub[0], choose_sd_sub[1], gridsub) )
                pdf_files.append(
                    '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}_gridsub{3}_pdf1d.fits'.format(
                    datamodel.project, choose_sd_sub[0], choose_sd_sub[1], gridsub) )
                lnp_files.append(
                    '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}_gridsub{3}_lnp.hd5'.format(
                    datamodel.project, choose_sd_sub[0], choose_sd_sub[1], gridsub) )

                gridpickle_files.append('{0}/grid_info_dict_gridsub{1}.pkl'.format(
                    datamodel.project, gridsub) )

                sd_sub_info.append([choose_sd_sub[0], choose_sd_sub[1]])
                gridsub_info.append(gridsub)
                
        # -- using source density info
        elif use_sd == True:

            phot_file_list = sorted( glob.glob(datamodel.obsfile.replace('.fits','_SD*_sub*.fits')) )

            for phot_file in phot_file_list:
                # get the sd/sub number
                dpos = phot_file.find('SD')
                spos = phot_file.find('sub')
                ppos = phot_file.rfind('.')
                curr_sd = phot_file[dpos+2:spos-1]
                curr_sub = phot_file[spos+3:ppos]
                
                # construct other file names
                for gridsub in gridsub_list:
                    photometry_files.append(phot_file)
                    modelsedgrid_files.append(
                        '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}_gridsub{3}_sed_trim.grid.hd5'.format(
                        datamodel.project, curr_sd, curr_sub, gridsub) )
                    noise_files.append(
                        '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}_gridsub{3}_noisemodel_trim.hd5'.format(
                        datamodel.project, curr_sd, curr_sub, gridsub) )

                    stats_files.append(
                        '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}_gridsub{3}_stats.fits'.format(
                        datamodel.project, curr_sd, curr_sub, gridsub) )
                    pdf_files.append(
                        '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}_gridsub{3}_pdf1d.fits'.format(
                        datamodel.project, curr_sd, curr_sub, gridsub) )
                    lnp_files.append(
                        '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}_gridsub{3}_lnp.hd5'.format(
                        datamodel.project, curr_sd, curr_sub, gridsub) )
                    
                    gridpickle_files.append('{0}/grid_info_dict_gridsub{1}.pkl'.format(
                        datamodel.project, gridsub) )

                    sd_sub_info.append([curr_sd, curr_sub])
                    gridsub_info.append(gridsub)

        # -- no source density splitting
        else:

            for gridsub in gridsub_list:
                photometry_files.append( datamodel.obsfile )
                modelsedgrid_files.append( '{0}/{0}_gridsub{1}_seds_trim.grid.hd5'.format(
                    datamodel.project, gridsub) )
                noise_files.append( '{0}/{0}_gridsub{1}_noisemodel_trim.hd5'.format(
                    datamodel.project, gridsub) )

                stats_files.append( '{0}/{0}_gridsub{1}_stats.fits'.format(
                    datamodel.project, gridsub) )
                pdf_files.append( '{0}/{0}_gridsub{1}_pdf1d.fits'.format(
                    datamodel.project, gridsub) )
                lnp_files.append( '{0}/{0}_gridsub{1}_lnp.hd5'.format(
                    datamodel.project, gridsub) )
    
                gridpickle_files.append('{0}/grid_info_dict_gridsub{1}.pkl'.format(
                    datamodel.project, gridsub) )

                gridsub_info.append(gridsub)
            
        
    # double check that all file lists are the same length
    n_file_list = [len(x) for x in
                    [photometry_files, modelsedgrid_files, noise_files,
                     stats_files, pdf_files, lnp_files]]
    if len(np.unique(n_file_list)) > 1:
        print("file list lengths don't match!")
        return None




    return {'photometry_files':photometry_files,
                'modelsedgrid_files':modelsedgrid_files,
                'noise_files':noise_files,
                'stats_files':stats_files,
                'pdf_files':pdf_files,
                'lnp_files':lnp_files,
                'gridpickle_files':gridpickle_files,
                'sd_sub_info':sd_sub_info,
                'gridsub_info':gridsub_info}


def fit_submodel(photometry_file, modelsedgrid_file, noise_file,
                     stats_file, pdf_file, lnp_file,
                     grid_info_dict=None, resume=False):
    """
    Code to run the SED fitting

    Parameters
    ----------
    photometry_file : string
        path+name of the photometry file

    modelsedgrid_file : string
        path+name of the physics model grid file

    noise_file : string
        path+name of the noise model file

    stats_file : string
        path+name of the file to contain stats output

    pdf_file : string
        path+name of the file to contain 1D PDF output

    lnp_file : string
        path+name of the file to contain log likelihood output

    grid_info_file : string (default=None)
        path+name for pickle file that contains dictionary with subgrid
        min/max/n_unique (required for a run with subgrids)

    resume : boolean (default=False)
        choose whether to resume existing run or start over


    Returns
    -------
    noisefile : string
        name of the created noise file
        
    """

    # read in the photometry catalog
    obsdata = datamodel.get_obscat(photometry_file, datamodel.filters)

    # check if it's a subgrid run by looking in the file name
    if 'gridsub' in modelsedgrid_file:
        subgrid_run = True
        print('loading grid_info_dict from ' + gridpickle_files[i])
        with open(grid_info_file, 'rb') as p:
            grid_info_dict = pickle.loads(p.read())
    else:
        subgrid_run = False
        

    # load the SED grid and noise model
    modelsedgrid = FileSEDGrid(modelsedgrid_file)
    noisemodel_vals = noisemodel.get_noisemodelcat(noise_file)


    if subgrid_run == True:
        try:
            fit.summary_table_memory(obsdata, noisemodel_vals,
                                         modelsedgrid,
                                         resume=resume,
                                         threshold=-10.,
                                         save_every_npts=100,
                                         lnp_npts=500,
                                         stats_outname=stats_file,
                                         pdf1d_outname=pdf_file,
                                         grid_info_dict=grid_info_dict,
                                         lnp_outname=lnp_file,
                                         do_not_normalize=True,
                                         surveyname=datamodel.surveyname)
            print("Done fitting on grid " + modelsedgrid_file)
        except Exception as e:
            if not args.ignore_missing_subresults:
                raise e

    if subgrid_run == False:
        
        fit.summary_table_memory(obsdata, noisemodel_vals,
                                         modelsedgrid,
                                         resume=resume,
                                         threshold=-10.,
                                         save_every_npts=100,
                                         lnp_npts=500,
                                         stats_outname=stats_file,
                                         pdf1d_outname=pdf_file,
                                         lnp_outname=lnp_file,
                                         surveyname=datamodel.surveyname)
        print("Done fitting on grid " + modelsedgrid_file)



if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--use_sd",
                        help="create source density dependent noise models",
                        action="store_true")
    parser.add_argument("--nsubs", type=int, default=1,
                        help='''Number of subgrids that the physics model
                        was split into''')
    parser.add_argument("--nprocs", type=int, default=1,
                        help='''Number of processes to use to process
                        the subgrids''')
    parser.add_argument("--choose_sd_sub", nargs='+', default=None,
                        help="Fit just this combo of SD+sub. Format: ['#-#','#']")
    parser.add_argument("--choose_subgrid", type=int, default=None,
                        help="Fit just this subgrid number")
    parser.add_argument("-r","--resume",
                        help='resume a fitting run',
                        action="store_true")

    args = parser.parse_args()


    run_fitting(use_sd=args.use_sd,
                nsubs=args.nsubs,
                nprocs=args.nprocs,
                choose_sd_sub=args.choose_sd_sub,
                choose_subgrid=args.choose_subgrid,
                resume=args.resume)
    
    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()

