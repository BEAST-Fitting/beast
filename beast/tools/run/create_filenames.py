# system imports
from __future__ import (absolute_import, division, print_function)
import os
import numpy as np
import glob

# BEAST imports
from beast.tools import verify_params
from beast.tools.run.helper_functions import get_modelsubgridfiles

import datamodel
import importlib


def create_filenames(use_sd=True, nsubs=1,
                     choose_sd_sub=None, choose_subgrid=None):
    """
    Helper function to make all of the filenames.  SED grid and noise model
    are trimmed versions.

    Parameters
    ----------
    use_sd : boolean (default=True)
        If True, create source density dependent noise models (determined by
        finding matches to datamodel.astfile with SD info)

    nsubs : int (default=1)
        number of subgrids used for the physics model

    choose_sd_sub : list of two strings (default=None)
        If this is set, the fitting will just be for this combo of SD+sub,
        rather than all of them.  Overrides use_sd.
        format of the list: ['#-#','#']

    choose_subgrid : int (default=None)
        If this is set, the fitting with just be for this subgrid index.
        If nsubs=1, this is ignored.

    Returns
    -------
    dictionary with the lists of filenames, plus the corresponding SD+sub and
    gridsub values for easy referencing

    """

    # before doing ANYTHING, force datamodel to re-import (otherwise, any
    # changes within this python session will not be loaded!)
    importlib.reload(datamodel)
    # check input parameters
    verify_params.verify_input_format(datamodel)


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
            modelsedgrid_files.append( '{0}/{0}_SD{1}_sub{2}_sed_trim.grid.hd5'.format(
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

            photometry_files = sorted( glob.glob(datamodel.obsfile.replace('.fits','_SD*_sub*.fits')) )

            for phot_file in photometry_files:
                # get the sd/sub number
                dpos = phot_file.find('SD')
                spos = phot_file.find('sub')
                ppos = phot_file.rfind('.')
                curr_sd = phot_file[dpos+2:spos-1]
                curr_sub = phot_file[spos+3:ppos]
                # construct other file names
                modelsedgrid_files.append( '{0}/{0}_SD{1}_sub{2}_sed_trim.grid.hd5'.format(
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
            modelsedgrid_files.append( '{0}/{0}_sed_trim.grid.hd5'.format(datamodel.project) )
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
        # use that to get the number of subgrids and make a list of them
        gridsub_list = np.arange(len(temp))
        # or a subset if set
        if choose_subgrid is not None:
            gridsub_list = [choose_subgrid]

        # -- SD+sub specified
        if choose_sd_sub is not None:

            for gridsub in gridsub_list:
                
                photometry_files.append( datamodel.obsfile.replace('.fits',
                    '_SD{0}_sub{1}.fits'.format(choose_sd_sub[0], choose_sd_sub[1])) )

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

                gridpickle_files.append('{0}/SD{1}_sub{2}/grid_info_dict.pkl'.format(
                    datamodel.project, choose_sd_sub[0], choose_sd_sub[1]) )

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
                    
                    gridpickle_files.append('{0}/SD{1}_sub{2}/grid_info_dict.pkl'.format(
                        datamodel.project, curr_sd, curr_sub) )

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
    
                gridpickle_files.append('{0}/grid_info_dict.pkl'.format(
                    datamodel.project) )

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
