# system imports
from __future__ import (absolute_import, division, print_function)
import os
import argparse
import glob


# BEAST imports
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel
from beast.physicsmodel.grid import FileSEDGrid
from beast.tools import verify_params
from beast.tools.run.helper_functions import (parallel_wrapper,
                                              get_modelsubgridfiles)


import datamodel
import importlib

#import pdb



def create_obsmodel(use_sd=True, nsubs=1, nprocs=1, subset=[None,None],
                        use_rate=True):
    """
    Create the observation models.  If nsubs > 1, this will find existing
    subgrids.  If use_sd is True, will also incorporate source density
    info.


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

    subset : list of two ints (default=[None,None])
        Only process subgrids in the range [start,stop].
        (only relevant if nsubs > 1)

    use_rate : boolean (default=True)
        Choose whether to use the rate or magnitude when creating the noise
        model.  This should always be True, but is currently an option to be
        compatible with the phat_small example (which has no rate info).
        When that gets fixed, please remove this option!

    """
    
    # before doing ANYTHING, force datamodel to re-import (otherwise, any
    # changes within this python session will not be loaded!)
    importlib.reload(datamodel)
    # check input parameters
    verify_params.verify_input_format(datamodel)


    # --------------------
    # figure out if there are source density bins
    # --------------------

    ast_file_list = sorted(glob.glob(datamodel.astfile.replace('.fits','*_SD*')))
    
    if (use_sd == True) and (len(ast_file_list) > 0):
        
        sd_list = []
        for ast_file in ast_file_list:
            dpos = ast_file.find('SD')
            ppos = ast_file.rfind('.')
            sd_list.append( ast_file[dpos+2:ppos] )
        print('sd list: ', sd_list)

    else:
        # if there are no ASTs with source densities, the flag should be "false"
        use_sd = False


    # --------------------
    # no subgrids
    # --------------------
    
    if nsubs == 1:

        modelsedgridfile = "{0}/{0}_seds.grid.hd5".format(datamodel.project)
        
        # if we're splitting by source density
        if use_sd == True:

            input_list = [(modelsedgridfile,curr_sd) for curr_sd in sd_list]
                
            parallel_wrapper(gen_obsmodel, input_list, nprocs=nprocs)


        # if we're not splitting by source density
        if use_sd == False:

            input_list = [(modelsedgridfile, None, use_rate)]
            
            parallel_wrapper(gen_obsmodel, input_list, nprocs=nprocs)




    # --------------------
    # use subgrids
    # --------------------

    if nsubs > 1:
    
        # get the list of physics model files
        outdir = os.path.join('.', datamodel.project)
        subgrid_names_file = os.path.join(outdir, 'subgrid_fnames.txt')
        modelsedgridfiles = get_modelsubgridfiles(subgrid_names_file)[ slice(subset[0], subset[1]) ]

        
        # if we're splitting by source density
        if use_sd == True:

            input_list = [(sedfile,curr_sd) for sedfile in modelsedgridfiles for curr_sd in sd_list]
                
            parallel_wrapper(gen_obsmodel, input_list, nprocs=nprocs)


        # if we're not splitting by source density
        if use_sd == False:

            input_list = [(sedfile, None) for sedfile in modelsedgridfiles]
            
            parallel_wrapper(gen_obsmodel, input_list, nprocs=nprocs)

     




def gen_obsmodel(modelsedgridfile, source_density=None, use_rate=True):
    """
    Code to create filenames and run the toothpick noise model

    Parameters
    ----------
    modelsedgridfile : string
        path+name of the physics model grid file

    source_density : string (default=None)
        set to None if there's no source density info, otherwise set to
        a string of the form "#-#"

    use_rate : boolean (default=True)
        Choose whether to use the rate or magnitude when creating the noise
        model.  This should always be True, but is currently an option to be
        compatible with the phat_small example (which has no rate info).
        When that gets fixed, please remove this option!

    Returns
    -------
    noisefile : string
        name of the created noise file
        
    """

    print('')
            
    # noise and AST file names
    noisefile = modelsedgridfile.replace('seds', 'noisemodel')
    astfile = datamodel.astfile
    
    # If we are treating regions with different
    # backgrounds/source densities separately, pick one of the
    # split ast files, and name noise file accordingly
    if source_density is not None:
        noisefile = noisefile.replace('noisemodel','noisemodel_SD'+source_density)
        astfile = datamodel.astfile.replace('.fits','_SD'+source_density.replace('_','-')+'.fits')


    # only create noise file if it doesn't already exist
    if not os.path.isfile(noisefile):

        print('creating '+noisefile)
            
        modelsedgrid = FileSEDGrid(modelsedgridfile)

        noisemodel.make_toothpick_noise_model(
                    noisefile,
                    astfile,
                    modelsedgrid,
                    absflux_a_matrix=datamodel.absflux_a_matrix,
                    use_rate=use_rate)

    else:
        print(noisefile+' already exists')
        
    return noisefile # (same as noisefile)





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
    parser.add_argument("--subset", type=int, nargs=2, default=[None, None],
                        help='''Only process subgrids in the range
                        [start, stop].''')

    args = parser.parse_args()


    create_obsmodel(use_sd=args.use_sd,
                        nsubs=args.nsubs,
                        nprocs=args.nprocs,
                        subset=args.subset)

    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()
