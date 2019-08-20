#!/usr/bin/env python
"""
Script to run the BEAST on the PHAT-like data.
Assumes that the datamodel.py file exists in the same directory as this script.
  And it must be called datamodel.py 
"""

# system imports
from __future__ import (absolute_import, division, print_function)
import os
import argparse

from astropy import constants as const

# BEAST imports
from beast.physicsmodel.create_project_dir import create_project_dir
from beast.physicsmodel.model_grid import (make_iso_table,
                                           make_spectral_grid,
                                           add_stellar_priors,
                                           make_extinguished_sed_grid)
from beast.physicsmodel.grid import FileSEDGrid
from beast.tools.run.helper_functions import parallel_wrapper

from beast.tools import verify_params
from beast.tools import subgridding_tools

import datamodel
import importlib

#import pdb


def create_physicsmodel(nsubs=1, nprocs=1, subset=[None,None]):
    """
    Create the physics model grid.  If nsubs > 1, this will make sub-grids.


    Parameters
    ----------
    nsubs : int (default=1)
        number of subgrids to split the physics model into

    nprocs : int (default=1)
        Number of parallel processes to use
        (currently only implemented for subgrids)

    subset : list of two ints (default=[None,None])
        Only process subgrids in the range [start,stop].
        (only relevant if nsubs > 1)

    """

    # before doing ANYTHING, force datamodel to re-import (otherwise, any
    # changes within this python session will not be loaded!)
    importlib.reload(datamodel)
    # check input parameters
    verify_params.verify_input_format(datamodel)


    # filename for the SED grid
    modelsedgrid_filename = "%s/%s_seds.grid.hd5"%(datamodel.project,
                                                   datamodel.project)

    # grab the current subgrid slice
    subset_slice = slice(subset[0], subset[1])


    # make sure the project directory exists
    create_project_dir(datamodel.project)


    
    # download and load the isochrones
    (iso_fname, oiso) = make_iso_table(datamodel.project,
                                        oiso=datamodel.oiso,
                                        logtmin=datamodel.logt[0],
                                        logtmax=datamodel.logt[1],
                                        dlogt=datamodel.logt[2],
                                        z=datamodel.z)

    if hasattr(datamodel, 'add_spectral_properties_kwargs'):
        extra_kwargs = datamodel.add_spectral_properties_kwargs
    else:
        extra_kwargs = None

    if hasattr(datamodel, 'velocity'):
        redshift = (datamodel.velocity / const.c).decompose().value
    else:
        redshift = 0

    # generate the spectral library (no dust extinction)
    (spec_fname, g_spec) = make_spectral_grid(
        datamodel.project,
        oiso,
        osl=datamodel.osl,
        redshift=redshift,
        distance=datamodel.distances,
        distance_unit=datamodel.distance_unit,
        extLaw=datamodel.extLaw,
        add_spectral_properties_kwargs=extra_kwargs)

    # add the stellar priors as weights
    #   also computes the grid weights for the stellar part
    (pspec_fname, g_pspec) = add_stellar_priors(datamodel.project, g_spec)


    # --------------------
    # no subgrids
    # --------------------
    
    if nsubs == 1:
        # generate the SED grid by integrating the filter response functions
        #   effect of dust extinction applied before filter integration
        #   also computes the dust priors as weights
        make_extinguished_sed_grid(
            datamodel.project,
            g_pspec,
            datamodel.filters,
            extLaw=datamodel.extLaw,
            av=datamodel.avs,
            rv=datamodel.rvs,
            fA=datamodel.fAs,
            rv_prior_model=datamodel.rv_prior_model,
            av_prior_model=datamodel.av_prior_model,
            fA_prior_model=datamodel.fA_prior_model,
            spec_fname=modelsedgrid_filename,
            add_spectral_properties_kwargs=extra_kwargs)


    # --------------------
    # use subgrids
    # --------------------

    if nsubs > 1:
        # Work with the whole grid up to there (otherwise, priors need a
        # rework - they don't like having only a subset of the parameter
        # space, especially when there's only one age for example)

        # Make subgrids, by splitting the spectral grid into equal sized pieces
        custom_sub_pspec = subgridding_tools.split_grid(pspec_fname, nsubs)

        file_prefix = '{0}/{0}_'.format(datamodel.project)


        # function to process the subgrids individually
        def gen_subgrid(i, sub_name):
            sub_g_pspec = FileSEDGrid(sub_name)
            sub_seds_fname = '{}seds.gridsub{}.hd5'.format(file_prefix, i)

            
            # generate the SED grid by integrating the filter response functions
            #   effect of dust extinction applied before filter integration
            #   also computes the dust priors as weights
            (sub_seds_fname, sub_g_seds) = make_extinguished_sed_grid(
                datamodel.project,
                sub_g_pspec,
                datamodel.filters,
                extLaw=datamodel.extLaw,
                av=datamodel.avs,
                rv=datamodel.rvs,
                fA=datamodel.fAs,
                rv_prior_model=datamodel.rv_prior_model,
                av_prior_model=datamodel.av_prior_model,
                fA_prior_model=datamodel.fA_prior_model,
                add_spectral_properties_kwargs=extra_kwargs,
                seds_fname=sub_seds_fname)

            return sub_seds_fname


        # run the above function
        par_tuples = [(i, sub_name)
                      for i, sub_name in enumerate(custom_sub_pspec)][subset_slice]

        parallel_wrapper(gen_subgrid, par_tuples, nprocs=nprocs)
        

        # Save a list of subgrid names that we expect to see
        required_names = ['{}seds.gridsub{}.hd5'.format(file_prefix, i)
                          for i in range(nsubs)]
        
        outdir = os.path.join('.', datamodel.project)
        subgrid_names_file = os.path.join(outdir, 'subgrid_fnames.txt')
        
        with open(subgrid_names_file, 'w') as fname_file:
            for fname in required_names:
                fname_file.write(fname + '\n')


if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--nsubs", type=int, default=1,
                        help='''Number of subgrids to split the physics
                        model into''')
    parser.add_argument("--nprocs", type=int, default=1,
                        help='''Number of processes to use to process
                        the subgrids''')
    parser.add_argument("--subset", type=int, nargs=2, default=[None, None],
                        help='''Only process subgrids in the range
                        [start, stop].''')

    args = parser.parse_args()


    create_physicsmodel(nsubs=args.nsubs,
                            nprocs=args.nprocs,
                            subset=args.subset)

    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()

