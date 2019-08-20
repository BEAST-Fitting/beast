#!/usr/bin/env python
"""
Script to run the BEAST on the PHAT-like data.
Assumes that the datamodel.py file exists in the same directory as this script.
  And it must be called datamodel.py
     used in make_model.py code in physicsmodel, more recoding needed to remove
     this dependency
"""

# system imports
from __future__ import (absolute_import, division, print_function)
import sys
import argparse
import string

from astropy import units

# BEAST imports
from beast.tools.run import (create_physicsmodel,
                             make_ast_inputs,
                             create_obsmodel,
                             run_fitting)
from beast.physicsmodel.grid import FileSEDGrid
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel
from beast.fitting import trim_grid
from beast.tools import verify_params


import datamodel


if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--physicsmodel",
                        help="Generate the physics model grid",
                        action="store_true")
    parser.add_argument("-a", "--ast", help="Generate an input AST file",
                        action="store_true")
    parser.add_argument("-o", "--observationmodel",
                        help="Calculate the observation model (bias and noise)",
                        action="store_true")
    parser.add_argument("-t", "--trim",
                        help="Trim the physics and observation model grids",
                        action="store_true")
    parser.add_argument("-f", "--fit", help="Fit the observed data",
                        action="store_true")
    parser.add_argument("-r", "--resume", help="Resume a fitting run",
                        action="store_true")
    args = parser.parse_args()

    # check input parameters, print what is the problem, stop run_beast
    verify_params.verify_input_format(datamodel)

    if args.physicsmodel:

        create_physicsmodel.create_physicsmodel(nsubs=1, nprocs=1)

    if args.ast:

        make_ast_inputs.make_ast_inputs(flux_bin_method=False)

    if args.observationmodel:
        print('Generating noise model from ASTs and absflux A matrix')

        create_obsmodel.create_obsmodel(use_sd=False,
                                            nsubs=1,
                                            nprocs=1,
                                            use_rate=False)

        # in the absence of ASTs, the splinter noise model can be used
        # instead of the toothpick model above
        #  **warning** not very realistic
        # import beast.observationmodel.noisemodel.splinter as noisemodel
        #
        # modelsedgridfile = datamodel.project + '/' + datamodel.project + \
        #    '_seds.grid.hd5'
        # modelsedgrid = FileSEDGrid(modelsedgridfile)
        #
        # noisemodel.make_splinter_noise_model(
        #    datamodel.noisefile,
        #    modelsedgrid,
        #    absflux_a_matrix=datamodel.absflux_a_matrix)

    if args.trim:
        print('Trimming the model and noise grids')

        # read in the observed data
        obsdata = datamodel.get_obscat(datamodel.obsfile,
                                       datamodel.filters)

        # get the modesedgrid on which to generate the noisemodel
        modelsedgridfile = datamodel.project + '/' + datamodel.project + \
            '_seds.grid.hd5'
        modelsedgrid = FileSEDGrid(modelsedgridfile)

        # read in the noise model just created
        noisemodel_vals = noisemodel.get_noisemodelcat(datamodel.noisefile)

        # trim the model sedgrid
        sed_trimname = '{0}/{0}_sed_trim.grid.hd5'.format(datamodel.project)
        noisemodel_trimname = '{0}/{0}_noisemodel_trim.hd5'.format(datamodel.project)

        trim_grid.trim_models(modelsedgrid, noisemodel_vals, obsdata,
                              sed_trimname, noisemodel_trimname, sigma_fac=3.)

    if args.fit:
        
        run_fitting.run_fitting(use_sd=False, nsubs=1, nprocs=1)


    if args.resume:

        run_fitting.run_fitting(use_sd=False, nsubs=1, nprocs=1,
                                    resume=True)


    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()
