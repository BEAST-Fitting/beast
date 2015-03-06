#!/usr/bin/env python
"""
Example of running the BEAST on M31 PHAT data.
Karl G. - 5 Nov 2014

------ previous comments by Morgan
Everything I need to generate a grid, make fake data and fit

I use the pipeline package I wrote in order to clean the syntax and allow more
flexibilities. In particular it will simplifies the managment of intermediate
results or broken jobs.

> python run [--models] [-?]

    --models         generates the models only (needed once)
    -?, --help       display this message

Make models
-----------
:func:`make_models` is equivalent to using individual tasks as follow:

project, noisefile, grid = project | t_project_dir
                                   | t_isochrones(**iso_kwargs)
                                   | t_spectra(**spec_kwargs)
                                   | t_seds(filters, **seds_kwargs)
                                   | t_gen_noise_model(astfile, **noise_kwargs)

Running the fit
---------------
:func:`run_fit`  is equivalent to using individual tasks as follow:

project, stat, obs, sedgrid = project | t_project_dir
                                      | t_get_obscat(**obscat_kwargs)
                                      | t_fit(g, noise, **fit_kwargs)
                                      | t_summary_table(g, **stat_kwargs)
"""

# system imports
from __future__ import print_function
import sys
import argparse
import time

# BEAST imports
from pipeline_small import run_fit, make_models
import datamodel_small as datamodel
import noisemodel 
import fit_memory
from merge_phat_asts import merge_phat_asts
from beast.core import prior_weights
from beast.core import trim_grid
from beast.core.grid import FileSEDGrid  

if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--models", help="Generate the model grid",
                        action="store_true")
    parser.add_argument("-n", "--noise", help="Calculate the noise model",
                        action="store_true")
    parser.add_argument("-t", "--trim", help="Trim the model and noise grids",
                        action="store_true")
    parser.add_argument("-f", "--fit", help="Fit the observed data",
                        action="store_true")
    args = parser.parse_args()

    if args.models:
        make_models()

    if args.noise:
        print('Generating noise model from ASTs and absflux A matrix')
 
        # get the modesedgrid on which to generate the noisemodel  
        modelsedgrid = FileSEDGrid('{project:s}/{project:s}_seds.grid.hd5'.format(project=datamodel.project))  
            
        # generate the AST noise model  
        noisemodel.make_toothpick_noise_model(datamodel.noisefile, datamodel.astfile, modelsedgrid, datamodel.absflux_a_matrix)  

    if args.trim:
        print('Trimming the model and noise grids')

        # get the modesedgrid on which to generate the noisemodel  
        modelsedgrid = FileSEDGrid('{project:s}/{project:s}_seds.grid.hd5'.format(project=datamodel.project))  

        # read in the noise model just created
        noisemodel_vals = noisemodel.get_noisemodelcat(datamodel.noisefile)

        # read in the ast file used to create the noise model
        astdata = noisemodel.PHAT_ToothPick_Noisemodel(datamodel.astfile, modelsedgrid.filters)            

        # read in the observed data
        obsdata = datamodel.get_obscat(datamodel.obsfile, datamodel.distanceModulus, datamodel.filters)

        # trim the model sedgrid
        sed_trimname = '{project:s}/{project:s}_seds_trim.grid.hd5'.format(project=datamodel.project)
        noisemodel_trimname = '{project:s}/{project:s}_noisemodel_trim.grid.hd5'.format(project=datamodel.project)

        trim_grid.trim_models(modelsedgrid, noisemodel_vals, obsdata, astdata, sed_trimname, noisemodel_trimname, sigma_fac=3.)

    if args.fit:
        start_time = time.clock()
    
        # the files for the trimmed model grid and noisemodel grid
        modelsedgrid = '{project:s}/{project:s}_seds_trim.grid.hd5'.format(project=datamodel.project)
        noisemodelfile = '{project:s}/{project:s}_noisemodel_trim.grid.hd5'.format(project=datamodel.project)
        statsfile = '{project:s}/{project:s}_memory_stats.fits'.format(project=datamodel.project)

        # read in the the AST noise model
        noisemodel_vals = noisemodel.get_noisemodelcat(noisemodelfile)

        # read in the observed data
        obsdata = datamodel.get_obscat(datamodel.obsfile, datamodel.distanceModulus, datamodel.filters)

        fit_memory.summary_table_memory(obsdata, noisemodel_vals, modelsedgrid,
                                        outname=statsfile)

        new_time = time.clock()
        print('time to fit: ',(new_time - start_time)/60., ' min')
