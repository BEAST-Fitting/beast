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
from pipeline_production import run_fit, make_models, compute_noise_and_trim_grid
#import datamodel_small as datamodel
import datamodel_production as datamodel
import noisemodel 
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
    parser.add_argument("brick", help="brick number")
    parser.add_argument("source_density", help="source density bin")
    parser.add_argument("sub_source_density", help="subset of the source density bin [A, B, C, ...]")
    args = parser.parse_args()

    # update project datamodel information for this brick and source density bin
    datamodel.project += '_b' + args.brick
    datamodel.obsfile = datamodel.project + '/b' + args.brick + '_obs' + \
                        '/b' + args.brick + '_' + args.source_density + \
                        '/b' + args.brick + '_4band_det_' + args.source_density + '_' + args.sub_source_density + '.fits'
    datamodel.astfile = datamodel.project + '/merged_asts_b' + args.brick + '/fake_stars_b' + args.brick + \
                        '_' + args.source_density + '_all.hd5'
    datamodel.noisefile = datamodel.project + '/' + datamodel.project + \
                          '_sd' + args.source_density + '_noisemodel.hd5'
    stats_filebase = datamodel.project + '/' + datamodel.project + \
                     '_sd' + args.source_density + '_' + args.sub_source_density 
    sed_trimname = stats_filebase + '_sed_trim.grid.hd5'
    noisemodel_trimname = stats_filebase + '_noisemodel_trim.hd5'

    print("***run information***")
    print("  project = " + datamodel.project)
    print("  obsfile = " + datamodel.obsfile)
    print("  astfile = " + datamodel.astfile)
    print("        noisefile = " + datamodel.noisefile)
    print("   trimed sedfile = " + sed_trimname)
    print("trimed noisefiles = " + noisemodel_trimname)
    print("   stats filebase = " + stats_filebase)

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

        # read in the observed data
        obsdata = datamodel.get_obscat(datamodel.obsfile, datamodel.distanceModulus, datamodel.filters)
        # trim the model sedgrid

        trim_grid.trim_models(modelsedgrid, noisemodel_vals, obsdata, sed_trimname, noisemodel_trimname, sigma_fac=3.)

    if args.fit:
        start_time = time.clock()

        # the files for the trimmed model grid and noisemodel grid
        # read in the the AST noise model
        noisemodel_vals = noisemodel.get_noisemodelcat(noisemodel_trimname)

        run_fit(datamodel.project, sed_trimname, noise=noisemodel_vals, obsfile=datamodel.obsfile,
                outname=stats_filebase)

        new_time = time.clock()
        print('time to fit: ',(new_time - start_time)/60., ' min')
