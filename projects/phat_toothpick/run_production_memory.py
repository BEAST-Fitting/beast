#!/usr/bin/env python
"""
Code to run the BEAST on the PHAT data in production mode.
  - specific to the BEAST toothpick version (v1, 5 Feb 2015)
  - specific to the format setup by Heddy A. for the PHAT bricks and source density regions
  - setup to allow for batch processing through commandline variables

based on previous code by Morgan F.
started modifications by Karl G. - 5 Nov 2014
finalized by Karl G. - 5 Feb 2015

see datamodel_production.py, pipeline_production.py, noisemodel.py, fit.py for more details
"""

# system imports
from __future__ import print_function
import sys
import os
import argparse
import time
import string

# BEAST imports
from pipeline_production import run_fit, make_models, compute_noise_and_trim_grid
#import datamodel_small as datamodel
import datamodel_production as datamodel
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
    parser.add_argument("-r", "--resume", help="Resume a run",
                        action="store_true")
    parser.add_argument("-a", "--faint", help="Faint set of 4 band detections (instead of bright)",
                        action="store_true")
    parser.add_argument("-i", "--index", help="Use indexes to generate sed and noisemodel trim files",
                        action="store_true")
    parser.add_argument("brick", help="brick number")
    parser.add_argument("source_density", help="source density bin")
    parser.add_argument("sub_source_density", help="subset of the source density bin [A, B, C, ...]")
    args = parser.parse_args()

    # update project datamodel information for this brick and source density bin
    datamodel.project = 'b' + args.brick
    if args.faint:
        datamodel.project += 'f'
        datamodel.obsfile = 'BEAST_production/' + datamodel.project + '/obscat/b' + args.brick + \
            '-6filt-cut-4band-gst-faint-SD-' + string.replace(args.source_density,'_','-') + \
            '-sub' + args.sub_source_density + '.fits'
    else:
        datamodel.obsfile = 'BEAST_production/' + datamodel.project + '/obscat/b' + args.brick + \
            '-6filt-cut-4band-gst-bright-SD-' + string.replace(args.source_density,'_','-') + \
            '-sub' + args.sub_source_density + '.fits'

    datamodel.astfile = 'BEAST_production/merged_asts/PHAT_fake_stars_SD_' + \
                        string.replace(args.source_density,'-','_' ) + '.fits'
    datamodel.noisefile = 'BEAST_production/BEAST_production_sd_' + \
                          string.replace(args.source_density,'_','-' ) + '_noisemodel.fits'

    stats_filebase = 'BEAST_production/' + datamodel.project + '/b' + args.brick + \
                     '_sd' + string.replace(args.source_density,'_','-') + '_sub' + args.sub_source_density 
    sed_trimname = stats_filebase + '_sed_trim.grid.hd5'
    noisemodel_trimname = stats_filebase + '_noisemodel_trim.hd5'
    index_name = stats_filebase+'_sed_trim.grid_indexs.fits'

    #modelsedgrid_filename = 'b15_late_jan15_test_small/b15_late_jan15_test_small_seds.grid.hd5'
    modelsedgrid_filename = 'BEAST_production/BEAST_production_seds.grid.hd5'

    print("***run information***")
    print("  project = " + datamodel.project)
    print("  obsfile = " + datamodel.obsfile)
    print("  astfile = " + datamodel.astfile)
    print("        noisefile = " + datamodel.noisefile)
    print("   trimed sedfile = " + sed_trimname)
    print("trimed noisefiles = " + noisemodel_trimname)
    print("   stats filebase = " + stats_filebase)

    if args.models:
        datamodel.project = 'BEAST_production'
        make_models()

    if args.noise:
        print('Generating noise model from ASTs and absflux A matrix')
 
        # get the modesedgrid on which to generate the noisemodel  
        modelsedgrid = FileSEDGrid(modelsedgrid_filename)  

        # generate the AST noise model  
        noisemodel.make_toothpick_noise_model(datamodel.noisefile, datamodel.astfile, modelsedgrid, datamodel.absflux_a_matrix)  

    if args.trim:
        print('Trimming the model and noise grids')
        start_time = time.clock()

        # check if the trimmed files already exist
        if (not os.path.isfile(sed_trimname)) | (not os.path.isfile(noisemodel_trimname)):
            # read in the observed data
            obsdata = datamodel.get_obscat(datamodel.obsfile, datamodel.distanceModulus, datamodel.filters)

            # get the modesedgrid on which to generate the noisemodel  
            modelsedgrid = FileSEDGrid(modelsedgrid_filename.format(project=datamodel.project))  

            # read in the noise model
            noisemodel_vals = noisemodel.get_noisemodelcat(datamodel.noisefile)

            # trim the model sedgrid
            trim_grid.trim_models(modelsedgrid, noisemodel_vals, obsdata, sed_trimname, noisemodel_trimname, sigma_fac=3.)
        else:
            print('trimming requested - but trimmed sed and noisemodel files already exist')
            print('using existing trimmed files')

        new_time = time.clock()
        print('time to trim: ',(new_time - start_time)/60., ' min')

    if args.fit:
        start_time = time.clock()

        if args.index:

        else:
            # the files for the trimmed model grid and noisemodel grid
            # read in the the AST noise model
            noisemodel_vals = noisemodel.get_noisemodelcat(noisemodel_trimname)

        # read in the observed data
        obsdata = datamodel.get_obscat(datamodel.obsfile, datamodel.distanceModulus, datamodel.filters)

        statsfile = stats_filebase + '_stats.fits'
        fit_memory.summary_table_memory(obsdata, noisemodel_vals, sed_trimname, resume=args.resume,
                                        threshold=-10., save_every_npts=250, lnp_npts=60,
                                        stats_outname=stats_filebase + '_stats.fits',
                                        pdf1d_outname=stats_filebase + '_pdf1d.fits',
                                        lnp_outname=stats_filebase + '_lnp.hd5')

        new_time = time.clock()
        print('time to fit: ',(new_time - start_time)/60., ' min')
