#!/usr/bin/env python
"""
Running the BEAST on PHAT data - Trunchen style!
Karl G. - 18 Dec 2015

"""

# system imports
from __future__ import print_function
import sys
import argparse
import time

# BEAST imports
from pipeline_small import make_models
import datamodel_small as datamodel
import noisemodel 
from beast.fitting import fit
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
    parser.add_argument("-f", "--fit", help="Fit the data with the model " +
                        "and noise grids",
                        action="store_true")
    parser.add_argument("-r", "--resume", help="Resume a run",
                        action="store_true")
    parser.add_argument("--nocovar", help="Do not use the full covariance " + 
                        "noisemodel (i.e., just use diagonals)",
                        action="store_true")
    parser.add_argument("--sp1", help="Special case (brick 09, field 02)",
                        action="store_true")
    parser.add_argument("--sp2", help="Special case (brick 03, field 15)",
                        action="store_true")
    args = parser.parse_args()

    # different sets of models (easier than making new scripts)
    if args.sp1:
        datamodel.project = 'b09_f02_jan16_small'
        datamodel.obsfile = 'data/brick09_field02_n1000.fits'
        datamodel.astfile = 'ASTs_6band/12057_M31-B09-F02.gst.fake.fits'
        datamodel.noisefile = datamodel.project + '/' + datamodel.project + \
                              '_noisemodel.hd5'
    elif args.sp2:
        datamodel.project = 'b03_f15_jan16_small'
        datamodel.obsfile = 'data/brick03_field15_n1000.fits'
        datamodel.astfile = 'ASTs_6band/12109_M31-B03-F15.gst.fake.fits'
        datamodel.noisefile = datamodel.project + '/' + datamodel.project + \
                              '_noisemodel.hd5'

    if args.models:
        make_models()

    if args.noise:
        print('Generating noise model from ASTs and absflux A matrix')
 
        # get the modesedgrid on which to generate the noisemodel  
        modelsedgrid = FileSEDGrid( \
            '{project:s}/{project:s}_seds.grid.hd5'.format( \
                project=datamodel.project))  
            
        # generate the AST noise model  
        noisemodel.make_trunchen_noise_model(datamodel.noisefile,
                                             datamodel.astfile,
                                             datamodel.ast_colnames,
                                             modelsedgrid,
                generic_absflux_a_matrix=datamodel.generic_absflux_a_matrix)  

    if args.trim:
        print('Trimming the model and noise grids')
        start_time = time.clock()

        # get the modesedgrid on which to generate the noisemodel  
        modelsedgrid = FileSEDGrid('{project:s}/{project:s}_seds.grid.hd5'.\
                                   format(project=datamodel.project))
           

        # read in the noise model just created
        noisemodel_vals = noisemodel.get_noisemodelcat(datamodel.noisefile)

        # read in the observed data
        obsdata = datamodel.get_obscat(datamodel.obsfile, 
                                       datamodel.distanceModulus,
                                       datamodel.filters)

        # trim the model sedgrid
        sed_trimname = '{project:s}/{project:s}_seds_trim.grid.hd5'.\
                       format(project=datamodel.project)
        nm_trimname = '{project:s}/{project:s}_noisemodel_trim.grid.hd5'.\
                      format(project=datamodel.project)

        trim_grid.trim_models(modelsedgrid, noisemodel_vals, obsdata,
                              sed_trimname, nm_trimname, sigma_fac=3.,
                              trunchen=True)

        new_time = time.clock()
        print('time to trim: ',(new_time - start_time)/60., ' min')

    if args.fit:
        start_time = time.clock()
    
        # the files for the trimmed model grid and noisemodel grid
        modelsedgrid = '{project:s}/{project:s}_seds_trim.grid.hd5'.\
                       format(project=datamodel.project)
        noisemodelfile = '{project:s}/{project:s}_noisemodel_trim.grid.hd5'.\
                         format(project=datamodel.project)
        if args.nocovar:
            print("Running with simple noise model")
            datamodel.use_full_cov_matrix = False
            statsfile = '{project:s}/{project:s}_memory_nocovar_stats.fits'.\
                format(project=datamodel.project)
        else:
            print("Running with full noise model")
            statsfile = '{project:s}/{project:s}_memory_stats.fits'.\
                format(project=datamodel.project)

        # read in the the AST noise model
        noisemodel_vals = noisemodel.get_noisemodelcat(noisemodelfile)

        # read in the observed data
        obsdata = datamodel.get_obscat(datamodel.obsfile,
                                       datamodel.distanceModulus,
                                       datamodel.filters)

        fit.summary_table_memory(obsdata, noisemodel_vals,
                                 modelsedgrid, resume=args.resume,
                                 threshold=-10., save_every_npts=100,
                                 lnp_npts=60,
                                 stats_outname=statsfile,
                                 pdf1d_outname=statsfile.replace('stats.fits',
                                                                 'pdf1d.fits'),
                                 lnp_outname=statsfile.replace('stats.fits',
                                                               'lnp.hd5'),
                                 use_full_cov_matrix=\
                                 datamodel.use_full_cov_matrix)

        new_time = time.clock()
        print('time to fit: ',(new_time - start_time)/60., ' min')
