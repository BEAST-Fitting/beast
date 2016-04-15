#!/usr/bin/env python
"""
 Program to generatge input SEDs for multi-band ASTs
 Uses the BEAST model SED grid as the basis


.. history::
    Written 1 Feb 2016 by Karl D. Gordon
"""

from __future__ import print_function

import argparse

import numpy as np

from astropy.table import Table

import datamodel_small as datamodel
import noisemodel 
from beast.core.grid import FileSEDGrid
from beast.core.vega import Vega

from beast.tools.pbar import Pbar

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("filebase", type=str,
                        help='filename base')
    parser.add_argument("-a", "--nasts", type=int, default=100,
                        help="number of models to use")
    parser.add_argument("-p", "--npermodel", type=int, default=20,
                        help="number of models to use")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="verbose output")
    args = parser.parse_args()
    filebase = args.filebase

    print_diag = args.verbose

    # get the modesedgrid on which to generate the noisemodel  
    modelsedgrid = FileSEDGrid(filebase+'_seds.grid.hd5')

    # get the models to sample
    flux = modelsedgrid.seds
    n_models, n_filters = flux.shape

    # get the vega fluxes for zero magnitude
    with Vega() as v:
        _, vega_fluxes, _ = v.getFlux(datamodel.filters)

    # convert the datamodel sensitivity limits in mags to fluxes
    sens_limits = vega_fluxes* \
        np.power(10.0,-0.4*np.array(datamodel.sens_limits_mag))
    # same for bright limits
    bright_limits = vega_fluxes* \
        np.power(10.0,-0.4*np.array(datamodel.bright_limits_mag))
    # make 5x lower to make sure to capture models that noise may make
    #  detectable
    sens_limits *= 0.2 

    # determine the number of models that are above the sensitivity limits
    #   in at least one band
    ngood = np.zeros((n_models))
    nbright = np.zeros((n_models))
    for i in range(n_filters):
        gindxs, = np.where(flux[:,i] > sens_limits[i])
        ngood[gindxs] += 1

        bindxs, = np.where(flux[:,i] < bright_limits[i])
        nbright[bindxs] += 1
        #print(i, len(bindxs), min(flux[:,i]), max(flux[:,i]), bright_limits[i])

    gm_indxs, = np.where((ngood >= 1) & (nbright >= 3))
    n_gmodels = len(gm_indxs)

    # sort order for fluxes
    #   particular order required for the AST code
    sort_band_indxs = np.array([3,4,5,6,0,1,2,7,8])
    brightness_sort_indx = 6

    # setup table for the AST inputs
    tnames = np.array([x.lower() for x in datamodel.basefilters])
    tnames = ['d1','d2','x','y'] + tnames[sort_band_indxs].tolist()
    print(tnames)
    ast_inputs = Table(names=tnames,
                       dtype=[int, int, float, float] + 
                       [float]*n_filters)

    # pick the models for ASTs
    # use a sorted list on one bad and than pick every nth source
    # n is set by the total number desired
    n_asts = int(args.nasts)
    n_asts_per_model = int(args.npermodel)
    n_ast_models = n_asts/n_asts_per_model
    n = n_gmodels/n_ast_models
    indxs = np.arange(0,n_gmodels,n)
    # sort indxes, hard coded for F814W
    sindxs = np.argsort(flux[gm_indxs,brightness_sort_indx])

    # final indexs to the full flux array 
    #    (reversed to put the bright sources at the top)
    fin_indxs = gm_indxs[sindxs[indxs]][::-1]

    # min/max for x, y coordinates
    # hard coded for now, need a better soluiton later (based on source density)
    x_range = [5600, 8800]
    y_range = [6100, 8150]

    # loop over the model SEDs and generate the copies of the AST inputs
    #   with random x,y positions within the limits
    for findx in Pbar(desc='Evaluating model').iterover(fin_indxs):
        #for findx in fin_indxs:
        #print(findx)
        x_vals = x_range[0] + \
            np.random.random_sample(n_asts_per_model)*(x_range[1]-x_range[0])
        y_vals = y_range[0] + \
            np.random.random_sample(n_asts_per_model)*(y_range[1]-y_range[0])
        
        for j in range(n_asts_per_model):
            ast_mags = -2.5*np.log10(flux[findx,:]/vega_fluxes)
            ast_mags = ast_mags[sort_band_indxs]

            row_data = [0,1] + [x_vals[j]] + [y_vals[j]] + ast_mags.tolist()
            ast_inputs.add_row(row_data)

    # save the table of simulated observations
    ast_inputs.write(filebase+'_ast_inputs.dat', 
                     format='ascii.commented_header')
    ast_inputs.write(filebase+'_ast_inputs.fits', overwrite=True)
