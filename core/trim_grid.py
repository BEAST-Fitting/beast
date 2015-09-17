"""
Trim the grid of models
================

For a given set of observations, there will be models that are so
bright or faint that they will always have ~0 probability of fitting
the data.  This program trims those models out of the SED grid
so that time is not spent calculating model points that are always 
zero probability.
"""
from __future__ import print_function

import numpy as np
import tables

from .grid import FileSEDGrid
from .grid import SpectralGrid
from ..external.eztables import Table

def trim_models(sedgrid, sedgrid_noisemodel, obsdata, sed_outname, noisemodel_outname,
                sigma_fac=3., n_detected=4, inFlux=True):
    
    # Store the brigtest and faintest fluxes in each band (for data and asts)
    n_filters = len(obsdata.filters)
    min_data = np.zeros(n_filters)
    max_data = np.zeros(n_filters)
    min_models = np.zeros(n_filters)
    max_models = np.zeros(n_filters)
    for k, filtername in enumerate(obsdata.filters):
        if inFlux:
            min_data[k] = np.amin(obsdata.data[:][obsdata.data.resolve_alias(filtername)]*obsdata.vega_flux[k])
            max_data[k] = np.amax(obsdata.data[:][obsdata.data.resolve_alias(filtername)]*obsdata.vega_flux[k])
        else:
            min_data[k] = np.amin(10 **(-0.4*obsdata.data[:][obsdata.data.resolve_alias(filtername)])*obsdata.vega_flux[k])
            max_data[k] = np.amax(10 **(-0.4*obsdata.data[:][obsdata.data.resolve_alias(filtername)])*obsdata.vega_flux[k])


        min_models[k] = np.amin(sedgrid.seds[:,k])
        max_models[k] = np.amax(sedgrid.seds[:,k])

    # first remove all models that have any band with fluxes below the faintest ASTs run
    # when the noisemodel was computed, models with fluxes below the faintest ASTs were tagged with a negative error/uncertainty
    # identify the models that have been detected in enough bands
    #   the idea here is that if the ASTs are not measured that means that *none* were recovered and this implies
    #   that no model with these values would be recovered and thus the probability should always be zero
    model_unc = sedgrid_noisemodel.root.error[:]
    above_ast = (model_unc > 0)
    sum_above_ast = np.sum(above_ast,axis=1)
    indxs, = np.where(sum_above_ast >= n_detected)

    #min_gmodel = np.zeros(n_filters)
    #for k in range(n_filters):
    #    min_gmodel[k] = np.amin(sedgrid.seds[indxs,k])


    #min_asts = [1.3879187495103194e-19,
    #            1.6550853881968002e-19,
    #            1.068305101054377e-20,
    #            2.0358040415431349e-20,
    #            1.0525632892236783e-19,
    #            7.2018979602833484e-20]

    #min_asts = [2.7215185446506946e-19,
    #            2.1813920602921763e-19,
    #            1.420197564598432e-19,
    #            6.6257679757877467e-20,
    #            7.3476470356322374e-19,
    #            1.9272005768395571e-19]

    #print(min_asts)
    #print(min_gmodel)

    #gvals = np.zeros(len(sedgrid.seds[:,0]),np.int64)  
    #for i in range(len(sedgrid.seds[:,0])):  
        # only look at models above the faintest asts run to avoide huge extrapolations  
        #  extrapolations to higher fluxes is likely ok given the small values at high fluxes  
    #    gindxs, = np.where((sedgrid.seds[i,:] >= min_asts) == True)  
    #    if len(gindxs) >= n_detected:  
    #        gvals[i] = 1  
    
    #oindxs, = np.where(gvals > 0)  

    #print(len(oindxs), len(indxs))

    #print(oindxs)
    #print(indxs)

    #exit()

    # cache the noisemodel values
    model_bias = sedgrid_noisemodel.root.bias[:]
    model_unc = np.fabs(sedgrid_noisemodel.root.error[:])  # needed to remove the negative values that designate a model below the faintest ASTs
    model_compl = sedgrid_noisemodel.root.completeness[:]

    if len(indxs) <= 0:
        print('no models are brighter than the minimum ASTs run')
        exit()

    #indxs = np.arange(len(model_bias))
    n_ast_indxs = len(indxs)

    # Find models with fluxes (with margin) between faintest and brightest data
    for k in range(n_filters):
        print('working on filter # = ', k)

        #min_gmodel = np.amin(sedgrid.seds[indxs,k])
        #print(min_gmodel, min_data[k], min_models[k])
        #max_gmodel = np.amax(sedgrid.seds[indxs,k])
        #print(max_gmodel, max_data[k], max_models[k])

        # Get upper and lower values for the models given the noise model 
        #  sigma_fac defaults to 3.
        model_val = sedgrid.seds[indxs,k] + model_bias[indxs,k]
        model_down = model_val - sigma_fac*model_unc[indxs,k]
        model_up = model_val + sigma_fac*model_unc[indxs,k]

        nindxs, = np.where((model_up >= min_data[k]) & (model_down <= max_data[k]))
        if len(nindxs) > 0:
            indxs = indxs[nindxs]

    if len(indxs) == 0:
        print('no models that are within the data range')
        exit()

    print('number of original models = ', len(sedgrid.seds[:,0]))
    print('number of ast trimmed models = ', n_ast_indxs)
    print(' number of trimmed models = ', len(indxs))

    #Save the grid
    print('Writing trimmed sedgrid to disk into {0:s}'.format(sed_outname))
    cols = {}
    for key in sedgrid.grid.keys():
        cols[key] = sedgrid.grid[key][indxs]

    cols['fullgrid_idx'] = indxs.astype(int) # New column to save the index of the model in the full grid
    g = SpectralGrid(sedgrid.lamb, seds=sedgrid.seds[indxs], grid=Table(cols), backend='memory')
    filternames = obsdata.filters
    g.grid.header['filters'] = ' '.join(filternames)

    g.writeHDF(sed_outname) # trimmed grid name

    # save the trimmed noise model
    print('Writing trimmed noisemodel to disk into {0:s}'.format(noisemodel_outname))
    with tables.openFile(noisemodel_outname, 'w') as outfile:
        outfile.createArray(outfile.root,'bias', model_bias[indxs])
        outfile.createArray(outfile.root,'error', model_unc[indxs])
        outfile.createArray(outfile.root,'completeness', model_compl[indxs])


