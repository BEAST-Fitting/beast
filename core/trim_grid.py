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

def trim_models(sedgrid, sedgrid_noisemodel, obsdata, astdata, sed_outname, noisemodel_outname,
                sigma_fac=3., n_detected=4):
    
    # Store the brigtest and faintest fluxes in each band (for data and asts)
    n_filters = len(obsdata.filters)
    min_data = np.zeros(n_filters)
    max_data = np.zeros(n_filters)
    min_asts = np.zeros(n_filters)
    max_asts = np.zeros(n_filters)
    min_models = np.zeros(n_filters)
    max_models = np.zeros(n_filters)
    for k, filtername in enumerate(obsdata.filters):
        # get the name of the column with the rate in it (normalized to the vega flux)
        nfiltername = filtername.split('_')[-1].lower() + '_rate'
        min_data[k] = np.amin(obsdata.data[:][nfiltername]*obsdata.vega_flux[k])
        max_data[k] = np.amax(obsdata.data[:][nfiltername]*obsdata.vega_flux[k])
        afiltername = filtername.split('_')[-1].upper() + '_IN'
        ofiltername = filtername.split('_')[-1].upper() + '_VEGA'
        indxs, = np.where(astdata.data[ofiltername] < 80.)
        _ast_fluxes = 10 ** (-0.4*astdata.data[afiltername][indxs])*obsdata.vega_flux[k]
        min_asts[k] = np.amin(_ast_fluxes)
        max_asts[k] = np.amax(_ast_fluxes)
        min_models[k] = np.amin(sedgrid.seds[:,k])
        max_models[k] = np.amax(sedgrid.seds[:,k])

    #print(min_data)
    #print(min_models)
    #print(min_asts)
    #exit()

    # cache the model bias and uncertainty
    model_bias = sedgrid_noisemodel.root.bias[:]
    model_unc = sedgrid_noisemodel.root.error[:]
    model_compl = sedgrid_noisemodel.root.completeness[:]

    # Get upper and lower values for the models given the noise model 
    #  sigma_fac defaults to 3.
    model_down = sedgrid.seds + model_bias - sigma_fac*model_unc
    model_up = sedgrid.seds + model_bias + sigma_fac*model_unc
    
    # first remove all models that have any band with fluxes below the faintest ASTs run
    #indxs = np.arange(len(sedgrid.seds[:,0]))
    gvals = np.zeros(len(sedgrid.seds[:,0]),np.int64)
    for i in range(len(sedgrid.seds[:,0])):
        # only look at models above the faintest asts run to avoide huge extrapolations
        #  extrapolations to higher fluxes is likely ok given the small values at high fluxes
        gindxs, = np.where((sedgrid.seds[i,:] >= min_asts) == True)
        if len(gindxs) >= n_detected:
            gvals[i] = 1

    indxs, = np.where(gvals > 0)
    
    if len(indxs) <= 0:
        print('no models are brighter than the minimum ASTs run')
        exit()

    n_ast_indxs = len(indxs)

    # Find models with fluxes (with margin) between faintest and brightest data
    for k in range(n_filters):
        nindxs, = np.where((model_up[indxs,k] >= min_data[k]) & (model_down[indxs,k] <= max_data[k]))
        if len(nindxs) > 0:
            indxs = indxs[nindxs]

    if len(indxs) == 0:
        print('no models that are within the data range')
        exit()

    print('number of original models = ', len(sedgrid.seds[:,0]))
    print('number of ast trimmed models = ', n_ast_indxs)
    print(' number of trimmed models = ', len(indxs))

    #print(min_data)
    #print(max_data)
    #for k in range(n_filters):
    #    print(np.min(sedgrid.seds[indxs,k]),np.max(sedgrid.seds[indxs,k]), np.min(model_down[indxs,k]), np.max(model_up[indxs,k]))
    #exit()

    #Save the grid
    print('Writing trimmed sedgrid to disk into {0:s}'.format(sed_outname))
    cols = {}
    for key in sedgrid.grid.keys():
        cols[key] = sedgrid.grid[key][indxs]

    cols['fullgrid_idx'] = indxs.astype(int) # New column to save the index of the model in the full grid
    g = SpectralGrid(sedgrid.lamb, seds=sedgrid.seds[indxs], grid=Table(cols), backend='memory')
    #filter_names = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W', 'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    #g.grid.header['filters'] = ' '.join(filter_names)
    filternames = obsdata.filters
    g.grid.header['filters'] = ' '.join(filternames)
    
    g.writeHDF(sed_outname) # trimmed grid name

    # save the trimmed noise model
    print('Writing trimmed noisemodel to disk into {0:s}'.format(noisemodel_outname))
    with tables.openFile(noisemodel_outname, 'w') as outfile:
        outfile.createArray(outfile.root,'bias', model_bias[indxs])
        outfile.createArray(outfile.root,'error', model_unc[indxs])
        outfile.createArray(outfile.root,'completeness', model_compl[indxs])


