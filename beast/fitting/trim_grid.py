"""
Trim the grid of models
=======================

For a given set of observations, there will be models that are so
bright or faint that they will always have ~0 probability of fitting
the data.  This program trims those models out of the SED grid
so that time is not spent calculating model points that are always 
zero probability.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import tables

from ..physicsmodel.grid import FileSEDGrid
from ..physicsmodel.grid import SpectralGrid
from ..external.eztables import Table

__all__ = ['trim_models']

def trim_models(sedgrid, sedgrid_noisemodel, obsdata, sed_outname,
                noisemodel_outname,
                sigma_fac=3., n_detected=4, inFlux=True,
                trunchen=False):
    """

    Keywords
    ---------
    sedgrid: str or grid.SEDgrid instance
    	model grid

    sedgrid_noisemodel: beast noisemodel instance
    	noise model data

    obsdata: Observation object instance
    	observation catalog

    sed_outname: str
	name for output sed file

    noisemodel_outname: str
    	name for output noisemodel file

    sigma_fac: float
    	factor for trimming the upper and lower range of grid so that
	the model range cuts off sigma_fac above and below the brightest
	and faintest models, respectively (default: 3.)

    n_detected: int
    	minimum number of bands where ASTs yielded a detection for
	a given model, if fewer detections than n_detected this model
	gets eliminated (default: 4)
    		
    inFlux: boolean
	if true data are in fluxes (default: True) 

    trunchen: boolean
    	if true use the trunchen noise model (default: False)


    returns
    --------
    N/A
    """
    
    # Store the brigtest and faintest fluxes in each band (for data and asts)
    n_filters = len(obsdata.filters)
    min_data = np.zeros(n_filters)
    max_data = np.zeros(n_filters)
    min_models = np.zeros(n_filters)
    max_models = np.zeros(n_filters)
    for k, filtername in enumerate(obsdata.filters):
        sfiltname = obsdata.data.resolve_alias(filtername)
        if inFlux:
            min_data[k] = np.amin(obsdata.data[sfiltname]*
                                  obsdata.vega_flux[k])
            max_data[k] = np.amax(obsdata.data[sfiltname]*
                                  obsdata.vega_flux[k])
        else:
            min_data[k] = np.amin(10 **(-0.4*obsdata.data[sfiltname])
                                  *obsdata.vega_flux[k])
            max_data[k] = np.amax(10 **(-0.4*obsdata.data[sfiltname])
                                  *obsdata.vega_flux[k])

        min_models[k] = np.amin(sedgrid.seds[:,k])
        max_models[k] = np.amax(sedgrid.seds[:,k])

    # first remove all models that have any band with fluxes below the
    #    faintest ASTs run
    # when the noisemodel was computed, models with fluxes below the
    #    faintest ASTs were tagged with a negative error/uncertainty
    # identify the models that have been detected in enough bands
    #   the idea here is that if the ASTs are not measured that means
    #   that *none* were recovered and this implies
    #   that no model with these values would be recovered and thus the
    #   probability should always be zero
    model_unc = sedgrid_noisemodel.root.error[:]
    above_ast = (model_unc > 0)
    sum_above_ast = np.sum(above_ast,axis=1)
    indxs, = np.where(sum_above_ast >= n_detected)

    # cache the noisemodel values
    model_bias = sedgrid_noisemodel.root.bias[:]
    model_unc = np.fabs(sedgrid_noisemodel.root.error[:])  
    model_compl = sedgrid_noisemodel.root.completeness[:]
    if trunchen:
        model_q_norm = sedgrid_noisemodel.root.q_norm[:]
        model_icov_diag = sedgrid_noisemodel.root.icov_diag[:]
        model_icov_offdiag = sedgrid_noisemodel.root.icov_offdiag[:]

    if len(indxs) <= 0:
        print('no models are brighter than the minimum ASTs run')
        exit()

    n_ast_indxs = len(indxs)

    # Find models with fluxes (with margin) between faintest and brightest data
    for k in range(n_filters):
        print('working on filter # = ', k)

        # Get upper and lower values for the models given the noise model 
        #  sigma_fac defaults to 3.
        model_val = sedgrid.seds[indxs,k] + model_bias[indxs,k]
        model_down = model_val - sigma_fac*model_unc[indxs,k]
        model_up = model_val + sigma_fac*model_unc[indxs,k]

        nindxs, = np.where((model_up >= min_data[k]) &
                           (model_down <= max_data[k]))
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
    for key in list(sedgrid.grid.keys()):
        cols[key] = sedgrid.grid[key][indxs]

    # New column to save the index of the model in the full grid
    cols['fullgrid_idx'] = indxs.astype(int) 
    g = SpectralGrid(sedgrid.lamb, seds=sedgrid.seds[indxs],
                     grid=Table(cols), backend='memory')
    filternames = obsdata.filters
    g.grid.header['filters'] = ' '.join(filternames)

    g.writeHDF(sed_outname) # trimmed grid name

    # save the trimmed noise model
    print('Writing trimmed noisemodel to disk into {0:s}'.\
          format(noisemodel_outname))
    with tables.open_file(noisemodel_outname, 'w') as outfile:
        outfile.create_array(outfile.root,'bias', model_bias[indxs])
        outfile.create_array(outfile.root,'error', model_unc[indxs])
        outfile.create_array(outfile.root,'completeness', model_compl[indxs])
        if trunchen:
            outfile.create_array(outfile.root,'q_norm', model_q_norm[indxs])
            outfile.create_array(outfile.root,'icov_diag',
                                model_icov_diag[indxs])
            outfile.create_array(outfile.root,'icov_offdiag',
                                model_icov_offdiag[indxs])
            


