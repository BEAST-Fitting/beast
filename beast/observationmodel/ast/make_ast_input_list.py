from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys, argparse
import numpy as np
from ..vega import Vega
import datamodel as datamodel
from ...physicsmodel.grid import FileSEDGrid
from astropy.table import Table
from astropy.io import ascii



def mag_limits(seds,limits,Nfilter=1):
    """
    Selects models which have at least N filter above the limits

    INPUTS:
    -------
    seds:   np.array
            Magnitude array from BEAST grid 

    limits: list
            List of limit magnitudes 

    Nfilter: integer
             In how many filters, you want a fake star to be brighter than the limit
 
    OUTPUT:
    -------
    idx:    np.array
            Array of integers contining the indices of allowed models

    """
    flag = seds.copy() 

    # flag is True if the models are brigter (=smaller number in mag) than the limits
    for i,limit in enumerate(limits):
        flag[:,i] = seds[:,i] < limit 
    
    # Keep index where model is brighter than the limit in N filters
    s = np.sum(flag,axis=1)
    idx, = np.where(s >= Nfilter)            
        
    return idx


def pick_models(sedgrid, mag_cuts, Nfilter=3, N_stars= 70, Nrealize = 20):
    """
    Creates a fake star catalog from a BEAST model grid

    INPUTS:
    -------
    sedgrid: beast.grid
               BEAST model grid from which the models are picked

    mag_cuts: list
               List of magnitude limits 

    Nfilter: Integer
             In how many filters, you want a fake star to be brighter than the limit
             (default = 3)

    N_stars: Integer
               Number of stellar models picked per a single log(age) (default=70)

    Nrealize: Integer
              Number of realization of each models (default = 20)
               
    return:
    -------
    ascii file: A list of selected models
    """

    filters = datamodel.filters

    with Vega() as v:               # Get the vega fluxes
        vega_f,vega_flux,lamb = v.getFlux(filters)

    sedsMags = -2.5 * np.log10(sedgrid.seds[:]/vega_flux)  # Convert to Vega mags  

    # Select the models above the magnitude limits in N filters
    idx = mag_limits(sedsMags, mag_cuts, Nfilter=Nfilter) 
    grid_cut = sedgrid.grid[idx] 
    mod_tot = len(idx)

    # Sample the model grid uniformly
    prime_params = np.column_stack((grid_cut['logA'],grid_cut['M_ini'],grid_cut['Av']))
    search_age = np.unique(prime_params[:,0])

    N_sample = N_stars
    models = []
    for iage in search_age:
        tmp, = np.where(prime_params[:,0] == iage)
        models.append(np.random.choice(tmp,N_sample))

         
    index = np.repeat(idx[np.array(models).reshape((-1))],Nrealize)      
    sedsMags = Table(sedsMags[index,:], names=filters)

    outfile = './' + datamodel.project + '/' + datamodel.project + '_inputAST.txt'

    ascii.write(sedsMags, outfile, overwrite=True, formats={k:'%.5f' for k in sedsMags.colnames}) 
