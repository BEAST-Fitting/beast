from __future__ import (absolute_import, division, print_function)
import os

import numpy as np

from . import grid
from . import creategrid
from .stars import isochrone
from .stars.isochrone import ezIsoch
from .grid_and_prior_weights import compute_age_mass_metallicity_weights

from ..tools.helpers import val_in_unit

__all__ = ['make_iso_table', 'make_spectral_grid', 'add_stellar_priors']

def make_iso_table(project, oiso=None, logtmin=6.0, logtmax=10.13, dlogt=0.05,
                   z=[0.0152]):
    """
    The isochrone tables are loaded (downloading if necessary)

    Parameters
    ----------
    project: str
        project name

    oiso: isochrone.Isochrone object
        contains the full isochrones information

    logtmin: float
        log-age min

    logtmax: float
        log-age max

    dlogt: float
        log-age step to request

    z: float or sequence
        list of metalicity values, where default (Z=0.152) is adopted Z_sun
        for PARSEC/COLIBRI models

    Returns
    -------
    outname: str
        file into which save the table of isochrones (any format eztables
        can handle)
    """
    iso_fname = '%s/%s_iso.csv' % (project, project)
    if not os.path.isfile(iso_fname):
        if oiso is None: 
            oiso = isochrone.PadovaWeb()
        
        t = oiso._get_t_isochrones(max(5.0, logtmin), min(10.13, logtmax),
                                   dlogt, z)
        t.header['NAME'] = '{0} Isochrones'.format('_'.join(iso_fname.split('_')[:-1]))
        print('{0} Isochrones'.format('_'.join(iso_fname.split('_')[:-1])))
    
        t.write(iso_fname)

    # read in the isochrone data from the file
    #   not sure why this is needed, but reproduces previous ezpipe method
    oiso = ezIsoch(iso_fname)
        
    return (iso_fname, oiso)

def make_spectral_grid(project, oiso, osl=None, bounds={}, distance=None,
                       verbose=True,
                       add_spectral_properties_kwargs=None, **kwargs):
    """
    The spectral grid is generated using the stellar parameters by
    interpolation of the isochrones and the generation of spectra into the
    physical units

    Parameters
    ----------
    project: str
        project name

    oiso: isochrone.Isochrone object
        set of isochrones to use

    osl: stellib.Stellib object
        Spectral library to use (default stellib.Kurucz)

    distance: float
        Distance at which models should be shifted
        0 means absolute magnitude.
        Expecting pc units

    add_spectral_properties_kwargs: dict
        keyword arguments to call :func:`add_spectral_properties`
        to add model properties from the spectra into the grid property table

    Returns
    -------

    outname: str
        file into which save the spectral grid
    """
    spec_fname = '%s/%s_spec_grid.hd5' % (project, project)

    if not os.path.isfile(spec_fname):
        osl = osl or stellib.Kurucz()

        #filter extrapolations of the grid with given sensitivities in logg
        #  and logT
        if 'dlogT' not in bounds:
            bounds['dlogT'] = 0.1
        if 'dlogg' not in bounds:
            bounds['dlogg'] = 0.3

        #make the spectral grid
        if verbose:
            print('Make spectra')
        g = creategrid.gen_spectral_grid_from_stellib_given_points(osl,
                                                                   oiso.data,
                                                               bounds=bounds)

        # get the distance
        if distance is not None:
            _distance = val_in_unit('distance', distance, 'pc').magnitude

        if verbose:    
            print('Adding spectral properties:', add_spectral_properties_kwargs
                  is not None)
        if add_spectral_properties_kwargs is not None:
            nameformat = add_spectral_properties_kwargs.\
                         pop('nameformat', '{0:s}') + '_nd'

        # write to disk
        # and apply the distance to the particular galaxy of interest
        # seds already at 10 pc, need multiplcation by the square of the ratio
        # to this distance
        if hasattr(g, 'writeHDF'):
            if distance is not None:
                g.seds = g.seds / (0.1 * _distance) ** 2
            if add_spectral_properties_kwargs is not None:
                g = creategrid.add_spectral_properties(g,
                                                       nameformat=nameformat,
                                            **add_spectral_properties_kwargs)
            g.writeHDF(spec_fname)
        else:
            for gk in g:
                if distance is not None:
                    gk.seds = gk.seds / (0.1 * _distance) ** 2
                if add_spectral_properties_kwargs is not None:
                    gk = creategrid.add_spectral_properties(gk,
                                              nameformat=nameformat,
                                              **add_spectral_properties_kwargs)
    
                gk.writeHDF(spec_fname, append=True)
    else:
        g = grid.FileSpectralGrid(spec_fname, backend='memory')
        
    return (spec_fname, g)

def add_stellar_priors(project, specgrid, verbose=True, **kwargs):
    """
    make_priors -- compute the weights for the stellar priors

    Parameters
    ----------
    project: str
        project name

    specgrid: grid.SpectralGrid object
        spectral grid to transform
        result from the make_spectra function

    returns
    -------

    outname: str
        file into which save the SED grid
    """
    priors_fname = '%s/%s_spec_w_priors.grid.hd5' % (project, project)
    if not os.path.isfile(priors_fname):

        if verbose:
            print('Make Prior Weights')

        compute_age_mass_metallicity_weights(specgrid.grid)

        #write to disk
        if hasattr(specgrid, 'writeHDF'):
            specgrid.writeHDF(priors_fname)
        else:
            for gk in specgrid:
                gk.writeHDF(priors_fname, append=True)

    g = grid.FileSpectralGrid(priors_fname, backend='memory')

    return (priors_fname, g)
