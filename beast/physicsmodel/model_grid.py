from __future__ import (absolute_import, division, print_function)
import os

import h5py
import numpy as np
from astropy import units

from . import grid
from . import creategrid
from .stars import isochrone
from .stars.isochrone import ezIsoch
from .grid_and_prior_weights import compute_age_mass_metallicity_weights
from ..external.eztables import Table

__all__ = ['make_iso_table', 'make_spectral_grid', 'add_stellar_priors',
           'make_extinguished_sed_grid']

def make_iso_table(project, oiso=None, logtmin=6.0, logtmax=10.13, dlogt=0.05,
                   z=[0.0152], iso_fname=None):
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
    fname: str
       name of saved file

    oiso: isochrone.Isochrone object
        contains the full isochrones information
    """
    if iso_fname is None:
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

def make_spectral_grid(project, oiso, osl=None, bounds={},
                       verbose=True, spec_fname=None, distance=10,
                       distance_unit=units.pc, filterLib=None,
                       add_spectral_properties_kwargs=None,
                       num_dist_subgrids=None, nprocs=None, **kwargs):
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

    distance: float or list of float
        distances at which models should be shifted, specified as a
        single number or as [min, max, step]

        0 means absolute magnitude.

    distance_unit: astropy length unit or mag
        distances will be evenly spaced in this unit
        therefore, specifying a distance grid in mag units will lead to
        a log grid

    spec_fname: str
        full filename to save the spectral grid into

    filterLib:  str
        full filename to the filter library hd5 file

    add_spectral_properties_kwargs: dict
        keyword arguments to call :func:`add_spectral_properties`
        to add model properties from the spectra into the grid property table

    Returns
    -------
    fname: str
       name of saved file

    g: grid.SpectralGrid object
        spectral grid to transform
    """
    if spec_fname is None:
        spec_fname = '%s/%s_spec_grid.hd5' % (project, project)

    subgrid_names = []
    if num_dist_subgrids is not None:
        for i in range(num_dist_subgrids):
            subgrid_name = spec_fname.replace(
                ".hd5", "_dist_sub{}.hd5".format(i))
            subgrid_names.append(subgrid_name)
    else:
        subgrid_names = [spec_fname]

    missing = False
    for n in subgrid_names:
        missing = missing or not os.path.isfile(n)

    if missing:
        # Clear the existing files, and then create new ones
        for n in subgrid_names + [spec_fname]:
            try:
                os.remove(n)
            except:
                pass

        osl = osl or stellib.Kurucz()

        # filter extrapolations of the grid with given sensitivities in
        # logg and logT
        if 'dlogT' not in bounds:
            bounds['dlogT'] = 0.1
        if 'dlogg' not in bounds:
            bounds['dlogg'] = 0.3

        # make the spectral grid
        if verbose:
            print('Make spectra')
        g = creategrid.gen_spectral_grid_from_stellib_given_points(osl,
                                                                   oiso.data,
                                                                   bounds=bounds)

        # Construct the distances array. Turn single value into
        # 1-element list if single distance is given.
        _distance = np.atleast_1d(distance)
        if len(_distance) == 3:
            mindist, maxdist, stepdist = _distance
            distances = np.arange(mindist, maxdist + stepdist, stepdist)
        elif len(_distance) == 1:
            distances = np.array(_distance)
        else:
            raise ValueError("distance needs to be (min, max, step) or single number")

        # calculate the distances in pc
        if distance_unit == units.mag:
            distances = np.power(10, distances / 5. + 1) * units.pc
        else:
            distances = (distances * distance_unit).to(units.pc)

        print('applying {} distances'.format(len(distances)))

        if verbose:
            print('Adding spectral properties:', add_spectral_properties_kwargs
                  is not None)
        if add_spectral_properties_kwargs is not None:
            nameformat = add_spectral_properties_kwargs.\
                pop('nameformat', '{0:s}') + '_nd'

        # Apply the distances to the stars. Seds already at 10 pc, need
        # multiplication by the square of the ratio to this distance.
        def apply_distance_and_spectral_props(g0, distances):
            g = creategrid.apply_distance_grid(g0, distances)
            if add_spectral_properties_kwargs is not None:
                g = creategrid.add_spectral_properties(g,
                                                       nameformat=nameformat,
                                                       filterLib=filterLib,
                                                       **add_spectral_properties_kwargs)
            return g

        def grid_processing_method(g0):
            if num_dist_subgrids is None or num_dist_subgrids <= 1:
                g = apply_distance_and_spectral_props(g0, distances)
                g.writeHDF(spec_fname, append=True)
                return spec_fname, g, subgrid_names

            # Else, save the grid for separate sets of distances, with
            # or without multiprocessing
            subs = []
            q = len(distances) // num_dist_subgrids # quotient
            r = len(distances) % num_dist_subgrids # remainder
            for i in range(num_dist_subgrids):
                # spread remainder acrogss the first r subs (they get q+1 distances each)
                if i < r:
                    start = i * (q+1)
                    stop = start + q + 1
                # the rest of the subs get q distances each
                else:
                    start = r * (q+1) + (i - r) * q
                    stop = start + q
                subs.append(distances[start:stop])

            if nprocs is None:
                # Do a regular for loop
                for i, dists in enumerate(subs):
                    g = apply_distance_and_spectral_props(g0, dists)
                    g.writeHDF(subgrid_names[i], append=True)
                    g.writeHDF(spec_fname, append=True)
            else:
                # Do a parallel for loop
                print("Parallel grid construction not implemented yet")

            return subgrid_names

        # Perform the processing defined above. When the grid is
        # generated in chunks, the result of the chunks will be appended
        # to all the subgrids, and the total grid. Watch out when trying
        # this in parallel though.
        if hasattr(g, 'writeHDF'):
            subgrid_names = grid_processing_method(g)
        else:
            for gk in g:
                subgrid_names = grid_processing_method(gk)

    g = grid.FileSpectralGrid(spec_fname, backend='memory')

    return (spec_fname, g, subgrid_names)


def add_stellar_priors(project, specgrid, verbose=True,
                       priors_fname=None,
                       **kwargs):
    """
    make_priors -- compute the weights for the stellar priors

    Parameters
    ----------
    project: str
        project name

    specgrid: grid.SpectralGrid object
        spectral grid to transform

    priors_fname: str
        full filename to which to save the spectral grid with priors

    returns
    -------
    fname: str
       name of saved file

    g: grid.SpectralGrid object
        spectral grid to transform
    """
    if priors_fname is None:
        priors_fname = '%s/%s_spec_w_priors.grid.hd5' % (project, project)
    if not os.path.isfile(priors_fname):

        if verbose:
            print('Make Prior Weights')

        compute_age_mass_metallicity_weights(specgrid.grid, **kwargs)

        #write to disk
        if hasattr(specgrid, 'writeHDF'):
            specgrid.writeHDF(priors_fname)
        else:
            for gk in specgrid:
                gk.writeHDF(priors_fname, append=True)

    g = grid.FileSpectralGrid(priors_fname, backend='memory')

    return (priors_fname, g)


def make_extinguished_sed_grid(project,
                               specgrid,
                               filters,
                               av=[0., 5, 0.1],
                               rv=[0., 5, 0.2],
                               fA=None,
                               av_prior_model={'name': 'flat'},
                               rv_prior_model={'name': 'flat'},
                               fA_prior_model={'name': 'flat'},
                               extLaw=None,
                               add_spectral_properties_kwargs=None,
                               absflux_cov=False,
                               verbose=True,
                               seds_fname=None,
                               filterLib=None,
                               **kwargs):
    """
    Create SED model grid integrated with filters and dust extinguished

    Parameters
    ----------
    project: str
        project name

    specgrid: grid.SpectralGrid object
        spectral grid to transform

    filters: sequence
        ordered sequence of filters to use to extract the photometry
        filter names are the full names in core.filters

    av: sequence
        sequence of Av values to sample

    av_prior_model: list
        list including prior model name and parameters

    rv: sequence
        sequence of Rv values to sample

    rv_prior_model: list
        list including prior model name and parameters

    fA: sequence (optional)
        sequence of fA values to sample (depending on extLaw definition)

    fA_prior_model: list
        list including prior model name and parameters

    extLaw: extinction.ExtLaw
        extinction law to use during the process

    add_spectral_properties_kwargs: dict
        keyword arguments to call :func:`add_spectral_properties`
        to add model properties from the spectra into the grid property table

    asbflux_cov: boolean
        set to calculate the absflux covariance matrices for each model
        (can be very slow!!!  But it is the right thing to do)

    seds_fname: str
        full filename to save the sed grid into

    filterLib:  str
        full filename to the filter library hd5 file

    returns
    -------
    fname: str
       name of saved file

    g: grid.SpectralGrid object
        spectral grid to transform
    """
    if seds_fname is None:
        seds_fname = '%s/%s_seds.grid.hd5' % (project, project)
    if not os.path.isfile(seds_fname):

        extLaw = extLaw or extinction.Cardelli()

        avs = np.arange(av[0], av[1] + 0.5 * av[2], av[2])
        rvs = np.arange(rv[0], rv[1] + 0.5 * rv[2], rv[2])

        if verbose:
            print('Make SEDS')

        if fA is not None:
            fAs = np.arange(fA[0], fA[1] + 0.5 * fA[2], fA[2])
            g = creategrid.make_extinguished_grid(specgrid,
                                                  filters,
                                                  extLaw,
                                                  avs,
                                                  rvs,
                                                  fAs,
                                                  av_prior_model=av_prior_model,
                                                  rv_prior_model=rv_prior_model,
                                                  fA_prior_model=fA_prior_model,
                                                  add_spectral_properties_kwargs=add_spectral_properties_kwargs,
                                                  absflux_cov=absflux_cov,
                                                  filterLib=filterLib)
        else:
            g = creategrid.make_extinguished_grid(specgrid, filters, extLaw,
                                                  avs,
                                                  rvs,
                                                  av_prior_model=av_prior_model,
                                                  rv_prior_model=rv_prior_model,
                                                  add_spectral_properties_kwargs=add_spectral_properties_kwargs,
                                                  absflux_cov=absflux_cov)

        #write to disk
        if hasattr(g, 'writeHDF'):
            g.writeHDF(seds_fname)
        else:
            for gk in g:
                gk.writeHDF(seds_fname, append=True)

    g = grid.FileSEDGrid(seds_fname, backend='hdf')

    return (seds_fname, g)


def split_grid(grid_fname, num_subgrids):
    """
    Splits a spectral or sed grid (they are the same class actually)
    according to grid point index (so basically, arbitrarily).

    Parameters
    ----------
    grid_fname: string
        file name of the existing grid to be split up

    num_subgrids: integer
        the number of parts the grid should be split into

    Returns
    -------
    list of string
        the names of the newly created subgrid files

    """

    # With h5py we can choose which data we want to load to memory by
    # providing a slice
    h5grid = h5py.File(grid_fname)
    lamb = h5grid['lamb']
    seds = h5grid['seds']
    gr = h5grid['grid']

    fnames = []

    num_seds = seds.shape[0]
    q = num_seds // num_subgrids
    r = num_seds % num_subgrids
    for i in range(num_subgrids):

        subgrid_fname = grid_fname.replace('.hd5', 'sub{}.hd5'.format(i))
        fnames.append(subgrid_fname)
        if os.path.isfile(subgrid_fname):
            print('{} already exists. Skipping.'.format(subgrid_fname))
            continue
        else:
            print('constructing subgrid ' + str(i))

        # First, do strides of q+1
        if i < r:
            start = i * (q + 1)
            stop = start + q + 1
        # After the remainder has been taken care of, do strides of q
        else:
            start = r * (q + 1) + (i - r) * q
            stop = start + q

        # Load a slice as a SpectralGrid object
        slc = slice(start, stop)
        g = grid.SpectralGrid(lamb, seds=seds[slc], grid=Table(gr[slc]),
                              backend='memory')

        # Save it to a new file
        g.writeHDF(subgrid_fname, append=False)

    return fnames


def merge_grids(seds_fname, sub_names):
    """
    Merges a set of grids into one big grid. The grids need to have the
    same columns

    Parameters
    ----------
    seds_fname: string
        path for the output file

    sub_names: list of strings
        paths for the input grids
    """

    if not os.path.isfile(seds_fname):
        for n in sub_names:
            print('Appending {} to {}'.format(n, seds_fname))
            g = grid.FileSEDGrid(n)
            g.writeHDF(seds_fname, append=True)
    else:
        print('{} already exists'.format(seds_fname))
