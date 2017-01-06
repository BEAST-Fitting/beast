"""
Grid Pipeline
=============

Everything needed to generate a grid in a few easy lines
   (see unittest function)

`ezpipe` is used, a pipeline package allowing a clean the syntax and
flexibility. In particular it simplifies the management of
intermediate results or broken jobs.

make models
-----------
the pipeline is sequence of tasks

for example:
    
tasks = ( t_isochrones(**iso_kwargs),  t_spectra(**spec_kwargs),
  t_seds(filters, **seds_kwargs) )

models = Pipeline('make_models', tasks)
models(project)

The pipeline is equivalent to:
grid = project | t_isochrones(**iso_kwargs) | t_spectra(**spec_kwargs) |
       t_seds(filters, **seds_kwargs)
"""

# system imports
from __future__ import print_function
import sys
import numpy as np

# BEAST imports
import grid
import creategrid
import prior_weights

from stars import stellib
from stars import isochrone
from stars.isochrone import ezIsoch

from dust import extinction

from ..external.ezpipe.helpers import RequiredFile, task_decorator
from ..external.ezpipe import Pipeline
from ..tools.helpers import val_in_unit

__all__ = ['make_iso_table','make_spectra', 'make_priors', 'make_seds']
#,
#           't_isochrones',  't_spectra', 't_priors', 't_seds' ]

def make_iso_table(outname, logtmin=6.0, logtmax=10.13, dlogt=0.05,
                   z=[0.019], trackVersion=2.3):
    """ Generate a proper table directly from the PADOVA website

    Parameters
    ----------

    outname: str
        file into which save the table of isochrones (any format
        eztables can handle)

    logtmin: float
        log-age min

    logtmax: float
        log-age max

    dlogt: float
        log-age step to request

    z: float or sequence
        list of metalicity values

    Returns
    -------
    outname: str
        file into which save the table of isochrones (any format eztables
        can handle)
    """
    oiso = isochrone.PadovaWeb()
    t = oiso._get_t_isochrones(max(6.0, logtmin), min(10.13, logtmax),
                               dlogt, z, trackVersion)
    t.header['NAME'] = '{0} Isochrones'.\
                       format('_'.join(outname.split('_')[:-1]))

    # Isochrone filtering, check that no AGB stars are removed
    if trackVersion < 2.7:  
        cond = '~((logL > 3.) & (M_act < 1.) & (log10(M_ini / M_act) > 0.1))'
    else:  
        cond = '~((M_ini < 12.) & (stage == 0))' # do not include Pre-MS  
    t = t.selectWhere('*', cond)

    t.write(outname)
    return outname


def make_spectra(outname, oiso, osl=None, bounds={}, distance=None,
                 add_spectral_properties_kwargs=None, **kwargs):
    """
     a spectral grid will be generated using the stellar parameters by
     interpolation of the isochrones and the generation of spectra into the
     physical units

    Parameters
    ----------

    outname: str
        file into which save the spectral grid

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
    osl = osl or stellib.Kurucz()

    #filter extrapolations of the grid with given sensitivities in logg
    #  and logT
    if 'dlogT' not in bounds:
        bounds['dlogT'] = 0.1
    if 'dlogg' not in bounds:
        bounds['dlogg'] = 0.3

    #make the spectral grid
    print('Make spectra')
    g = creategrid.gen_spectral_grid_from_stellib_given_points(osl,
                                                               oiso.data,
                                                               bounds=bounds)

    # get the distance
    if distance is not None:
        _distance = val_in_unit('distance', distance, 'pc').magnitude

    print('Adding spectral properties:', add_spectral_properties_kwargs
          is not None)
    if add_spectral_properties_kwargs is not None:
        nameformat = add_spectral_properties_kwargs.\
                     pop('nameformat', '{0:s}') + '_nd'

    #write to disk
    # and apply the distance to the particular galaxy of interest
    # seds already at 10 pc, need multiplcation by the square of the ratio
    # to this distance
    if hasattr(g, 'writeHDF'):
        if distance is not None:
            g.seds = g.seds / (0.1 * _distance) ** 2
        if add_spectral_properties_kwargs is not None:
            g = creategrid.add_spectral_properties(g, nameformat=nameformat,
                                            **add_spectral_properties_kwargs)
        g.writeHDF(outname)
    else:
        for gk in g:
            if distance is not None:
                gk.seds = gk.seds / (0.1 * _distance) ** 2
            if add_spectral_properties_kwargs is not None:
                gk = creategrid.add_spectral_properties(gk,
                                            nameformat=nameformat,
                                            **add_spectral_properties_kwargs)
            gk.writeHDF(outname, append=True)

    return outname

def make_priors(outname, specgrid, **kwargs):
    """make_priors -- compute the weights for the priors

    Parameters
    ----------

    outname: str
        file into which save the final SED grid (any format
        grid.SpectralGrid handles)

    specgrid: grid.SpectralGrid object
        spectral grid to transform
        result from the make_spectra function

    returns
    -------

    outname: str
        file into which save the SED grid
    """

    print('Make Prior Weights')

    prior_weights.compute_age_mass_metallicity_prior_weights(specgrid.grid)

    #write to disk
    if hasattr(specgrid, 'writeHDF'):
        specgrid.writeHDF(outname)
    else:
        for gk in specgrid:
            gk.writeHDF(outname, append=True)
    return outname

def make_seds(outname, specgrid, filters, av=[0., 5, 0.1], rv=[0., 5, 0.2],
              fbump=None, extLaw=None, add_spectral_properties_kwargs=None,
              absflux_cov=False,
              **kwargs):
    """make_seds -- Create SED model grid integrated into filters and
    extinguished using the rest of the parameters

    Parameters
    ----------

    outname: str
        file into which save the final SED grid (any format
        grid.SpectralGrid handles)

    specgrid: grid.SpectralGrid object
        spectral grid to transform
        result from the make_spectra function

    filters: sequence
        ordered sequence of filters to use to extract the photometry
        filter names are the full names in core.filters

    av: sequence
        sequence of Av values to sample

    rv: sequence
        sequence of Rv values to sample

    fbump: sequence (optional)
        sequence of fbump values to sample (depending on extLaw definition)

    extLaw: extinction.ExtLaw
        extinction law to use during the process

    add_spectral_properties_kwargs: dict
        keyword arguments to call :func:`add_spectral_properties`
        to add model properties from the spectra into the grid property table

    asbflux_cov: boolean
        set to calculate the absflux covariance matrices for each model
        (can be very slow!!!  But it is the right thing to do)

    returns
    -------

    outname: str
        file into which save the SED grid
    """

    extLaw = extLaw or extinction.Cardelli()

    avs = np.arange(av[0], av[1] + 0.5 * av[2], av[2])
    rvs = np.arange(rv[0], rv[1] + 0.5 * rv[2], rv[2])

    print('Make SEDS')

    if fbump is not None:
        fbumps = np.arange(fbump[0], fbump[1] + 0.5 * fbump[2], fbump[2])
        g = creategrid.make_extinguished_grid(specgrid, filters, extLaw,
                                              avs, rvs, fbumps,
                add_spectral_properties_kwargs=add_spectral_properties_kwargs,
                                              absflux_cov=absflux_cov)
    else:
        g = creategrid.make_extinguished_grid(specgrid, filters, extLaw,
                                              avs, rvs,
                add_spectral_properties_kwargs=add_spectral_properties_kwargs,
                                              absflux_cov=absflux_cov)

    #write to disk
    if hasattr(g, 'writeHDF'):
        g.writeHDF(outname)
    else:
        for gk in g:
            gk.writeHDF(outname, append=True)
    return outname

# =================== Pipeline Tasks ==========================

@task_decorator(logger=sys.stdout)
def t_isochrones(project, **iso_kwargs):
    """t_isochrones -- Task that generates the isochrone object

    Parameters
    ----------

    project: str
        token of the project this task belongs to

    **iso_kwargs: dict
        arguments to pass to make_iso_table

    returns
    -------
    project: str
       project token that needs to be kept along the task

    oiso: isochrone.Isochrone object
        Isochrone object instance
    """
    iso_fname = '{0}_iso.csv'.format(project)
    iso_source = RequiredFile(iso_fname, make_iso_table, iso_fname,
                              **iso_kwargs)
    oiso = ezIsoch(iso_source())
    #clean header
    oiso.data.header = {'NAME': '{0} isochrones'.format(project)}
    return (project, oiso)


@task_decorator(logger=sys.stdout)
def t_spectra(project, oiso, **spec_kwargs):
    """t_spectra -- Task that generates the spectral grid

    Parameters
    ----------

    project: str
        token of the project this task belongs to

    oiso: isochrone.Isochrone object
        Isochrone object instance

    **spec_kwargs:
        any arguments forwarded to make_spectra

    returns
    -------
    project: str
       project token that needs to be kept along the task

    g: grid.SpectralGrid instance
        spectral grid instance
    """
    spec_fname = '{0}_spec.grid.hd5'.format(project)
    spec_source = RequiredFile(spec_fname, make_spectra, spec_fname,
                               oiso, **spec_kwargs)
    g = grid.FileSpectralGrid(spec_source(), backend='memory')
    
    return project, g

@task_decorator(logger=sys.stdout)
def t_priors(project, specgrid, **priors_kwargs):
    """t_priors -- Task that updates the weights to include the priors

    Parameters
    ----------

    project: str
        token of the project this task belongs to

    specgrid: grid.SpectralGrid instance
        spectral grid instance

    filters: sequence
        sequence of filter standard names

    **priors_kwargs:
        any arguments forwarded to compute_priors

    returns
    -------
    project: str
       project token that needs to be kept along the task

    g: grid.SpectralGrid instance
        spectral grid instance
    """
    priors_fname = '{0}_spec_w_priors.grid.hd5'.format(project)
    priors_source = RequiredFile(priors_fname, make_priors, priors_fname,
                                 specgrid, **priors_kwargs)
    g = grid.FileSpectralGrid(priors_source(), backend='memory')
    
    return project, g

@task_decorator(logger=sys.stdout)
def t_seds(project, specgrid, filters, **seds_kwargs):
    """t_seds -- Generate an sed grid from convolution with filters (incl. dust
    attenuation)

    Parameters
    ----------
    project: str
        token of the project this task belongs to

    specgrid: grid.SpectralGrid instance
        spectral grid instance

    filters: sequence
        sequence of filter standard names

    **seds_kwargs:
        any arguments forwarded to make_seds

    returns
    -------
    project: str
       project token that needs to be kept along the task

    g: grid.SpectralGrid instance
        SED grid instance
    """

    seds_fname = '{0}_seds.grid.hd5'.format(project)
    seds_source = RequiredFile(seds_fname, make_seds, seds_fname,
                               specgrid, filters, **seds_kwargs)
    g = grid.FileSEDGrid(seds_source(), backend='hdf')
    return project, g

# =================== Pipeline Example ==========================

def unittest():
    """
    Global parameters
    =================

    """
    #project token name, every task will use it as id
    project = 'unittest'

    #define the filters to use to generate an SED grid. i
    # Names are standard in the package (libs/filters.hd5/content).
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']

    # define the sampling in log(age): min, max, step
    logt = [6.0, 10.13, 0.05]

    # set the stellar library
    osl = stellib.Kurucz()

    # set the extinction law and parameters sampling (min, max, step)
    extLaw = extinction.Cardelli()
    avs = [0., 5., 0.2]
    rvs = [3.1, 3.1, 0.1]
    fbumps = None  # cardelli = fixed bump (no variations)

    # generate task respective dictionaries
    iso_kwargs = dict(logtmin=logt[0], logtmax=logt[1], dlogt=0.05, z=0.019)
    spec_kwargs = dict(osl=osl)
    seds_kwargs = dict(extLaw=extLaw, av=avs, rv=rvs, fbump=fbumps)

    # ======================= Pipeline ==============================
    # actual pipeline making models
    tasks = ( t_isochrones(**iso_kwargs),  t_spectra(**spec_kwargs),
              t_seds(filters, **seds_kwargs) )
    models = Pipeline('make_models', tasks)

    # run
    models(project)
