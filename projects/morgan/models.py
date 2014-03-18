"""
Grid Pipeline
=============

Everything I need to generate a grid in a few easy lines (see unittest function)

I use `ezpipe`, a pipeline package I wrote in order to clean the syntax and
allow more flexibilities. In particular it will simplifies the management of
intermediate results or broken jobs.

make models
-----------
the pipeline is sequence of tasks
tasks = ( t_isochrones(**iso_kwargs),  t_spectra(**spec_kwargs), t_seds(filters, **seds_kwargs) )

models = Pipeline('make_models', tasks)
models(project)

The pipeline is equivalent to:
grid = project | t_isochrones(**iso_kwargs) | t_spectra(**spec_kwargs) | t_seds(filters, **seds_kwargs)
"""

# system imports
import sys
import numpy as np

# BEAST imports
from beast.core import grid
from beast.core import creategrid
from beast.core import stellib
from beast.core import extinction
from beast.core import ezpadova
from beast.core.odict import odict
from beast.core.isochrone import ezIsoch
from beast.external.ezpipe.helpers import RequiredFile, task_decorator
from beast.external.ezpipe import Pipeline
from beast.external.eztables import Table

__all__ = [ 't_isochrones',  't_spectra', 't_seds' ]


def make_iso_table(outname, logtmin=6.0, logtmax=10.13, dlogt=0.05, z=0.019, **kwargs):
    """ Generate a proper table directly from the PADOVA website

    keywords
    --------

    outname: str
        file into which save the table of isochrones (any format eztables can handle)

    logtmin: float
        log-age min

    logtmax: float
        log-age max

    dlogt: float
        log-age step to request

    z: float
        unique metalicity value

    returns
    -------
    outname: str
        file into which save the table of isochrones (any format eztables can handle)
    """
    # grab a set of solar isochrones
    iso_table = ezpadova.get_t_isochrones(max(6.0, logtmin), min(10.13, logtmax), dlogt, z)
    iso_table.header['NAME'] = '{0} Isochrones'.format('_'.join(outname.split('_')[:-1]))

    #clean painful naming
    d = odict()
    drop = ['C/O', 'M_hec', 'int_IMF', 'period', 'pmode'] + "U UX B BX V R I J H K L L' M".split()
    for k in iso_table.keys():
        if not k in drop:
            if k == 'log(age/yr)':
                d['logA'] = iso_table[k]
            elif k == 'logL/Lo':
                d['logL'] = iso_table[k]
            elif k == 'logTe':
                d['logT'] = iso_table[k]
            elif k == 'logG':
                d['logg'] = iso_table[k]
            else:
                d[k] = iso_table[k]
    #add metal
    d['Z'] = np.ones(iso_table.nrows) * z

    #make a Table
    t = Table(d)
    for k, v, in iso_table.header.iteritems():
        t.header[k] = v

    # polish the header
    t.setUnit('logA', 'yr')
    t.setComment('logA', 'Age')
    t.setUnit('logT', 'K')
    t.setComment('logT', 'Effective temperature')
    t.setUnit('logL', 'Lsun')
    t.setComment('logL', 'Luminosity')
    t.setUnit('M_ini', 'Msun')
    t.setComment('M_ini', 'Initial Mass')
    t.setUnit('M_act', 'Msun')
    t.setComment('M_act', 'Current Mass, M(t)')
    t.setUnit('logMdot', 'Msun/yr')
    t.setComment('logMdot', 'Mass loss')
    t.setUnit('logg', 'cm/s**2')
    t.setComment('logg', 'Surface gravity')
    t.setComment('Z', 'Metallicity')

    #set proper aliases
    t.set_alias('logTe', 'logT')
    t.set_alias('logG', 'logg')

    t.write(outname)
    return outname


def make_spectra(outname, oiso, osl=None, **kwargs):
    """
     a spectral grid will be generated using the stellar parameters by
     interpolation of the isochrones and the generation of spectra into the
     physical units

    keywords
    --------

    outname: str
        file into which save the spectral grid

    oiso: isochrone.Isochrone object
        set of isochrones to use

    osl: stellib.Stellib object
        Spectral library to use (default stellib.Kurucz)

    returns
    -------

    outname: str
        file into which save the spectral grid
    """
    osl = osl or stellib.Kurucz()

    #filter extrapolations of the grid with given sensitivities in logg and logT
    bounds = dict(dlogT=0.1, dlogg=0.3)

    #make the spectral grid
    g = creategrid.gen_spectral_grid_from_stellib_given_points(osl, oiso.data, bounds=bounds)

    #write to disk
    g.writeHDF(outname)

    return outname


def make_seds(outname, specgrid, filters, av=[0., 5, 0.1], rv=[0., 5, 0.2], fbump=None, extLaw=None, **kwargs):
    """make_seds -- Create SED model grid integrated into filters and
    extinguished using the rest of the parameters

    keywords
    --------

    outname: str
        file into which save the final SED grid (any format grid.SpectralGrid handles)

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

    returns
    -------

    outname: str
        file into which save the SED grid
    """

    extLaw = extLaw or extinction.Cardelli()

    avs = np.arange(av[0], av[1] + 0.5 * av[2], av[2])
    rvs = np.arange(rv[0], rv[1] + 0.5 * rv[2], rv[2])
    if fbump is not None:
        fbumps = np.arange(fbump[0], fbump[1] + 0.5 * fbump[2], fbump[2])
        g = creategrid.make_extinguished_grid(specgrid, filters, extLaw, avs, rvs, fbumps)
    else:
        g = creategrid.make_extinguished_grid(specgrid, filters, extLaw, avs, rvs)
    g.writeHDF(outname)
    return outname


# =================== Pipeline Tasks ==========================

@task_decorator(logger=sys.stdout)
def t_isochrones(project, **iso_kwargs):
    """t_isochrones -- Task that generates the isochrone object

    keywords
    --------

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
    iso_source = RequiredFile(iso_fname, make_iso_table, iso_fname, **iso_kwargs)
    oiso = ezIsoch(iso_source())
    #clean header
    oiso.data.header = {'NAME': '{0} isochrones'.format(project)}
    return (project, oiso)


@task_decorator(logger=sys.stdout)
def t_spectra(project, oiso, **spec_kwargs):
    """t_spectra -- Task that generates the spectral grid

    keywords
    --------

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
    spec_source = RequiredFile(spec_fname, make_spectra, spec_fname, oiso, **spec_kwargs)
    g = grid.FileSpectralGrid(spec_source(), backend='memory')
    return project, g


@task_decorator(logger=sys.stdout)
def t_seds(project, specgrid, filters, **seds_kwargs):
    """t_seds -- Generate an sed grid from convolution with filters (incl. dust
    attenuation)

    keywords
    --------
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
    seds_source = RequiredFile(seds_fname, make_seds, seds_fname, specgrid, filters, **seds_kwargs)
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
    tasks = ( t_isochrones(**iso_kwargs),  t_spectra(**spec_kwargs), t_seds(filters, **seds_kwargs) )
    models = Pipeline('make_models', tasks)

    # run
    models(project)
