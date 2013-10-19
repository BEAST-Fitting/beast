"""
Everything I need to generate a grid


TODO: make better documentation
TODO: maybe take all subfunctions into individual package
"""

# system imports
import sys
import numpy as np
#import tables

# BEAST imports
#from beast.core.observations import Observations
from beast.core.isochrone import ezIsoch
#from beast.core.vega import Vega, from_Vegamag_to_Flux
#from beast.core.anased import computeLogLikelihood
from beast.core import grid
from beast.core import creategrid
from beast.core import stellib
from beast.core import extinction

# morgan imports
import ezpadova
from ezpipe.helpers import RequiredFile, task_decorator
from ezpipe import Pipeline
from beast.external.eztables import Table
from beast.external.eztables.core.odict import odict

__all__ = [ 't_isochrones',  't_spectra', 't_seds' ]


def make_iso_table(outname, logtmin=6.0, logtmax=10.13, dlogt=0.05, z=0.019):
    """ Generate a proper table directly from the PADOVA website """
    # grab a set of solar isochrones
    iso_table = ezpadova.get_t_isochrones(max(6.0, logtmin), min(10.13, logtmax), dlogt, z)
    iso_table.header['NAME'] = '{} Isochrones'.format('_'.join(outname.split('_')[:-1]))

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


def make_spectra(outname, oiso, extLaw=None, osl=None):
    """
     a spectral grid will be generated using the stellar parameters by
     interpolation of the isochrones and the generation of spectra into the
     physical units
    """
    osl = osl or stellib.Kurucz()

    #filter extrapolations of the grid with given sensitivities in logg and logT
    bounds = dict(dlogT=0.1, dlogg=0.3)

    #make the spectral grid
    creategrid.gen_spectral_grid_from_stellib_given_points(outname, osl, oiso, oiso.data, bounds=bounds)

    return outname


def make_seds(outname, specgrid, filters, av=[0., 5, 0.1], rv=[0., 5, 0.2], fbump=None, extLaw=None):
    extLaw = extLaw or extinction.Cardelli()

    avs = np.arange(av[0], av[1] + 0.5 * av[2], av[2])
    rvs = np.arange(rv[0], rv[1] + 0.5 * rv[2], rv[2])
    if fbump is not None:
        fbumps = np.arange(fbump[0], fbump[1] + 0.5 * fbump[2], fbump[2])
        g = creategrid.make_extinguished_grid(specgrid, filters, extLaw, avs, rvs, fbumps)
    else:
        g = creategrid.make_extinguished_grid(specgrid, filters, extLaw, avs, rvs)
    g.write(outname)
    return outname


"""
Pipeline
========
I will use the pipeline package I wrote in order to clean the syntax and allow
more flexibilities. In particular it will simplifies the managment of
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


# =================== Pipeline Tasks ==========================

@task_decorator(logger=sys.stdout)
def t_isochrones(project, **iso_kwargs):
    iso_fname = '{}_iso.csv'.format(project)
    iso_source = RequiredFile(iso_fname, make_iso_table, iso_fname, **iso_kwargs)
    oiso = ezIsoch(iso_source())
    #clean header
    oiso.data.header = {'NAME': '{} isochrones'.format(project)}
    return (project, oiso)


@task_decorator(logger=sys.stdout)
def t_spectra(project, oiso, **spec_kwargs):
    spec_fname = '{}_spec.grid.fits'.format(project)
    spec_source = RequiredFile(spec_fname, make_spectra, spec_fname, oiso, **spec_kwargs)
    g = grid.FileSpectralGrid(spec_source())
    return project, g


@task_decorator(logger=sys.stdout)
def t_seds(project, specgrid, filters, **seds_kwargs):
    seds_fname = '{}_seds.grid.fits'.format(project)
    seds_source = RequiredFile(seds_fname, make_seds, seds_fname, specgrid, filters, **seds_kwargs)
    g = grid.FileSEDGrid(seds_source())
    return project, g


def unittest():
    """
    Global parameters
    =================

    """
    project = 'unittest'

    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']

    logt = [6.0, 10.13, 0.05]

    osl = stellib.Kurucz()

    extLaw = extinction.Cardelli()
    avs = [0., 5., 0.2]
    rvs = [3.1, 3.1, 0.1]
    fbumps = None

    iso_kwargs = dict(logtmin=logt[0], logtmax=logt[1], dlogt=0.05, z=0.019)
    spec_kwargs = dict(osl=osl)
    seds_kwargs = dict(extLaw=extLaw, av=avs, rv=rvs, fbump=fbumps)

    # ======================= Pipeline ==============================
    # actual pipeline making models
    tasks = ( t_isochrones(**iso_kwargs),  t_spectra(**spec_kwargs), t_seds(filters, **seds_kwargs) )
    models = Pipeline('make_models', tasks)

    models(project)
