"""
Everything I need to generate a grid, make fake data and fit

I use the pipeline package I wrote in order to clean the syntax and allow more
flexibilities. In particular it will simplifies the managment of intermediate
results or broken jobs.

make models
-----------
the pipeline is sequence of tasks
tasks = ( t_isochrones(**iso_kwargs),  t_spectra(**spec_kwargs), t_seds(filters, **seds_kwargs) )

models = Pipeline('make_models', tasks)
job, (project, grid) = models(project)

The pipeline is equivalent to using individual tasks as follow:
project, grid = project | t_isochrones(**iso_kwargs) | t_spectra(**spec_kwargs) | t_seds(filters, **seds_kwargs)

However it allows more management options, including logs, memoizations etc.

TODO: make better documentation
TODO: make all static code go into a different module
"""


# BEAST imports
#from beast.core import stellib
#from beast.core import extinction
from beast.external.ezpipe import Pipeline
from beast.external.ezpipe.helpers import task_decorator
from beast.core.observations import Observations
from beast.core.vega import Vega
from fit import t_fit, t_summary_table

from models import t_isochrones, t_spectra, t_seds
import tables

from GRID_PATH_DEF import *  # TO UPDATE


import os


"""
Model Pipeline
==============

Create a model grid:

    1. download isochrone(**pars)
    2. make spectra(osl)
    3. make seds(filters, **av_pars)

each step outputs results that are stored into <project>_<...>.<csv|fits>

Just call "models(project)" to make a model grid

TODO: make a function that takes user pars and return the pipeline instance
"""

# calling sequences
iso_kwargs = dict(logtmin=logt[0], logtmax=logt[1], dlogt=logt[2], z=z)
spec_kwargs = dict(osl=osl)
seds_kwargs = dict(extLaw=extLaw, av=avs, rv=rvs, fbump=fbumps)


@task_decorator()
def t_project_dir(project, *args, **kwargs):
    #dir = '/astro/dust_kg/harab/beast_last/projects/morgan/'
    outdir = project
    print outdir
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise Exception('Output directory "{0}" already exists but is not a directory'.format(outdir))
    else:
        os.mkdir(outdir)
    return '{0:s}/{1:s}'.format(outdir,project)


# actual pipeline making models
tasks_models = ( t_project_dir, t_isochrones(**iso_kwargs),  t_spectra(**spec_kwargs), t_seds(filters, **seds_kwargs) )
models = Pipeline('make_models', tasks_models)

job, (p, g) = models(project)

"""
Data interface
==============
Run the grid on specific data
TO DO: make figures
"""

with Vega() as v:
    vega_f, vega_mag, lamb = v.getMag(filters)


class PHATcatalog(Observations):
    """PHAT 6 filter photometry
    This class implements a direct access to the PHAT measured fluxes.

    ..note::
        it does not implement uncertainties as in this model, the noise is
        given through artificial star tests
    """
    def __init__(self, inputFile, distanceModulus=distanceModulus, filters=filters):
        """ Construct the interface """
        desc = 'PHAT star: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        #some bad values smaller than expected
        # in RATE = flux units
        self.setBadValue(6e-11)

        #hard code mapping directly with the interface to PHAT
        for k in filters:
            self.data.set_alias(k, k.split('_')[-1].lower() + '_rate')

    def getFlux(self, num, units=False):
        """returns the absolute flux of an observation from the number of
        counts

        Parameters
        ----------
        num: int
            index of the star in the catalog to get measurement from

        units: bool
            if set returns the fluxes with a unit capsule

        Returns
        -------
        flux: ndarray[dtype=float, ndim=1]
            Measured integrated flux values throughout the filters in erg/s/cm^2
        """

        flux = Observations.getFlux(self, num) * self.vega_flux
        if units is True:
            return flux * unit['erg/s/cm^2']
        else:
            return flux

    def getFluxerr(self, num, units=False):
        """returns the error on the absolute flux of an observation from the
        number of counts (not used in the analysis)

        Parameters
        ----------
        num: int
            index of the star in the catalog to get measurement from

        units: bool
            if set returns the fluxes with a unit capsule

        Returns
        -------
        fluxerr: ndarray[dtype=float, ndim=1]
            Measured integrated flux uncertainties in erg/s/cm^2
        """

        fluxerr = Observations.getFluxerr(self, num) * self.vega_flux
        if units is True:
            return fluxerr * unit['erg/s/cm^2']
        else:
            return fluxerr

    def setFilters(self, filters):
        """ set the filters and update the vega reference for the conversions

        Parameters
        ----------
        filters: sequence
            list of filters using the internally normalized namings
        """
        self.filters = filters

        #Data "rates" are normalized to Vega already, fits are not using vega

        # for optimization purpose: pre-compute
        #   getting vega mags, require to open and read the content of one file.
        #   since getObs, calls getFlux, for each star you need to do this expensive
        #   op.
        with Vega() as v:
            _, vega_flux, _ = v.getFlux(filters)

        self.vega_flux = vega_flux


obs = PHATcatalog(obsfile, distanceModulus)
obs.setFilters(filters)
for k in filters:
    obs.data.set_alias(k, k.split('_')[-1].lower() + '_rate')

ast = tables.openFile(astfile)

fit_kwargs = dict( threshold=-10 )
stat_kwargs = dict( keys=None, method=None )

tasks_fit = ( t_project_dir, t_fit(obs, g, ast, **fit_kwargs), t_summary_table(g, **stat_kwargs) )
fit_phat = Pipeline('fit_phat', tasks_fit)

job, (p, stat, obs, sedgrid) = fit_phat(project)

#TODO: make figures
