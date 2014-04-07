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

from models import t_isochrones, t_spectra, t_seds

from coarse_grid import *

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
iso_kwargs = dict(logtmin=logt[0], logtmax=logt[1], dlogt=logt[2], z=z)#logmhmin=logmh[0], logmhmax=logmh[1], dlogmh=logmh[2])
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

from beast.core.observations import Observations
from beast.core.vega import Vega, from_Vegamag_to_Flux
from fit import t_fit, t_summary_table
import numpy as np

with Vega() as v:
    vega_f, vega_mag, lamb = v.getMag(filters)

class PHATcatalog(Observations):
    """PHAT new photometry in M31"""
    def __init__(self, inputFile, distanceModulus=distanceModulus):
        desc = 'PHAT star: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        self.setBadValue(50.0)  # some bad values smaller than expected
        self.minError = 0.0001
        self.floorError = 0.01  # constant error term to match the IDL fitter

    @from_Vegamag_to_Flux(lamb,vega_mag)
    def getObs(self,num):
        """ Using the decorator @from_Rate_to_Flux
            Hence, results are in flux
            Returns the fluxes, errors and mask of an observation.
        """
        return Observations.getObs(self, num)

    def getObsinMag(self, num):
        """ Returns the original catalog magnitudes """
        return Observations.getObs(self, num)

    def getErrors(self, num, filters):
        """ Redifined to impose a minimal error """
        err = np.array([ self.data[tt + 'err'][num] for tt in filters])
        if self.floorError > 0.:
            err = np.sqrt(err ** 2 + self.floorError ** 2)
        if self.minError > min(err):
            err[ err < self.minError ] = self.minError
        return err


obs = PHATcatalog(obsfile, distanceModulus)
obs.setFilters(filters)
for k in filters:
    obs.data.set_alias(k, k.split('_')[-1].lower() + '_vega')
    obs.data.set_alias(k + 'err', k.split('_')[-1].lower() + '_err')

fit_kwargs = dict( threshold=-20 )
stat_kwargs = dict( keys=None, method=None )

tasks_fit = ( t_project_dir, t_fit(obs, g, **fit_kwargs), t_summary_table(g, **stat_kwargs) )
fit_phat = Pipeline('fit_phat', tasks_fit)

job, (p, stat, obs, sedgrid) = fit_phat(project)

