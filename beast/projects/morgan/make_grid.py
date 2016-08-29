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
from beast.core import stellib
from beast.core import extinction
from beast.external.ezpipe import Pipeline
from beast.external.ezpipe.helpers import task_decorator

from models import t_isochrones, t_spectra, t_seds

import os


"""
Global parameters
=================

Parameters that are required to make models

TODO: filters will come from observation object but currently sensitivity tests makes fake data
TODO: make a dictionary that can be easily updated
"""
project = 'mf200'

filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
           'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']

distanceModulus = 24.3

logt = [6.0, 10.13, 0.05]
#note: I do not sample mass, and use the isochrone def instead.

#TODO: isochrone only download 1 z, need to make a loop and merge outputs
z = 0.019

#osl = stellib.Tlusty() + stellib.Kurucz()
osl = stellib.Kurucz()

extLaw = extinction.Cardelli()
avs = [0., 5., 0.2]
rvs = [3.1, 3.1, 0.1]
fbumps = None


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
    outdir = project
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise Exception('Output directory "{0}" already exists but is not a directory'.format(outdir))
    else:
        os.mkdir(outdir)
    return '{0:s}/{0:s}'.format(outdir)


# actual pipeline making models
tasks_models = ( t_project_dir, t_isochrones(**iso_kwargs),  t_spectra(**spec_kwargs), t_seds(filters, **seds_kwargs) )
models = Pipeline('make_models', tasks_models)

job, (p, g) = models(project)

"""
Sensitivity Tests Pipeline
==========================

Fit a fake catalog of stars:

    1. Make fake catalog
    2. Fit the stars
    3. Extract statistics
TODO:    4. Make figures

TODO: Fit task uses the N_likelihood, should replace for SN and update Observation class
"""
from fake import t_fakedata
from fit import t_fit, t_summary_table

fake_kwargs = dict( nstars=1000, ferr=0.05, nsamp=5)
fit_kwargs = dict( threshold=-10 )
stat_kwargs = dict( keys=None, method=None )

tasks_fit = ( t_project_dir, t_fakedata(g, **fake_kwargs),  t_fit(g, **fit_kwargs), t_summary_table(g, **stat_kwargs) )
fit_fake = Pipeline('fit_fake', tasks_fit)

job, (p, stat, obs, sedgrid) = fit_fake(project)

#TODO: make figures
