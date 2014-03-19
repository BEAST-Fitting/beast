"""
Based on make_grid.py

Testing how Stellib classes react to merging through CompositeStellib.
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

No need for fine gridding we are just testing the functions
"""

project = 'mf_osl'

filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
           'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']

distanceModulus = 24.3

logt = [6.0, 10.13, 0.1]
z = 0.019
# set the extinction law and parameters sampling (min, max, step)
extLaw = extinction.Cardelli()
avs = [0., 1., 1.]
rvs = [3.1, 3.1, 0.1]
fbumps = None  # cardelli = fixed bump (no variations)

# Make a composite library
osl = stellib.Tlusty() + stellib.Kurucz()
#extrapolation flexibility
bounds = dict(dlogT=0.1, dlogg=0.3)


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
spec_kwargs = dict(osl=osl, bounds=bounds)
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
tasks_models = (t_project_dir,
                t_isochrones(**iso_kwargs),
                t_spectra(**spec_kwargs),
                t_seds(filters, **seds_kwargs) )

models = Pipeline('make_models', tasks_models)

job, (p, g) = models(project)
