"""
Model Pipeline
==============

Create a model grid:

    1. download isochrone(**pars)
    2. make spectra(osl)
    3. make seds(filters, **av_pars)

each step outputs results that are stored into <project>_<...>.<csv|fits>

TODO: make a function that takes user pars and return the pipeline instance
"""
from __future__ import print_function
from beast.external.ezpipe import Pipeline
from beast.external.ezpipe.helpers import task_decorator
from beast.external.eztables import Table
from beast.tools.helpers import val_in_unit
from beast.tools.pbar import Pbar
from beast.external.ezunits import unit

import os

import datamodel_small as datamodel
import noisemodel
from beast.fitting.models import t_isochrones, t_spectra, t_seds, t_priors

@task_decorator()
def t_get_obscat(project, obsfile=datamodel.obsfile,
                 distanceModulus=datamodel.distanceModulus,
                 filters=datamodel.filters,
                 *args, **kwargs):
    """ task that generates a data catalog object with the correct arguments

    Parameters
    ----------
    obsfile: str, optional (default datamodel.obsfile)
        observation file

    distanceModulus: float, optional (default datamodel.distanceModulus)
        distance modulus to correct the data from (in magitude)

    filters: sequence(str), optional, datamodel.filters
        seaquence of filters of the data

    returns
    -------
    project: str
        project id

    obs: PHATFluxCatalog
        observation catalog
    """
    obs = datamodel.get_obscat(obsfile, distanceModulus, filters,
                               *args, **kwargs)
    return project, obs


@task_decorator()
def t_project_dir(project, *args, **kwargs):
    """ Task that creates the project directory if necessary

    Parameters
    ----------
    project: str
        project name

    Returns
    -------
    dirname: str
        <project>/<project>

    Raises
    ------
    Exception
        if already exists a file that is not a directory
    """
    outdir = project
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise Exception('Output directory "{0}" already exists but is not a directory'.format(outdir))
    else:
        os.mkdir(outdir)
    return '{0:s}/{1:s}'.format(outdir, project)


def make_models(*args, **kwargs):
    """ generates models from scratch

    1. creates the project directory,
    2. download isochrones,
    3. compute dust-free spectra from the isochrones
    4. compute the age-mass-metallicity prior weights
    5. apply dust and generate photometry to obtain the set of seds

    Equivalent to using individual tasks as follow:
    project, noisefile, grid = project | t_isochrones(**iso_kwargs)
                                       | t_spectra(**spec_kwargs)
                                       | t_priors()
                                       | t_seds(filters, **seds_kwargs)

    returns
    -------
    job: int
        job id

    (p, g): project and ModelGrid
        project identification
        Modelgrid instance constaining the collection of SEDs
    """
    # calling sequences
    iso_kwargs = dict(logtmin=datamodel.logt[0],
                      logtmax=datamodel.logt[1],
                      dlogt=datamodel.logt[2],
                      z=datamodel.z,
                      trackVersion=datamodel.trackVersion)

    dmod = val_in_unit('distance Modulus', datamodel.distanceModulus, 'mag').magnitude
    distance = 10 ** ( (dmod / 5.) + 1 ) * unit['pc']

    spec_kwargs = dict(osl=datamodel.osl, distance=distance)

    seds_kwargs = dict(extLaw=datamodel.extLaw,
                       av=datamodel.avs,
                       rv=datamodel.rvs,
                       fbump=datamodel.fbumps
                       )

    if hasattr(datamodel, 'add_spectral_properties_kwargs'):
        seds_kwargs['add_spectral_properties_kwargs'] = datamodel.add_spectral_properties_kwargs
        spec_kwargs['add_spectral_properties_kwargs'] = datamodel.add_spectral_properties_kwargs

    # make models if not there yet
    tasks_models = (t_project_dir,
                    t_isochrones(**iso_kwargs),
                    t_spectra(**spec_kwargs),
                    t_priors(),
                    t_seds(datamodel.filters, **seds_kwargs))

    models = Pipeline('make_models', tasks_models)
    job, (p, g) = models(datamodel.project)

    return job, (p, g)


