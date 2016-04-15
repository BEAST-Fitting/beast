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
from beast.tools.helpers import chunks
from beast.tools.helpers import val_in_unit
from beast.tools.pbar import Pbar
from beast.external.ezunits import unit

import os

import datamodel
import noisemodel
from models import t_isochrones, t_spectra, t_seds, t_priors
from noisemodel import t_gen_noise_model


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


def prepare_individual_inputs(obsfile, chunksize=14000):
    """ Prepare N chuncks of observation input to be run in parallel

    Parameters
    ----------
    obsfile: fname
        input file containing observations

    chunksize: int
        number of sources per chunk of data
        (default number is based on maximized HDF5 node numbers/speed ratio)

    Returns
    -------
    obsfiles: sequence
        list of created files
        Files are generated with the given number of sources per individual catalog.
        Namings respects this convention: `<initial name>.<partk>.<initial_extension>`

    .. note::
        this function uses `beast.external.eztables` to read the catalog with
        the largest flexibility on the input format

    .. todo::
        loading the catalog could be optimized for memory usage.
        In practice this could also just be replaced by database queries.
    """
    if chunksize <= 0:
        return [obsfile]

    obs = Table(obsfile)
    # name will be <initial name>.<partk>.<initial_extension>
    outname = obsfile.split('.')
    outname = ('.'.join(outname[:-1]), outname[-1])

    obsfiles = []

    fpart = 0
    for chunk_slice in Pbar(desc='Preparing input catalogs').iterover(chunks(range(obs.nrows), chunksize)):
        l_obs = obs[chunk_slice]
        l_file = '{0:s}.part{1:d}.{2:s}'.format(outname[0], fpart, outname[1])
        Table(l_obs).write(l_file)
        obsfiles.append(l_file)
        fpart += 1

    return obsfiles


def merge_individual_outputs(obsfile=datamodel.obsfile, project=datamodel.project):
    """
    Parameters
    ----------
    obsfile: fname
        input file containing observations

    project: str
        project name that will identify the individual parts

    Returns
    -------
    obs: Table
        final catalog Table
    """
    from glob import glob
    lst = glob('./{project:s}/{project:s}_stats_part*.fits'.format(project=project))
    N = len(lst)
    fname = './{project:s}/{project:s}_stats_part{chunk:d}.fits'
    t = Table(fname.format(project=project, chunk=0))
    for c in range(1, N):
        _fname = fname.format(project=project, chunk=c)
        print('Processing table: {0:s}'.format(_fname))
        t.stack(Table(_fname).data)
    t.write('./{project:s}/{project:s}_stats.fits'.format(project=project))
    print('Merging with input catalog')
    obs = Table(datamodel.obsfile)
    for k in t.keys():
        obs.addCol(k, t[k])
    obs.write('./{project:s}/{project:s}_merged.fits'.format(project=project))

    return obs


def make_models(*args, **kwargs):
    """ generates models from scratch

    1. creates the project directory,
    2. download isochrones,
    3. compute dust-free spectra from the isochrones
    4. apply dust and generate photometry to obtain the set of seds
    5. generate the noise model

    Equivalent to using individual tasks as follow:
    project, noisefile, grid = project | t_isochrones(**iso_kwargs)
                                       | t_spectra(**spec_kwargs)
                                       | t_priors()
                                       | t_seds(filters, **seds_kwargs)
                                       | t_gen_noise_model(astfile, **noise_kwargs)

    returns
    -------
    job: int
        job id

    (p, n, g): project, noisefile and ModelGrid
        project identification
        noise model file
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

    spec_kwargs = dict(osl=datamodel.osl, distance=distance,
                       extLaw=datamodel.FGextLaw, av=datamodel.FGextLaw)

    seds_kwargs = dict(extLaw=datamodel.extLaw,
                       av=datamodel.avs,
                       rv=datamodel.rvs,
                       fbump=datamodel.fbumps,
                       absflux_cov=datamodel.moddep_absflux_cov)

    if hasattr(datamodel, 'add_spectral_properties_kwargs'):
        seds_kwargs['add_spectral_properties_kwargs'] = datamodel.add_spectral_properties_kwargs
        spec_kwargs['add_spectral_properties_kwargs'] = datamodel.add_spectral_properties_kwargs

    #noise_kwargs = dict(covariance=datamodel.absflux_a_matrix)

    # make models if not there yet
    tasks_models = (t_project_dir,
                    t_isochrones(**iso_kwargs),
                    t_spectra(**spec_kwargs),
                    t_priors(),
                    t_seds(datamodel.filters, **seds_kwargs))

    models = Pipeline('make_models', tasks_models)
    job, (p, g) = models(datamodel.project)

    return job, (p, g)

def compute_noise_and_trim_grid(*args, **kwargs):
    """ compute the noise model for a specific catalog and
        trim bright/faint models to speed up the calcuations

    1. compute the noise model given the ASTs and absflux A matrix
    2. trim the model sed grid of models that are so bright
       or faint that will always get a zero fit probability (speed optimization)

    Equivalent to using individual tasks as follow:
    project, noisefile, grid = project | t_gen_noise_model(datamodel.astfile, **noise_kwargs)
                                       | t_trim_model_sed(????)

    returns
    -------
    job: int
        job id

    (p, n, g): project, noisemodel, ModelGrid
        project identification
        noise model file
        Modelgrid instance constaining the collection of SEDs
    """
    # calling sequences

    noise_kwargs = dict(covariance=datamodel.absflux_a_matrix)

    # make models if not there yet
    tasks_noise_and_trim = (t_gen_noise_model(datamodel.astfile, **noise_kwargs))

    noise_and_trim = Pipeline('compute_noise_and_trim_grid', tasks_noise_and_trim)
    job, (p, n, g) = noise_and_trim(datamodel.project)

    return job, (p, n, g)
