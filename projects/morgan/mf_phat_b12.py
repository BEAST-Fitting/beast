"""
Everything I need to generate a grid, make fake data and fit

"""
# BEAST imports
from beast.core import stellib
from beast.core import extinction
from beast.core.observations import Observations
from beast.core.vega import Vega, from_Vegamag_to_Flux_SN_errors
from beast.external.ezpipe import Pipeline
from beast.external.ezpipe.helpers import task_decorator
from beast.tools.helpers import chunks
from beast.tools.pbar import Pbar
from beast.external.eztables import Table

from models import t_isochrones, t_spectra, t_seds
from fit import t_fit, t_summary_table

import os
import numpy as np


#---------------------------------------------------------
# User inputs                                   [sec:conf]
#---------------------------------------------------------
# Parameters that are required to make models
# and to fit the data
#---------------------------------------------------------

project = 'mf_phat_b12'
obsfile = 'mf_phat_b12/mf_phat_b12_input.fits'

#all filters
filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
           'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']

distanceModulus = 24.47

logt = [6.0, 10.13, 0.15]
z = [0.03, 0.019, 0.008, 0.004]


osl = stellib.Tlusty() + stellib.Kurucz()

#Extinction's law definition
extLaw = extinction.RvFbumpLaw()
avs = [0., 10.01,0.15]
rvs = [2.0, 6.1,0.5]
fbumps = [0.,1.,0.25]


#---------------------------------------------------------
# Data interface                                [sec:data]
#---------------------------------------------------------
# mapping between input cat and std_names
#---------------------------------------------------------
data_mapping = {'f275w_vega': 'HST_WFC3_F275W',
                'f275w_err': 'HST_WFC3_F275Werr',
                'f336w_vega': 'HST_WFC3_F336W',
                'f336w_err': 'HST_WFC3_F336Werr',
                'f475w_vega': 'HST_ACS_WFC_F475W',
                'f475w_err': 'HST_ACS_WFC_F475Werr',
                'f814w_vega': 'HST_ACS_WFC_F814W',
                'f814w_err': 'HST_ACS_WFC_F814Werr',
                'f110w_vega': 'HST_WFC3_F110W',
                'f110w_err': 'HST_WFC3_F110Werr',
                'f160w_vega': 'HST_WFC3_F160W',
                'f160w_err': 'HST_WFC3_F160Werr' }


#Data are in Vega magnitudes
#  Need to use Vega
with Vega() as v:
    vega_f, vega_mag, lamb = v.getMag(filters)


# derive the global class and update what's needed
class Data(Observations):
    """PHAT new photometry in M31"""
    def __init__(self, inputFile, distanceModulus=distanceModulus):
        desc = 'PHAT star: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        self.setBadValue(50.0)  # some bad values smaller than expected
        self.minError = 0.001
        self.floorError = 0.01  # constant error term to match the IDL fitter

    @from_Vegamag_to_Flux_SN_errors(lamb, vega_mag)
    def getObs(self, num):
        """ Using the decorator @from_Vegamag_to_Flux
            Hence, results are in flux (not in flux/flux_vega)
            Returns the fluxes, errors and mask of an observation.
        """
        mags = self.getMags(num, self.filters)
        errs = self.getErrors(num, self.filters)
        if self.badvalue is not None:
            mask = (mags >= self.badvalue)
        else:
            mask = np.zeros(len(mags), dtype=bool)

        #faking non symmetric errors
        return mags, errs, errs, mask

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


def get_obscat(obsfile=obsfile, distanceModulus=24.3, *args, **kwargs):
    obs = Data(obsfile, distanceModulus)
    obs.setFilters(filters)
    for k, v in data_mapping.items():
        obs.data.set_alias(v, k)
    return obs


@task_decorator()
def t_get_obscat(project, obsfile=obsfile, distanceModulus=24.3, *args, **kwargs):
    obs = get_obscat(obsfile, distanceModulus, *args, **kwargs)
    return project, obs


@task_decorator()
def t_project_dir(project, *args, **kwargs):
    outdir = project
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise Exception('Output directory "{0}" already exists but is not a directory'.format(outdir))
    else:
        os.mkdir(outdir)
    return '{0:s}/{0:s}'.format(outdir)


#---------------------------------------------------------
# Model Pipeline                                [sec:pipe]
#---------------------------------------------------------
# Create a model grid:
#     1. download isochrone(**pars)
#     2. make spectra(osl)
#     3. make seds(filters, **av_pars)
#
# Do the actual fit
#     4. load catalog(obsfile, distanceModulus)
#     5. Fit the stars
#     6. Extract statistics
#---------------------------------------------------------
def prepare_individual_inputs(obsfile, chunksize):
    """ Prepare N chuncks of data input to be run in parallel """
    if chunksize <= 0:
        return [obsfile]

    obs = Table(obsfile)
    # name will be <initial name>.<partk>.<initial_extension>
    outname = obsfile.split('.')
    outname = ('.'.join(outname[:-1]), outname[-1])

    obsfiles = []

    fpart = 0
    for chunk_slice in Pbar().iterover(chunks(range(obs.nrows), chunksize)):
        l_obs = obs[chunk_slice]
        l_file = '{0:s}.part{1:d}.{2:s}'.format(outname[0], fpart, outname[1])
        Table(l_obs).write(l_file)
        obsfiles.append(l_file)
        fpart += 1

    return obsfiles


def make_models():
    # calling sequences
    iso_kwargs = dict(logtmin=logt[0], logtmax=logt[1], dlogt=logt[2], z=z)
    spec_kwargs = dict(osl=osl)
    seds_kwargs = dict(extLaw=extLaw, av=avs, rv=rvs, fbump=fbumps)

    # make models if not there yet
    tasks_models = ( t_project_dir, t_isochrones(**iso_kwargs),  t_spectra(**spec_kwargs), t_seds(filters, **seds_kwargs) )
    models = Pipeline('make_models', tasks_models)
    job, (p, g) = models(project)
    return job, (p, g)


def run_fit(project, g, obsfile=obsfile):
    fit_kwargs = dict( threshold=-10 )
    stat_kwargs = dict( keys=None, method=['best'] )
    obscat_kwargs = dict(obsfile=obsfile, distanceModulus=distanceModulus)
    # do the real job
    tasks_fit = ( t_project_dir, t_get_obscat(**obscat_kwargs),  t_fit(g, **fit_kwargs), t_summary_table(g, **stat_kwargs) )
    fit_data = Pipeline('fit', tasks_fit)
    job, (p, stat, obs, sedgrid) = fit_data(project)
    return job, (p, stat, obs, sedgrid)


def run_chunk_fit(project, g, chunk, obsfile=obsfile):
    import glob

    obs_base = obsfile.split('.')
    obs_ext = obs_base[-1]
    obs_base = obs_base[:-1]

    lst = glob.glob('.'.join(obs_base) + '.part*.' + obs_ext)
    if len(lst) == 0:
        raise ValueError('cannot find any chunk. Did you run prepare_individual_inputs?')

    if chunk >= len(lst):
        print('Chunk not valid')

    l_obs = lst[chunk]
    print('running chunk {0:s}'.format(l_obs))

    #forcing output names
    outname = project[:]
    l_file = '{0:s}.part{1:d}'.format(outname, chunk)

    fit_kwargs = dict( threshold=-10, outname=l_file )
    stat_kwargs = dict( keys=None, method=['best'], outname=l_file)
    obscat_kwargs = dict(obsfile=l_obs, distanceModulus=distanceModulus)
    # do the real job
    tasks_fit = ( t_project_dir, t_get_obscat(**obscat_kwargs),  t_fit(g, **fit_kwargs), t_summary_table(g, **stat_kwargs) )
    fit_data = Pipeline('fit', tasks_fit)
    job, (p, stat, obs, sedgrid) = fit_data(project)
    return job, (p, stat, obs, sedgrid)


if __name__ == '__main__':
    job, (p, g) = make_models()
    job, (p, stat, obs, sedgrid) = run_fit(p, g, obsfile, distanceModulus)
