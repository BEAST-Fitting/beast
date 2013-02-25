"""
Handling SED Cluster analysis of NGC4214

CONTENT:

    Imports.................................[sec:imports]
        Needed packages incl. sedfitter imports

    Globals.................................[sec:globals]
        Data file and filters

    Data interface.............................[sec:data]
        Deriving Observation class
            The class will take care of converting
            vegamags to fluxes, and bad value limits
        file column mappings

    Model definitions..........................[sec:model]

        Main ingredients.................[sec:models:main]
            define extinction, stellib, isochrones to use

        Grid sampling definition.....[sec:models:sampling]
            define model grid sampling

        Grid generation................[sec:models:create]
            create the grid if it does not exist
            according to all the assumptions

     Fitting.....................................[sec:fit]
        fitting routines

    outputs..................................[sec:outputs]
        extract minimal informations per star
"""

#---------------------------------------------------------
# Imports                                    [sec:imports]
#---------------------------------------------------------

import numpy as np
import tables
from sedfitter.observations import Observations
from sedfitter.ezisoch import ezIsoch
from sedfitter.vega import Vega, from_Vegamag_to_Flux
from sedfitter.ezpipeline import RequiredFile
from sedfitter.anased import computeLogLikelihood
from sedfitter.eztables import Table
from sedfitter.output import summaryStats
from sedfitter import grid
from sedfitter import creategrid
from sedfitter import stellib
from sedfitter import extinction
from sedfitter import progressbar

#---------------------------------------------------------
# Globals                                    [sec:globals]
#---------------------------------------------------------


obsfile = 'data/N4214_4band_detects.fits'

filters = ['HST_WFC3_F225W', 'HST_WFC3_F336W', 'HST_WFC3_F438W',
           'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']

distanceModulus = 27.41

lnp_outname = 'ngc4214_lnp.hd5'
stat_outname = 'ngc4214_stat.hd5'

#---------------------------------------------------------
# Data interface                                [sec:data]
#---------------------------------------------------------

#Data are in Vega magnitudes
#  Need to use Vega
with Vega() as v:
	vega_f, vega_mag, lamb = v.getMag(filters)


# derive the global class and update what's needed
class Data(Observations):
    """ PHAT catalog for clusters in M31 """
    def __init__(self, inputFile, distanceModulus=distanceModulus):
        desc = 'NGC4214 cluster: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        self.setBadValue(50.0)  # some bad values smaller than expected
        self.minError = 0.001
        self.floorError = 0.01  # constant error term

    @from_Vegamag_to_Flux(lamb, vega_mag)
    def getObs(self, num):
        """ Using the decorator @from_Vegamag_to_Flux
            Hence, results are in flux (not in flux/flux_vega)
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

obs = Data(obsfile, distanceModulus)
obs.setFilters(filters)
for k in filters:
    obs.data.set_alias(k, k.split('_')[-1] + '_VEGA')
    obs.data.set_alias(k + 'err', k.split('_')[-1] + '_ERR')


#---------------------------------------------------------
# Model definitions                            [sec:model]
#---------------------------------------------------------

# Main ingredients                       [sec:models:main]
#---------------------------------------------------------
isofile = 'libs/iso.proposal.fits'
spectral_grid_fname = 'libs/ngc4214.spectral.iso.fits'
sed_grid_fname = 'libs/ngc4214.sed.grid.fits'

extLaw = extinction.RvFbumpLaw()
osl = stellib.Kurucz()
oiso = ezIsoch(isofile)

# Grid sampling definition           [sec:models:sampling]
#---------------------------------------------------------

# variable to ensure that range is fully covered in using np.arange
__tiny_delta__ = 0.001

ages   = 10 ** np.arange(6., 9. + __tiny_delta__, 0.1)
masses = 10 ** np.arange(0.5, 20 + __tiny_delta__, 0.1)
Z      = [(0.02)]
avs    = np.arange(0.0, 5.0 + __tiny_delta__, 0.1)
rvs    = np.arange(1.0, 6.0 + __tiny_delta__, 0.5)
fbumps = np.asarray([1.0])

griddef = (ages, masses, Z, avs, rvs, fbumps)


# Grid generation                      [sec:models:create]
#---------------------------------------------------------
def create_sed_grid(sed_grid_fname=sed_grid_fname,
                    filter_names=filters,
                    oiso=oiso,
                    osl=osl,
                    extLaw=extLaw,
                    spectral_grid_fname=spectral_grid_fname,
                    griddef=griddef):
    """
     a spectral grid will be generated using the stellar parameters by
     interpolation of the isochrones and the generation of spectra into the
     physical units

     a photometric grid precomputing attenuation values as well
        sed_grid_fname = spectral_grid_fname.replace('spectral', 'seds')
    """
    ages, masses, Z, avs, rvs, fbumps = griddef

    #filter extrapolations of the grid with given sensitivities in logg and logT
    bounds = dict(dlogT=0.1, dlogg=0.3)

    #make the spectral grid
    creategrid.gen_spectral_grid_from_stellib(spectral_grid_fname, osl, oiso, ages=ages, masses=masses, Z=Z, bounds=bounds)

    # make the grid
    extgrid = creategrid.make_extinguished_grid(spectral_grid_fname, filter_names, extLaw, avs, rvs, fbumps)

    # save grid to file
    extgrid.write(sed_grid_fname, clobber=True)


def get_sedgrid():

    with RequiredFile(sed_grid_fname,
                    create_sed_grid,
                    filter_names=filters,
                    oiso=oiso,
                    osl=osl,
                    extLaw=extLaw,
                    spectral_grid_fname=spectral_grid_fname,
                    griddef=griddef) as __sed_grid_fname__:

            return grid.FileSEDGrid(__sed_grid_fname__)


#---------------------------------------------------------
# Fitting                                 [sec:fit]
#---------------------------------------------------------

def fit_model_seds_pytables(obs, sedgrid, threshold=-40, outname=lnp_outname):

    """
    Fit model seds with noise for sensitivity tests
    INPUTS:
        obs         Observtions     Observation object to analyze
        sedgrid     SEDgrid         stellar model SEDs (luminosities)

    KEYWORDS:
        threshold   float          toss out grid points where lnp - lnp_max < threshold
        outname     string         output file directory for results

    TODO: Clean up the disk output structures
    HDF is not a bad choice since it compress and keep all in one file!
    """
    filters = obs.getFilters()

    with tables.openFile(outname, 'w') as outfile:
        #Save wavelengths in root, remember #n_stars = root._v_nchildren -1
        outfile.createArray(outfile.root, 'grid_waves', sedgrid.lamb)
        outfile.createArray(outfile.root, 'obs_filters', filters)

        #loop over the obs and do the work
        with progressbar.PBar(len(obs), txt="Calculating lnp") as pbar:

            for tn, (sed, err, mask) in obs.enumobs():

                lnp = computeLogLikelihood(sed, err, sedgrid.seds, normed=False, mask=mask)

                #Need ragged arrays rather than uniform table
                star_group = outfile.createGroup('/', 'star_%d'  % tn, title="star %d" % tn)
                indx = np.where((lnp - max(lnp[np.isfinite(lnp)])) > -40.)
                outfile.createArray(star_group, 'input', np.array([sed, err, mask]).T)
                outfile.createArray(star_group, 'idx', np.array(indx[0], dtype=np.int32))
                outfile.createArray(star_group, 'lnp', np.array(lnp[indx[0]], dtype=np.float32))
                #commit changes
                outfile.flush()

                pbar.update(tn, force=True)  # Forcing because it can be long to show the first ETA


def do_fit():

    sedgrid = get_sedgrid()
    fit_model_seds_pytables(obs, sedgrid, threshold=-40, outname=lnp_outname)


#---------------------------------------------------------
# Outputs                                   [sec:outputs]
#---------------------------------------------------------
def get_expectation(resfname=None, sedgrid=None, filters=None, Q='logA logM Z Av Rv f_bump logT logg logL'):
    import numexpr
    # get the nD stellar/dust SED grid
    if type(sedgrid) == str:
        _sedgrid = grid.FileSEDGrid(sedgrid)
    else:
        _sedgrid = sedgrid

    if type(Q) == str:
        _Q = Q.split()
    else:
        _Q = Q

    with tables.openFile(resfname) as hdf:
        n_stars = hdf.root._v_nchildren -2
        res = np.empty((n_stars, len(_Q)), dtype=float)
        with progressbar.PBar(n_stars, txt='expectations') as pbar:
            for tn in range(n_stars):
                star_node = hdf.getNode('/star_%d'  % tn)
                lnp = star_node.lnp[:]
                idx = star_node.idx[:]
                prob = numexpr.evaluate('exp(lnp)', local_dict={'lnp': lnp})
                norm = numexpr.evaluate('sum(prob)', local_dict={'prob': prob})
                for e, qk in enumerate(_Q):
                    vals = _sedgrid.grid[qk][idx]
                    res[tn, e] = numexpr.evaluate('sum(prob/norm*vals)', local_dict={'vals': vals, 'prob': prob, 'norm': norm})
                pbar.update(tn, force=True)

        return res, _Q


def __get_Partial_Stats__(lst, resfname=None, sedgrid=None, filters=None, Q='logA logM Z Av Rv f_bump logT logg logL'):
    if lst is None:
        return
    # get the nD stellar/dust SED grid
    if type(sedgrid) == str:
        _sedgrid = grid.FileSEDGrid(sedgrid)
    else:
        _sedgrid = sedgrid

    r = []
    with tables.openFile(resfname) as hdf:
        for tn in lst:
            if tn is not None:
                star_node = hdf.getNode('/star_%d'  % tn)
                x = summaryStats(_sedgrid, star_node.idx[:], star_node.lnp[:], filters, Q)
                r.append(x.values())
    d = {}
    r = np.array(r)
    for e, k in enumerate(x.keys()):
        d[k] = r[:, e]
    return d


def grouper(n, iterable, fillvalue=None):
    from itertools import izip_longest
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def summaryTable(outname, resfname, sedgrid, obs, Q='logA logM Z Av Rv f_bump logT logg logL', nthreads=0, pool=None):

    if nthreads <= 0:
        r = __get_Partial_Stats__(range(len(obs)), resfname, sedgrid, obs, Q)
    else:
        from functools import partial
        if pool is None:
            from multiprocessing import Pool
            pool = Pool(nthreads)
        print 'pool'
        #cut into nearly equal pieces among the threads
        idx = range(len(obs))
        lst = grouper(nthreads, idx)
        print 'grouper'
        p_task = partial(__get_Partial_Stats__, resfname=resfname, sedgrid=sedgrid, filters=obs.filters, Q=Q )
        r = pool.map(p_task, lst)

    #res.set_name('SED analysis table')

    #for k in obs.data.keys()[::-1]:
    #    res.addCol(k, obs.data[k], position=0)

    #res.write(outname)
    return r


def do_summaryTable(outname=stat_outname, Q='logA logM Z Av Rv f_bump logT logg logL', **kwargs):
    return summaryTable(outname, lnp_outname, get_sedgrid(), obs, Q=Q, **kwargs)
