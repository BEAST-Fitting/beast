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

        Grid sampling definition.........[sec:models:grid]
            define model grid sampling

        Grid generation................[sec:models:create]
            create the grid if it does not exist
            according to all the assumptions

     Fitting.....................................[sec:fit]
        fitting routines

    outputs..................................[sec:outputs]
        extract minimal informations per star

    figures..................................[sec:figures]
        currently diagnostic plots

    pipeline................................[sec:pipeline]
        fast aliases checking the previous steps
        all in one call (assuming all params are ok)
        >>> run()
"""

#---------------------------------------------------------
# Imports                                    [sec:imports]
#---------------------------------------------------------
try:
    import androfig
    __FIG__ = True
except ImportError:
    __FIG__ = False

import numpy as np
import tables
from sedfitter.core.observations import Observations
from sedfitter.core.isochrone import ezIsoch, padova2010
from sedfitter.core.vega import Vega, from_Vegamag_to_Flux
from sedfitter.tools.ezpipeline import RequiredFile
from sedfitter.core.anased import computeLogLikelihood
from sedfitter.external.eztables import AstroTable
from sedfitter.core import grid
from sedfitter.core import creategrid
from sedfitter.core import stellib
from sedfitter.core import extinction
from sedfitter.tools import progressbar
from sedfitter.tools import binningAxis

#---------------------------------------------------------
# Globals                                    [sec:globals]
#---------------------------------------------------------


#obsfile = 'data/N4214_3band_detects_arbitrary_sub.fits'
obsfile = 'data/fakeobs_1.fits'

#project = 'mf08_225_336_438_rv31'
#project = 'mf09_814_110_160_rv31'
project = 'fakeobs_1'

#filters = ['HST_WFC3_F225W', 'HST_WFC3_F336W', 'HST_WFC3_F438W',
#           'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
#filters = ['HST_WFC3_F225W', 'HST_WFC3_F336W', 'HST_WFC3_F438W']
#filters = ['HST_WFC3_F438W', 'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
#filters = ['HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
filters = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

distanceModulus = 27.41

subproject = 'flatten'

isofile             = 'sedfitter/libs/iso.proposal.fits'
#isofile             = 'sedfitter/libs/padova2010.iso.fits'
lnp_outname         = project+ '/lnp.{}.hd5'.format(project)
stat_outname        = project+'/stat.{}.hd5'.format(project)
res_outname         = project+'/res.{}.fits'.format(subproject)
spectral_grid_fname = 'ngc4214.zsmc.spectral.grid.fits'
sed_grid_fname      = 'ngc4214.zsmc.sed.grid.fits'
#spectral_grid_fname = 'sedfitter/libs/kurucz.proposal.spectral.grid.fits'
#sed_grid_fname = 'sedfitter/libs/kurucz.proposal.sed.grid.fits'

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
        self.floorError = 0.05  # constant error term

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

#obs = Data(obsfile, distanceModulus)
obs = Observations(obsfile)
obs.setFilters(filters)
#for k in filters:
#    obs.data.set_alias(k, k.split('_')[-1] + '_VEGA')
#    obs.data.set_alias(k + 'err', k.split('_')[-1] + '_ERR')


#---------------------------------------------------------
# Model definitions                            [sec:model]
#---------------------------------------------------------

# Main ingredients                       [sec:models:main]
#---------------------------------------------------------

extLaw = extinction.RvFbumpLaw()
osl = stellib.Kurucz()
oiso = ezIsoch(isofile)
#oiso = padova2010()

# Grid sampling definition               [sec:models:grid]
#---------------------------------------------------------

# variable to ensure that range is fully covered in using np.arange
__tiny_delta__ = 0.001

ages   = 10 ** np.arange(6., 9. + __tiny_delta__, 0.1)
masses = 10 ** np.arange(-0.1, 2. + __tiny_delta__, 0.01)
Z      = np.asarray([0.004])
avs    = np.arange(0.0, 5.0 + __tiny_delta__, 0.1)
rvs    = np.arange(2.8, 3.4 + __tiny_delta__, 0.1)
#rvs    = np.asarray([3.1])
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

def fit_model_seds_pytables(obs, sedgrid, threshold=-60, outname=lnp_outname):

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



#---------------------------------------------------------
# Outputs                                   [sec:outputs]
#---------------------------------------------------------
def get_expectation(resfname=lnp_outname,
                    sedgrid=None,
                    obs=obs,
                    valuetype='expect',
                    Q='logA logM Z Av Rv f_bump logT logg logL', flatten=False):
    import numexpr
    if valuetype not in ['expect', 'best']:
        raise ValueError("Valuetype must be in [expect, best], got %s" % valuetype)

    # get the nD stellar/dust SED grid
    if type(sedgrid) == str:
        _sedgrid = grid.FileSEDGrid(sedgrid)
    elif sedgrid is None:
        _sedgrid = get_sedgrid()
    else:
        _sedgrid = sedgrid

    if type(Q) == str:
        _Q = Q.split()
    else:
        _Q = Q

    binAx_dict = {}
    for q in _Q:
        span = _sedgrid.grid[q].max() - _sedgrid.grid[q].min()
        bin_edges = np.linspace(_sedgrid.grid[q].min()-(0.05*span), _sedgrid.grid[q].max()+(0.05*span), 20)
        if q=='logM' or q=='logT' or q=='logL':
            binAx_dict[q] = binningAxis.binningAxis(_sedgrid, q, bin_edges, var_spacing='rect', fixed_key='logA')
        else:
            binAx_dict[q] = binningAxis.binningAxis(_sedgrid, q, bin_edges, var_spacing='fixed')
    
    with tables.openFile(resfname) as hdf:
        n_stars = hdf.root._v_nchildren - 2
        # +2 for the peak value and integrated probability
        res = np.empty((n_stars, len(_Q) + 2), dtype=float)
        with progressbar.PBar(n_stars, txt='expectations') as pbar:
            for tn in range(n_stars):
                star_node = hdf.getNode('/star_%d'  % tn)
                lnp = star_node.lnp[:]
                idx = star_node.idx[:]
                prob = numexpr.evaluate('exp(lnp)', local_dict={'lnp': lnp})
                norm = numexpr.evaluate('sum(prob)', local_dict={'prob': prob})
                if valuetype == 'expect':
                    for e, qk in enumerate(_Q):
                        #vals = _sedgrid.grid[qk][idx]
                        #res[tn, e] = numexpr.evaluate('sum(prob/norm*vals)', local_dict={'vals': vals, 'prob': prob, 'norm': norm})
                        res[tn, e] = binAx_dict[qk].project(prob, inds=idx, to_return='expectation', flatten=flatten)
                else:  # best
                    #vals = _sedgrid.grid[idx[lnp.argmax()]]
                    for e, qk in enumerate(_Q):
                        res[tn, e] = binAx_dict[qk].project(prob, inds=idx, to_return='max')
                res[tn, -1] = prob.max() / norm
                res[tn, -2] = binAx_dict[_Q[0]].project(prob, inds=idx, to_return='total', flatten=False)
                pbar.update(tn, force=True)

    d = {}
    for e, k in enumerate(_Q):
        d[k] = res[:, e]

    d['Pmax'] = res[:, -1]
    d['Ptot'] = res[:, -2]
    t = AstroTable(d)  # add RA, DEC tools
    #adding original values
    for k in obs.data.keys()[::-1]:
        try:
            t.addCol(k, obs.data[k], position=0)
        except:
            pass
    return t


#---------------------------------------------------------
# figures                                   [sec:figures]
#---------------------------------------------------------
def diag_figs(t, Q='logA logM Av Rv logT logL log10(Pmax)', figs=[]):

    import pylab as plt

    def ezrc(fontSize=22., lineWidth=2., labelsize=None):
        """
        ezrc - Define params to make pretty fig
        """
        if labelsize is None:
            labelsize = fontSize + 5
        from pylab import rc, rcParams
        rc('figure', figsize=(8, 6))
        rc('lines', linewidth=lineWidth)
        rc('font', size=fontSize, family='serif', weight='small')
        rc('axes', linewidth=lineWidth, labelsize=labelsize)
        rc('legend', borderpad=0.1, markerscale=1., fancybox=False)
        rc('text', usetex=True)
        rc('image', aspect='auto')
        rc('ps', useafm=True, fonttype=3)
        rcParams['xtick.major.size'] = 10
        rcParams['xtick.minor.size'] = 5
        rcParams['ytick.major.size'] = 10
        rcParams['ytick.minor.size'] = 5
        rcParams['font.sans-serif'] = 'Helvetica'
        rcParams['font.serif'] = 'Helvetica'
        rcParams['text.latex.preamble'] = '\usepackage{pslatex}'

    def colorify(data, vmin=None, vmax=None, cmap=plt.cm.Spectral):
        """ Associate a color map to a quantity vector """
        import matplotlib.colors as colors

        _vmin = vmin or min(data)
        _vmax = vmax or max(data)
        cNorm = colors.normalize(vmin=_vmin, vmax=_vmax)

        scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=cmap)
        colors = map(scalarMap.to_rgba, data)
        return colors, scalarMap, cNorm

    def getcmap(q, default=plt.cm.jet_r):
            _map = { 'logT': plt.cm.jet,
                    'logA': plt.cm.jet,
                    #'logM':
                    'Av': plt.cm.YlOrRd,
                    'Rv': plt.cm.YlOrRd,
                    #'logL':
                    #'log10(Pmax)':
                    }
            return _map.get(q, default)

    def gauss_kern(size, sizey=None):
        """ Returns a normalized 2D gauss kernel array for convolutions """
        size = int(size)
        if not sizey:
            sizey = size
        else:
            sizey = int(sizey)
        x, y = np.mgrid[-size: size + 1, -sizey: sizey + 1]
        g = np.exp( -(x ** 2 / float(size) + y ** 2 / float(sizey)))
        return g / g.sum()

    def cmdfig(x, y, c, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()
        _x = t.evalexpr(x)
        _y = t.evalexpr(y)
        _c = t.evalexpr(c)
        idx = (_x < 20.) & (_x > -20.) & (_y < 50)
        sc = ax.scatter(_x[idx], _y[idx], c=_c[idx], cmap=getcmap(c), **kwargs)
        cb = plt.colorbar(sc, ax=ax)
        cb.set_label(c.replace('_', ' '))
        ax.set_xlabel(x.replace('_VEGA', ' '))
        ax.set_ylabel(y.replace('_VEGA', ' '))
        ax.set_ylim(ax.get_ylim()[::-1])

    def radec_kpc_fig(ra, dec, c, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()
        _c = t.evalexpr(c)
        ra0 = ra.mean()
        dec0 = dec.mean()
        x = obs.distance * 10. * np.tan( np.deg2ra(ra - ra0) )
        y = obs.distance * 10. * np.tan( np.deg2ra(dec - dec0) )
        sc = ax.scatter(x, y, c=_c, cmap=getcmap(c), **kwargs)
        cb = plt.colorbar(sc, ax=ax)
        cb.set_label(c.replace('_', ' '))
        ax.set_xlabel('dx [kpc]')
        ax.set_ylabel('dy [kpc]')

    def radecfig(ra, dec, c, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()
        _c = t.evalexpr(c)
        sc = ax.scatter(ra, dec, c=_c, cmap=getcmap(c), **kwargs)
        cb = plt.colorbar(sc, ax=ax)
        cb.set_label(c.replace('_', ' '))
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')

    def radecim(ra, dec, c, res=4.17e-4, margin=2, reducefn=np.sum, ax=None, **kwargs):
        #res = binning in deg/pix
        bx = np.arange(ra.min() - (margin + 0.5) * res, ra.max() + (margin + 0.5) * res, res )
        by = np.arange(dec.min() - (margin + 0.5) * res, dec.max() + (margin + 0.5) * res, res )
        gsize = int((ra.max() - ra.min())/res+1)
        if ax is None:
            ax = plt.gca()
        _c = t.evalexpr(c)
        sc = ax.hexbin(ra, dec, C=_c, gridsize=gsize, reduce_C_function=np.mean, linewidth=0.)
        cb = plt.colorbar(sc, ax=ax)
        cb.set_label(c.replace('_', ' '))
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')

    def sedsfig(c, ax=None, step=10, **kwargs):
        filts = [ k.split('_')[-1] for k in obs.filters]
        lambs = np.array([ float(k[1: -1]) if (float(k[1: -1]) > 200.) else float(k[1: -1]) * 10 for k in filts])

        select = [ k + '_VEGA' for k in filts ]
        st = np.array(t[select].tolist())[::step]

        _c = t.evalexpr(c)
        _colors, cmap, cnorm = colorify(_c, cmap=getcmap(c))
        _colors = _colors[::step]

        if ax is None:
            ax = plt.gca()
        for k in range(len(st)):
            ax.plot( np.log10(lambs) + 1., -0.4 * st[k], '-', color=_colors[k], alpha=0.5)
        ylim = ax.get_ylim()
        ax.set_xlabel('log(Wavelength/AA)')
        ax.set_ylabel('log(Flux/Vega)')
        ax.set_xlim(2.3, 3.3)
        ax.set_ylim(ylim)
        scalarcmap = ax.pcolorfast(np.asarray(_colors), norm=cnorm)
        cb = plt.colorbar(scalarcmap, ax=ax)
        cb.set_label(c.replace('_VEGA', '').replace('_', ' '))

    if type(Q) == str:
        _Q = Q.split()
    else:
        _Q = Q

    ezrc(14, 1, 15)
    shapes = { 1: (1, 1), 2: (1, 2), 3: (1, 3), 4: (2, 2), 5: (2, 3), 6: (2, 3), 7: (3, 3), 8: (3, 3), 9: (3, 3) }
    _shape = shapes[len(_Q)]

    s = 8

    _figs = []
    if 0 in figs:
        _figs.append(plt.figure(figsize=(_shape[1] * 5, _shape[0] * 4)))
        _axes = [ plt.subplot(_shape[0], _shape[1], j + 1) for j in range(len(_Q)) ]
        cmd_x = 'F438W_VEGA-F814W_VEGA'
        cmd_y = 'F814W_VEGA'

        for e, k in enumerate(_Q):
            cmdfig( cmd_x, cmd_y, k, ax=_axes[e], edgecolor='None', s=s )

    if 1 in figs:
        _figs.append(plt.figure(figsize=(_shape[1] * 5, _shape[0] * 4)))
        _axes = [ plt.subplot(_shape[0], _shape[1], j + 1) for j in range(len(_Q)) ]
        cmd_x = 'F336W_VEGA-F438W_VEGA'
        cmd_y = 'F438W_VEGA'

        for e, k in enumerate(_Q):
            cmdfig( cmd_x, cmd_y, k, ax=_axes[e], edgecolor='None', s=s )

    if 2 in figs:
        _figs.append(plt.figure(figsize=(_shape[1] * 5, _shape[0] * 4)))
        _axes = [ plt.subplot(_shape[0], _shape[1], j + 1) for j in range(len(_Q)) ]
        cmd_x = 'F814W_VEGA-F160W_VEGA'
        cmd_y = 'F814W_VEGA'

        for e, k in enumerate(_Q):
            cmdfig( cmd_x, cmd_y, k, ax=_axes[e], edgecolor='None', s=s )

    if 10 in figs:
        _figs.append(plt.figure(figsize=(_shape[1] * 5, _shape[0] * 4)))
        _axes = [ plt.subplot(_shape[0], _shape[1], j + 1) for j in range(len(_Q)) ]
        for e, k in enumerate(_Q):
            radecfig( t.getRA(), t.getDEC(), k, ax=_axes[e], s=s, edgecolor='None' )

    if 11 in figs:
        _figs.append(plt.figure(figsize=(_shape[1] * 5, _shape[0] * 4)))
        _axes = [ plt.subplot(_shape[0], _shape[1], j + 1) for j in range(len(_Q)) ]
        for e, k in enumerate(_Q):
            radec_kpc_fig( t.getRA(), t.getDEC(), k, ax=_axes[e], s=s, edgecolor='None' )

    if 12 in figs:
        _figs.append(plt.figure(figsize=(_shape[1] * 5, _shape[0] * 4)))
        _axes = [ plt.subplot(_shape[0], _shape[1], j + 1) for j in range(len(_Q)) ]
        for e, k in enumerate(_Q):
            radecim( t.getRA(), t.getDEC(), k, ax=_axes[e], s=s, edgecolor='None' )

    if 20 in figs:
            skip = ['log10(Pmax)']
            _figs.append(plt.figure(figsize=(_shape[1] * 5, _shape[0] * 4)))
            _axes = [ plt.subplot(_shape[0], _shape[1], j + 1) for j in range(len(_Q)) if _Q[j] not in skip ]
            for e, k in enumerate(_Q):
                if k not in skip:
                    sedsfig(k, ax=_axes[e])

    return _figs


#---------------------------------------------------------
# pipeline                                  [sec:pipeline]
#---------------------------------------------------------
def gen_log_description():
    from time import ctime
    with open(project + '.log', 'w') as log:
        log.write('SED fitter log\n')
        log.write('==============\n')

        log.write('\n')
        log.write('Project\t {}\n'.format(project))
        log.write('Started at {}\n'.format(ctime()))

        log.write('\n')
        log.write('Data file\t {}\n'.format(obsfile))
        log.write('\tnObs \t {}\n'.format(len(obs)))
        log.write('\tDistance Modulus\t {}\n'.format(distanceModulus))
        log.write('\tbad values \t >{}\n'.format(obs.badvalue))
        log.write('\tmin Error\t {}\n'.format(obs.minError))
        log.write('\tfloor Error\t {}\n'.format(obs.floorError))
        log.write('\tFilters\t {}\n'.format(filters))
        log.write('\tMapping:\n')
        for k, v in obs.data._aliases.iteritems():
            log.write('\t\t{:<18} --> {:<18}\n'.format(k, v))

        log.write('\n')
        log.write('Spectral Grid\t {}\n'.format(spectral_grid_fname))
        log.write('\t Stellar Library\t {}\n'.format(osl.name))
        log.write('\t Isochrone source\t {}\n'.format(isofile))
        log.write('\t Ages\t from {} to {} by {}\n'.format(np.log10(ages.min()), np.log10(ages.max()), np.log10(ages[1]) - np.log10(ages[0])))
        log.write('\t Masses\t from {} to {} by {}\n'.format(np.log10(masses.min()), np.log10(masses.max()), np.log10(masses[1]) - np.log10(masses[0])))
        log.write('\t Z\t from {} to {} by{}\n'.format(Z.min(), Z.max(), Z[1] - Z[0] if (len(Z) > 1) else 0))

        log.write('\n')
        log.write('SED grid\t {}\n'.format(sed_grid_fname))
        log.write('\t Extinction Law\t {}\n'.format(extLaw.name))
        log.write('\t Av\t from {} to {} by {}\n'.format(avs.min(), avs.max(), avs[1] - avs[0] if (len(avs) > 1) else 0))
        log.write('\t Rv\t from {} to {} by {}\n'.format(rvs.min(), rvs.max(), rvs[1] - rvs[0] if (len(rvs) > 1) else 0))
        log.write('\t Fbump\t from {} to {} by {}\n'.format(fbumps.min(), fbumps.max(), fbumps[1] - fbumps[0] if (len(fbumps) > 1) else 0 ))

        log.write('\n')
        log.write('Likelihood output\t {}\n'.format(lnp_outname))
        log.write('\n')
        log.write('Expectations table\t {}\n'.format(res_outname))


def read_log_description(logfile):

    def skip(f, n):
        for k in range(n):
            f.readline()

    def readline(f, skipempty=True):
        l = f.readline().replace('\n', '')
        if skipempty:
            while (l == ''):
                l = f.readline().replace('\n', '')
        return l

    with open(logfile, 'r') as log:

        skip(log, 3)
        project   = readline(log).split('\t')[-1].replace(' ', '')
        starttime = readline(log).split('\t')[-1][1:]
        obsfile   = readline(log).split('\t')[-1].replace(' ', '')
        nObs      = int(readline(log).split('\t')[-1].replace(' ', ''))
        dMod      = float(readline(log).split('\t')[-1].replace(' ', ''))
        badVal    = readline(log).split('\t')[-1].replace(' ', '')
        minErr    = float(readline(log).split('\t')[-1].replace(' ', ''))
        floorErr  = float(readline(log).split('\t')[-1].replace(' ', ''))
        filters   = readline(log).split('\t')[-1].replace(' ', '').replace("'", '')[1:-1].split(',')
        skip(log, 1)
        l = readline(log)
        maps      = []
        while (l != ''):
            maps.append(l.split('\t')[-1].replace(' ', '').replace('-->', ' ').split())
            l = readline(log, skipempty=False)

        spec_grid_fname = readline(log).split('\t')[-1].replace(' ', '')
        oslname         = readline(log).split('\t')[-1].replace(' ', '')
        isofname        = readline(log).split('\t')[-1].replace(' ', '')
        agesampling     = readline(log).split('\t')[-1].split()
        masssampling    = readline(log).split('\t')[-1].split()
        Zsampling       = readline(log).split('\t')[-1].split()
        sed_grid_fname  = readline(log).split('\t')[-1].replace(' ', '')
        extLaw_name     = readline(log).split('\t')[-1].replace(' ', '')
        avsampling      = readline(log).split('\t')[-1].split()
        rvsampling      = readline(log).split('\t')[-1].split()
        fbsampling      = readline(log).split('\t')[-1].split()

        lnp_fname       = readline(log).split('\t')[-1].replace(' ', '')
        res_fname       = readline(log).split('\t')[-1].replace(' ', '')

        return dict( project=project, starttime=starttime, obsfile=obsfile,
                    nObs=nObs, dMod=dMod, badVal=badVal, minErr=minErr,
                    floorErr=floorErr, filters=filters,
                    spec_grid_fname=spec_grid_fname, oslname=oslname,
                    isofname=isofname, agesampling=agesampling,
                    masssampling=masssampling, Zsampling=Zsampling,
                    sed_grid_fname=sed_grid_fname, avsampling=avsampling,
                    rvsampling=rvsampling, fbsampling=fbsampling,
                    lnp_fname=lnp_fname, res_fname=res_fname, maps=maps,
                    extLaw_name=extLaw_name)


def do_fit():
    sedgrid = get_sedgrid()
    fit_model_seds_pytables(obs, sedgrid, threshold=-40, outname=lnp_outname)


def do_expectations(Q='Av Rv logT logL', flatten=False):
    with RequiredFile(lnp_outname,
                    do_fit) as __lnp_outname__:
        t = get_expectation(resfname=__lnp_outname__, valuetype='expect', Q=Q, flatten=flatten)
        t.write(res_outname, clobber=True, append=False)


def do_figs(project=project, res_outname=res_outname, figs=range(10), Q='logA logM Av logT logL log10(Pmax)', fmt='png', **kwargs):
    if project is None:
        project = res_outname.split('.')[-2]
    if type(fmt) == str:
        _fmt = [fmt]
    else:
        _fmt = fmt

    t = AstroTable(res_outname)
    _figs = diag_figs(t, figs=figs, Q=Q)
    print 'saving figures'
    for e, k in enumerate(_figs):
        for fk in _fmt:
            k.savefig('{}/{}_{}.{}'.format(project, project, e, fk))
    return _figs


def do_figs_previous_project(logfile, figs=range(10) + [20], Q='logA logM Av logT logL log10(Pmax)', ROOT='products/', fmt='png'):
    if ROOT[-1] != '/':
        ROOT += '/'

    log = read_log_description(logfile)

    res_outname = ROOT + log['res_fname']
    project = log['project']
    _figs = do_figs(res_outname=res_outname, figs=figs, Q=Q, project=project, fmt=fmt)

    if __FIG__:
            androfig.__figlist__ = list(_figs)
            androfig.show()

def do_fake_figs(inpath=obsfile,outpaths=[res_outname],labels=[project],Q='Av Rv logM logA',show=True, save=None):
    import fake_figs
    fake_figs.plot_keys(Q, inpath, outpaths, labels, show=show)
    if save is not None:
        fake_figs.plt.savefig(save)

def run():
    gen_log_description()
    if not __FIG__:
        do_expectations()
    else:
        with RequiredFile(res_outname, do_expectations) as __res_outname__:
            Q = 'logA logM Rv Av logT logL log10(Pmax)'
            androfig.__figlist__ = list( do_figs(res_outname=__res_outname__, figs=range(10) + [20], Q=Q, project=project))
            androfig.show()

if __name__ == '__main__':
    #do_fit()
    do_expectations(Q = 'logL', flatten=True)
    #do_fake_figs( Q= 'Av Rv f_bump logL logT', outpaths=[project+'/res.no_flatten.fits',project+'/res.flatten.fits'], labels=['no_flatten', 'flatten'])
    #do_figs(Q='logT logL Av Rv 1/(1+1/Ptot)')
