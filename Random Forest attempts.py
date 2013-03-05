# -*- coding: utf-8 -*-

import os
import sys
sys.path.append(os.getenv('HOME') + '/bin/python/libs')

import pylab as plt
import numpy as np
from scipy import interpolate
from numpy import log10

import grid
import isochrone
import stellib
import eztables

from sklearn.ensemble import RandomForestClassifier

#===============================================================================
#============== MAP GENERALIZATION TO PARALLEL PROCESSING ======================
# TODO:
#   find a way to make a decorator without pickling issues
#===============================================================================
from functools import update_wrapper


class parallelTask(object):
    def __init__(self, func=None, nthreads=0, pool=None, managed={}):
        self.nthreads = nthreads
        self.__set_pool__(pool)
        self.func = func
        if func is not None:
            update_wrapper(self, self.func)

    def __set_pool__(self, pool):
        if (self.nthreads != 0) & (pool is None):
            from multiprocessing import Pool, cpu_count
            if self.nthreads < 0:
                self.nthreads = cpu_count()
            self.pool = Pool(self.nthreads)
            self.map = self.pool.map
        elif pool is not None:
            self.map = pool.map
        else:
            self.map = map

    def __call__(self, *args, **kwargs):
        if self.func is None:
            self.func = args[0]
            update_wrapper(self, self.func)
            return self
        else:
            return self.map(func, *args, **kwargs)

#===============================================================================
#============== FIGURE SETUP FUNCTIONS =========================================
#===============================================================================


def slides(fontSize=22., lineWidth=2.):
    """
    slides - Define params to make pretty fig for slides
    """
    from pylab import rc, rcParams
    rc('figure', figsize=(8, 6))
    rc('lines', linewidth=lineWidth)
    rc('font', size=fontSize, family='serif', weight='small')
    rc('axes', linewidth=lineWidth, labelsize=fontSize + 5)
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


def theme(ax=None, minorticks=False):
    """ update plot to make it nice and uniform """
    from matplotlib.ticker import AutoMinorLocator
    from pylab import rcParams, gca, tick_params
    if minorticks:
        if ax is None:
            ax = gca()
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    tick_params(which='both', width=rcParams['lines.linewidth'])


def steppify(x, y):
    """ Steppify a curve (x,y). Useful for manually filling histograms """
    import numpy as np
    dx = 0.5 * (x[1:] + x[:-1])
    xx = np.zeros( 2 * len(dx), dtype=float)
    yy = np.zeros( 2 * len(y), dtype=float)
    xx[0::2], xx[1::2] = dx, dx
    yy[0::2], yy[1::2] = y, y
    xx = np.concatenate(([x[0] - (dx[0] - x[0])], xx, [x[-1] + (x[-1] - dx[-1])]))
    return xx, yy


def colorify(data, vmin=None, vmax=None, cmap=plt.cm.Spectral):
    """ Associate a color map to a quantity vector """
    import matplotlib.colors as colors

    _vmin = vmin or min(data)
    _vmax = vmax or max(data)
    cNorm = colors.normalize(vmin=_vmin, vmax=_vmax)

    scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=cmap)
    colors = map(scalarMap.to_rgba, data)
    return colors, scalarMap

#===============================================================================
#===============================================================================
#===============================================================================


def plot_importance(g, Forest):
    importances = Forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in Forest.estimators_], axis=0)
    if isinstance(g, grid.FileSEDGrid):
        #x = np.array([275, 336, 475, 814, 1100, 1600], dtype=float)
        lbl = 'F275W F336W F475W F814W F110W F160W'.split()

        plt.figure()
        ax = plt.subplot(111)
        ax.fill_between(log10(g.lamb), importances + std, importances - std, alpha=0.2)
        ax.plot(log10(g.lamb), importances)
        ax.set_xlabel(r'$\log_{10}(\lambda) [\rm{\AA}]$')
        ax.set_ylabel(r'Importance')
        ax1 = plt.twiny(ax)
        ax1.set_xlim(ax.get_xlim())
        plt.xticks(log10(g.lamb), lbl, rotation=45)
        ax1.set_xlabel(r'Filters')
        ax1.grid(linewidth=1, linestyle='-', color='0.', alpha=0.5)
    else:
        plt.figure()
        ax = plt.subplot(111)
        ax.set_title("Wavelength importances")
        #plt.errorbar(log10(g.lamb), importances, yerr=std)
        ax.plot(log10(g.lamb), importances)
        ax.set_xlabel(r'$\log_{10}(\lambda) [\rm{\AA}]$')
        ax.set_ylabel(r'Importance')
        ax.set_ylim(-1e-3, plt.ylim()[1])


def plot_predict(Forest, data_X, data_Y, X, Y, noise=0.001):
    predict_Y = Forest.predict(data_X)
    plt.figure(figsize=(15, 5))
    ax = plt.subplot(1, 2, 1)
    nx = np.random.normal(0, noise, len(predict_Y))
    ny = np.random.normal(0, noise, len(predict_Y))
    ax.scatter(data_Y + nx, predict_Y + ny, s=100, c=(predict_Y - data_Y ))
    lim = ax.get_xlim()
    ax.plot(lim, lim, ls=':')
    ax.set_xlabel('data Y')
    ax.set_ylabel('predicted Y')

    ax = plt.subplot(1, 2, 2)
    ax.scatter(X, Y, c=(predict_Y - data_Y ), s=100)


def plot_seds(g, data_Y):

    #phat filters
    x = np.array([275, 336, 475, 814, 1100, 1600], dtype=float)
    lbl = 'F275W F336W F475W F814W F110W F160W'.split()

    plt.figure()
    if isinstance(g, grid.FileSEDGrid):
        _colors, scalarMap = colorify(data_Y, cmap=plt.cm.RdYlBu_r)

        ax = plt.subplot(111)
        ax.set_color_cycle(_colors)
        ax.semilogy(log10(x), g.seds.T, alpha=0.2)

        ax.set_xlim(log10(250), log10(1700))
        ax.set_ylabel(r'Flux [$\rm{erg/s/\AA}]$')
        ax.set_xlabel(r'$\log_{10}(\lambda) [\rm{\AA}]$')

        ax1 = plt.twiny(ax)
        ax1.set_xlim(250, 1700)
        ax1.set_xscale('log')
        plt.xticks(x, lbl, rotation=45)
        ax1.set_xlabel(r'Filters')
        ax1.grid(linewidth=1, linestyle='-', color='0.', alpha=0.5)
    else:
        _colors, scalarMap = colorify(data_Y, cmap=plt.cm.Reds)

        ax = plt.subplot(111)
        ax.set_color_cycle(_colors)
        ax.loglog(g.lamb, g.seds.T)

        ax.set_xlim(100, 4000)
        ax.set_ylim(1e20, g.seds.max())
        ax.set_ylabel(r'Flux [$\rm{erg/s/\AA}]$')
        ax.set_xlabel(r'$\log_{10}(\lambda) [\rm{\AA}]$')

        ax1 = plt.twiny(ax)
        ax1.set_xlim(ax.get_xlim())
        ax1.set_xscale('log')
        plt.xticks(x, lbl, rotation=45)
        ax1.set_xlabel(r'Filters')
        ax1.grid(linewidth=1, linestyle='-', color='0.', alpha=0.5)


def test1():
    slides(20, 2)  # Look better

    #inputs
    g = grid.FileSEDGrid('kurucz2004.seds.grid.fits')
    #g = grid.FileSpectralGrid('kurucz2004.spectral.grid.fits')
    Forest = RandomForestClassifier(n_estimators=100, compute_importances=True, n_jobs=-1, random_state=0)

    #training
    Forest.verbose = 0
    N = 1000
    data_Y = g.logT[::N]  # g.seds[:,2]/g.seds.max()
    data_X = np.log10(g.seds / g.seds.max())[::N]
    #data_X[:,-1] = 0.
    data_X[np.isinf(data_X) | np.isnan(data_X) ] = 0.

    Forest = Forest.fit(data_X, data_Y)

    #plots
    slides(20, 2)
    plot_importance(g, Forest)
    plot_predict(Forest,
                data_X * np.random.normal(1., 0.01, data_X.shape),  # prediction made on this variable
                data_Y,
                g.logA[::N] * np.random.normal(1., 0.005, len(data_X)),
                g.logM[::N] * np.random.normal(1., 0.005, len(data_X)),
                noise=0.001)
    #plot_seds(g, data_Y)


def test2():
    #inputs
    g = grid.FileSEDGrid('kurucz2004.seds.grid.fits')
    Forest = RandomForestClassifier(n_estimator=100, compute_importances=True, n_jobs=-1, random_state=0)

    #training
    Forest.verbose = 0
    N = 1000
    data_Y = g.logA[::N]
    data_X = 10 ** (-0.4 * g.seds[::N])
    data_X = log10(data_X / data_X.max())
    #data_X[:,-1] = 0.
    data_X[np.isinf(data_X) | np.isnan(data_X) ] = 0.

    Forest = Forest.fit(data_X, data_Y)

    #plots
    slides(20, 2)
    plot_importance(g, Forest)
    plot_predict(Forest,
                data_X * np.random.normal(1., 0.01, data_X.shape),  # prediction made on this variable
                data_Y,
                g.logA[::N] * np.random.normal(1., 0.005, len(data_X)),
                g.logM[::N] * np.random.normal(1., 0.005, len(data_X)),
                noise=0.001)
    #plot_seds(g, data_Y)


def test3():
    slides(20, 2)
    st = stellib.Kurucz()
    iso = isochrone.padova2010()
    g1 = grid.FileSEDGrid('kurucz2004.spectral.grid.fits')

    ax = plt.subplot(111)
    ax.plot(st.Teff, st.logg, '.')
    ax.plot(iso.data['logT'], iso.data['logg'], ',')
    for k in np.unique(g1.logA):
        plt.figure()
        ax = plt.subplot(111)
        ik = iso._get_isochrone(k, metal=0.02)
        ax.plot(ik.data['logT'], ik.data['logg'])
    ax.plot(g1.logT, g1.logg, '.')


def isointerp(pk, iso, dloga=0.1, dlogm=0.2, fZ=5., predict='logA logT logg logM logL Z'.split()):
    """ Interpolate isochrones in the (logA, logM, Z) space
    at a given vector point pk

    INPUTS:
        pk  ndarray[float, ndim=1]  (logA, logM, Z) point
        iso isochrone.Isochrone     isochrone object instance

    KEYWORDS:
        dloga   float   logA range around pk used in the interpolation
        dlogm   float   logM range around pk used in the interpolation
        fZ      float   [zk / fZ, fZ * zk] range around pk used in the interpolation
        predict list    list of isochrone properties returned by the interpolator
    """
    ak, mk, zk = pk

    t = eztables.Table(iso.data)

    # Check isochrones cover pk
    if (ak > max(iso.ages)) | (ak < min(iso.ages)):
                return [np.nan] * len(predict)
    if (zk > max(iso.Z)) | (mk < min(iso.Z)):
                return [np.nan] * len(predict)
    if (mk > max(iso.data['logM'])) | (mk < min(iso.data['logM'])):
                return [np.nan] * len(predict)

    # Compute norm_inf
    d_a = abs(t['logA'] - ak)
    d_m = abs(t['logM'] - mk)
    d_Z = abs(t['Z'] - zk)

    d = np.zeros(d_a.shape, dtype=float)
    for k in (d_a, d_m, d_Z):
        d += k

    # if pk is exactly in the isochrones then return without interpolation
    # else do a ND linear interpolation in (logA, logM, Z)
    if (d.min() == 0.):
        return t[d.argmin()][predict].tolist()
    else:
        # In order to reduce the computations, isochrone points are restricted
        # to a ball/boxel defined by the keyword values
        idx = (d_a <= dloga) & (d_m <= dlogm) & (t['Z'] >= (zk / fZ)) & (t['Z'] <= (zk * fZ))
        if True in idx:
            #restrict to a few points
            data = t.data[idx]
            #interpolation points
            points = np.array([data[k] for k in 'logA logM Z'.split()]).T
            #grid values
            values = np.array([ data[k] for k in predict ]).T
            # Sometimes it does not work, which can be because there are not
            # enought points around pk to keep a correct interpolation.
            # If the point is not in the convex boxel, interpolation returns
            # nans
            try:
                return interpolate.LinearNDInterpolator(points, values)([pk]).flatten().tolist()
            except:
                return [np.nan] * len(predict)
        else:
            return [np.nan] * len(predict)


__tiny_delta__ = 0.001
ages   = 10 ** np.linspace(6, 7, 11)
masses = 10 ** np.linspace(-0.8, 2, 100)
Z      = (0.02, 0.004)
# grid spacing for dust
# variable to ensure that range is fully covered in using numpy.arange
#avs            = numpy.arange(0.0, 5.0 + __tiny_delta__, 0.1)
#rvs            = numpy.arange(1.0, 6.0 + __tiny_delta__, 0.5)
#fbumps         = numpy.arange(0.0, 1. + __tiny_delta__, 0.1)

iso = isochrone.padova2010()

_ages, _masses, _Z = np.ix_(log10(ages), log10(masses), Z)
itpts = np.nditer([_ages, _masses, _Z])
predict = 'logA logT logg logM logL Z'.split()
t = eztables.Table(iso.data)
dt = t[predict][0].dtype

niter = len(ages) * len(Z)
ndata = niter * len(masses)
pts = np.asarray([ k for k in itpts ])


def func(pk):
    return isointerp(pk, iso)

parafunc = parallelTask(func, nthreads=25)
r = np.rec.fromrecords(parafunc(pts), dtype=dt)

"""
plt.figure(figsize=(10, 10))
ax = plt.subplot(221)
ax.plot(t['logT'], t['logL'], ',', alpha=0.5)
idx = (  (t['logA'] <= log10(max(ages)) )  & (t['logA'] >= log10(min(ages)) )
       & (t['logM'] <= log10(max(masses))) & (t['logM'] >= log10(min(masses)))
       & (t['Z']    <= max(Z))      & (t['Z'] >= min(Z) ) )
ax.plot(t['logT'][idx], t['logL'][idx], '.')
xlim = ax.get_xlim()
ylim = ax.get_ylim()

ax = plt.subplot(222)
ax.scatter(r['logT'], r['logL'], c=r['logA'], edgecolor='None', s=20)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax = plt.subplot(223)
ax.scatter(r['logT'], r['logL'], c=r['Z'], edgecolor='None', s=20)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax = plt.subplot(224)
ax.scatter(r['logT'], r['logL'], c=r['logM'], edgecolor='None', s=20)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
"""
