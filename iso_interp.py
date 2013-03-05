# -*- coding: utf-8 -*-

import os
import sys
sys.path.append(os.getenv('HOME') + '/bin/python/libs')

import pylab as plt
import numpy as np
from scipy import interpolate
from numpy import log10

import isochrone
import eztables


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
            return self.map(self.func, *args, **kwargs)

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


def isointerp(pk, iso, dloga=0.1, dlogm=0.3, fZ=10., predict='logA logT logg logM logL Z'.split()):
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
    a_lim = [log10(min(iso.ages)), log10(max(iso.ages))]
    m_lim = [min(t['logM']), max(t['logM'])]
    z_lim = [min(iso.Z), max(iso.Z)]
    if ((ak > a_lim[1]) | (ak < a_lim[0])):
        print 'log(age) {} outside {}'.format(ak, a_lim)
        return [np.nan] * len(predict)
    if ((zk > z_lim[1]) | (zk < z_lim[0])):
        print 'Z {} outside {}'.format(zk, z_lim)
        return [np.nan] * len(predict)
    if ((mk > m_lim[1]) | (mk < m_lim[0])):
        print 'log(mass) {} outside {}'.format(mk, m_lim)
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
            print 'empty Ball around pk = {}'.format(pk)
            return [np.nan] * len(predict)


def isointerp1(pk, iso, dloga=0.1, dlogm=0.3, fZ=10., predict='logA logT logg logM logL Z'.split()):
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


    New idea:
        * ak, mk, zk = pk
        * find i / a[i] <= ak < a[i+1]
        * find j / z[i] <= zk < a[i+1]

    """
    ak, mk, zk = pk

    t = eztables.Table(iso.data)

    # Check isochrones cover pk
    a_lim = [log10(min(iso.ages)), log10(max(iso.ages))]
    m_lim = [min(t['logM']), max(t['logM'])]
    z_lim = [min(iso.Z), max(iso.Z)]
    if ((ak > a_lim[1]) | (ak < a_lim[0])):
        print 'log(age) {} outside {}'.format(ak, a_lim)
        return [np.nan] * len(predict)
    if ((zk > z_lim[1]) | (zk < z_lim[0])):
        print 'Z {} outside {}'.format(zk, z_lim)
        return [np.nan] * len(predict)
    if ((mk > m_lim[1]) | (mk < m_lim[0])):
        print 'log(mass) {} outside {}'.format(mk, m_lim)
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
            print 'empty Ball around pk = {}'.format(pk)
            return [np.nan] * len(predict)


def doplots(t, r):
    slides(20, 2)
    plt.figure(figsize=(10, 10))
    ax = plt.subplot(221)
    ax.plot(t['logT'], t['logg'], 'k,', alpha=0.1, zorder=0)
    ax.scatter(r['logT'], r['logg'], c=r['logA'], edgecolor='None', s=20, zorder=10)
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_xlabel('log(T)')
    ax.set_ylabel('log(g)')
    ax.set_title('log(T) - log(g) - log(A)')

    ax1 = plt.subplot(222, sharex=ax)
    ax1.plot(t['logT'], t['logL'], 'k,', alpha=0.1, zorder=0)
    ax1.scatter(r['logT'], r['logL'], c=r['logA'], edgecolor='None', s=20, zorder=10)
    ax1.set_xlabel('log(T)')
    ax1.set_ylabel('log(L)')
    ax1.set_title('log(T) - log(L) - log(t)')

    ax2 = plt.subplot(223)
    ax2.plot(t['logA'], t['logM'], 'k,', label='origin', alpha=0.2, zorder=0)
    ax2.plot(r['logA'], r['logM'], ',', label='interp', zorder=10)
    ax2.set_xlabel('log(t)')
    ax2.set_ylabel('log(M)')

    ax3 = plt.subplot(224, sharex=ax1, sharey=ax1)
    ax3.plot(t['logT'], t['logL'], 'k,', alpha=0.1, zorder=0)
    ax3.scatter(r['logT'], r['logL'], c=r['logM'], edgecolor='None', s=20, zorder=10)
    ax3.set_xlabel('log(T)')
    ax3.set_ylabel('log(L)')
    ax3.set_title('log(T) - log(L) - log(M)')

    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()


# get the isochrones
iso = isochrone.padova2010()

# Define the grid

__tiny_delta__ = 0.001
#ages   = 10 ** np.linspace(log10(iso.ages.min()), log10(iso.ages.max()), 11)
ages   = 10 ** np.arange(log10(iso.ages.min()), log10(iso.ages.max()) + __tiny_delta__, 0.05)
#masses = 10 ** np.linspace(iso.data['logM'].min(), iso.data['logM'].max(), 100)
#Z      = (0.02, 0.007, 0.004)

#ages   = 10 ** np.linspace(6, 7, 11)
#masses = 10 ** np.linspace(-0.8, 1.99, 200)
masses = 10 ** np.arange(-0.8, 1.99 + __tiny_delta__, 0.001)
Z      = (0.02, 0.004)
# grid spacing for dust
# variable to ensure that range is fully covered in using numpy.arange
#avs            = numpy.arange(0.0, 5.0 + __tiny_delta__, 0.1)
#rvs            = numpy.arange(1.0, 6.0 + __tiny_delta__, 0.5)
#fbumps         = numpy.arange(0.0, 1. + __tiny_delta__, 0.1)


# Do the actual isochrone interpolation:
# i.e get 'logA logT logg logM logL Z' from the initial isochrone grid
predict = 'logA logT logg logM logL Z'.split()

# generate ND grid points from age, mass and Z
_ages, _masses, _Z = np.ix_(log10(ages), log10(masses), Z)
itpts = np.nditer([_ages, _masses, _Z])
# needs to keep orders => using eztables
t = eztables.Table(iso.data)
dt = t[predict][0].dtype

ndata = len(ages) * len(Z) * len(masses)
pts = np.asarray([ k for k in itpts ])


def func(pk):
    return isointerp(pk, iso)
parafunc = parallelTask(func, nthreads=25)


def working_3dinterp():
    # partial function that will be parallelized

    from decorators import timeit

    with timeit('Interpolation of {} points using {} threads'.format(ndata, parafunc.nthreads)):
        r = np.rec.fromrecords(parafunc(pts), dtype=dt)
    return r

#r = working_3dinterp()
#doplots(t, r)
