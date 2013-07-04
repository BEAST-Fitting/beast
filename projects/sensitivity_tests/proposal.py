"""
All the figures for proposals
"""

from sedfitter import figure as ezfig
from sedfitter import eztables
from sedfitter import grid
import pylab as plt
import tables
import numpy as np
import numexpr


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

#prepare figure output, one panel per key
ezrc(16, 2, 18)


def get_pdf(g, hdf, key, xp=None):
    """ Build the pdf of a quantity key evaluated at given points xp """

    model_vals = g.grid[key]

    if xp is None:
        xp = 5

    if not hasattr(xp, '__iter__'):
        bins = np.linspace(model_vals.min(), model_vals.max(), xp)
    else:
        bins = xp

    #number of fake stars:
    n_stars = hdf.root._v_nchildren - 1  # - 1 since one node is wavelength
    #number of noise samples per star (assuming all identical)
    m_tests = hdf.root.fakeStar_0._v_nchildren
    results = np.zeros((n_stars * m_tests, 2), dtype=float)
    for i in range(n_stars):
        #get input value
        model_idx = hdf.getNode('/fakeStar_%d' % (i))._v_attrs['fakein']
        input_val = model_vals[model_idx]
        for j in range(m_tests):
            lnp = hdf.getNode('/fakeStar_%d/fake_%d/lnp' % (i, j)).read()
            indx = hdf.getNode('/fakeStar_%d/fake_%d/idx' % (i, j)).read()
            val = model_vals[indx]
            prob = numexpr.evaluate('exp(lnp)', {'lnp': lnp})
            norm = numexpr.evaluate('sum(prob)', {'prob': prob})
            E = numexpr.evaluate('sum(prob/norm * val)', {'prob': prob, 'norm': norm, 'val': val})
            results[2 * i + j, 0] = input_val
            results[2 * i + j, 1] = E

    bx = np.digitize(results[:, 0], bins)
    a = np.asarray([np.mean(results[bx == k], axis=0) for k in set(bx)])

    return a, results


def plot_keys(grid_fname, fname, keys):
    g = grid.FileSEDGrid(grid_fname)

    if type(keys) == str:
        _keys = keys.split()
    else:
        _keys = keys

    if type(fname) == str:
        _fname = fname.split()
    else:
        _fname = fname

    #prepare figure output, one panel per key
    ezrc(16, 2, 18)
    shapes = { 1: (1, 1), 2: (1, 2), 3: (1, 3), 4: (2, 2), 5: (2, 3), 6: (2, 3), 7: (3, 3), 8: (3, 3), 9: (3, 3) }
    _shape = shapes[len(_keys)]
    plt.figure(figsize=(_shape[1] * 5, _shape[0] * 4))
    _axes = [ plt.subplot(_shape[0], _shape[1], j) for j in range(len(_keys)) ]

    for fk in _fname:
        with tables.openFile(fk, 'r') as hdf:
            for e, k in enumerate(_keys):
                print 'plotting %s from %s' % (k, fk)
                a, r = get_pdf(g, hdf, k, 6)
                _axes[e].plot(r[:, 0], r[:, 1] - r[:, 0], '.')
                _axes[e].plot(a[:, 0], a[:, 1] - a[:, 0], 'o-')

    for e, k in enumerate(_keys):
        _axes[e].set_xlabel(k)
        _axes[e].set_ylabel(r'$\Delta {}$'.format(k))
    plt.subplots_adjust(wspace=0.3)
