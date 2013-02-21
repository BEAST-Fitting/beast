import numpy as np
import pylab as plt
from sedfitter import eztables


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


def steppify(x, y):
    """ Steppify a curve (x,y). Useful for manually filling histograms """
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


def plot_keys(keys, outdir, outnames, show=True):
    """ Plot average offset per unique values in keys
    INPUTS:
        keys        list    columns to process
        outdir      str     directory where to find the summary fits tables
        outnames    list    list of names for tables

    TODO: replace inputs by a list of files and list of keys
    """
    _outdir = outdir
    if _outdir[-1] != '/':
        _outdir += '/'

    if type(keys) == str:
        _keys = keys.split()
    else:
        _keys = keys

    #prepare figure output, one panel per key
    ezrc(16, 2, 18)
    shapes = { 2: (1, 2), 3: (1, 3), 4: (2, 2), 5: (2, 3), 6: (2, 3), 7: (3, 3), 8: (3, 3), 9: (3, 3) }
    _shape = shapes[len(_keys)]
    plt.figure(figsize=(_shape[1] * 5, _shape[0] * 4))
    _axes = [ plt.subplot(_shape[0], _shape[1], j + 1) for j in range(len(_keys)) ]

    for i in range(len(outnames)):
        summary_table = eztables.Table(_outdir + 'summary_' + outnames[i] + '.fits')
        for j, key in enumerate(_keys):
            #rec_vals = summary_table.data[key+'_recovered']  # recovered values
            true_vals = summary_table.data[key]              # true values
            #rec_vals = rec_vals - true_vals                  # offset
            rec_vals = summary_table.evalexpr('{}_recovered - {}'.format(key, key))  # offset
            uniq_vals = np.unique(true_vals)  # unique true values
            avg_offs = np.zeros(uniq_vals.size)  # Mean of recovered params for given input param
            for k in range(avg_offs.size):
                sel = np.where(true_vals == uniq_vals[k])
                avg_offs[k] = rec_vals[sel].mean()
            #figure
            _axes[j].plot(uniq_vals, avg_offs, label=(outnames[i]).replace('_', ' '), marker='+', markersize=10)

    for j, key in enumerate(_keys):
        _axes[j].set_xlabel(key.replace('_', ''))
        _axes[j].set_ylabel(key.replace('_', '') + ' (out - in)')
        if j == 0:
            l = _axes[j].legend(loc='best')
            l.draw_frame(False)
        _axes[j].plot(_axes[j].get_xlim(), [0, 0], linestyle='dotted', color='black')
    plt.subplots_adjust(wspace=0.3)
    if show:
        plt.show()
