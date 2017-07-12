#!/usr/bin/env python
"""
Show diagnostic plots for a square region in two parameters
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.colors import LogNorm
from matplotlib import rc
import numpy as np

# local imports
from .beastplotlib import (fancify_colname, initialize_parser,
                          plot_generic, set_params)

def make_region_diag_plots(stats, xparam='logT', yparam='logL',
                           suffix='Exp', pxrange=[3.8,4.0], pyrange=[1.0,6.0],
                           figsize=(15,10)):
    '''Makes a set of 4 diagnostic 2D histograms for BEAST output.

    Parameters
    ----------
    stats : astropy Table
        stats table from BEAST run
    xparam: str, optional
        Parameter (column) name
        Default to 'logT'
    yparam: str, optional
        Parameter (column) name
        Default to 'logL'
    pxrange: float, float, optional
        x range for region
        default to [0.0, 1.0]
    pyrange: float, float, optional
        y range for region
        default to [0.0, 1.0]
    suffix : str, optional
        Column type ('Exp', 'Best', 'p16', 'p50', 'p84').
        Defaults to 'Exp'.
    figsize : tuple of ints, optional
        Size of figure to return in inches (width, height).
        Defaults to (10, 5.5).

    Returns
    -------
    fig : matplotlib figure object
        Figure with diagnostic plots
    '''
    base_cnames = [xparam, yparam]
    cnames = ['{}_{}'.format(n, suffix) for n in base_cnames]

    # cut out the data in the region defined
    indxs, = np.where(np.logical_and(pxrange[0] <= stats[cnames[0]],
                                     stats[cnames[0]] <= pxrange[1]))
    indxs2, = np.where(np.logical_and(pyrange[0] <= stats[cnames[1]][indxs],
                                      stats[cnames[1]][indxs] <= pyrange[1]))
    if len(indxs2) <= 0:
        print('no data in selection window')
        exit()
    else:
        stats_region = stats[indxs[indxs2]]

    base_pnames = ['logT','logL', 'Av', 'Rv']
    pnames = ['{}_{}'.format(n, suffix) for n in base_pnames]
    plot_pairs = [[pnames[0], pnames[1]],
                  [pnames[0], pnames[2]],
                  [pnames[2], pnames[3]],
                  ['RA', 'DEC']]

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    ax = axes.ravel()
    for i, pair in enumerate(plot_pairs):
        plot_generic(stats_region, pair[0], pair[1], fig, ax[i],
                     plot_kwargs={'norm':LogNorm()})
    fig.tight_layout()
    return fig

if __name__ == '__main__':
    parser = initialize_parser()
    parser.add_argument('filename', type=str,
                        help='Path to FITS file with output stats')
    suffixes = ['Exp', 'Best', 'p50']
    parser.add_argument('--suffix', action='store', default='Exp',
                        choices=suffixes,
                        help='Choose column type to plot. \
                        Must be one of: "{}"'.format('", "'.join(suffixes))
                        )
    params = ['logA','M_ini','Z','Av','Rv','f_A','logT','logL']
    parser.add_argument('--xparam', action='store', default='logT',
                        choices=params,
                        help='Choose column for xaxis \
                        Must be one of: "{}"'.format('", "'.join(params))
                        )
    parser.add_argument('--xrange', type=float, nargs=2, default=[0.0,1.0],
                        help='x range for selection')
    parser.add_argument('--yparam', action='store', default='logL',
                        choices=params,
                        help='Choose column for yaxis \
                        Must be one of: "{}"'.format('", "'.join(params))
                        )
    parser.add_argument('--yrange', type=float, nargs=2, default=[0.0,1.0],
                        help='y range for selection')
    args = parser.parse_args()
    if args.tex:
        plt.rc({'usetex':True})
    basename = args.filename.replace('.fits', '_diagnostics')

    set_params(lw=2, fontsize=16, usetex=False)
    stats = Table.read(args.filename)
    fig = make_region_diag_plots(stats, suffix=args.suffix,
                              xparam=args.xparam, yparam=args.yparam,
                              xrange=args.xrange, yrange=args.yrange)

    if args.savefig:
        fig.savefig('{}.{}'.format(basename, args.savefig))
    else:
        plt.show()
