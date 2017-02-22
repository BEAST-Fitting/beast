#!/usr/bin/env python
"""
Show diagnostic plots from BEAST runs
Meant as a quick check that the results are reasonable

.. history::
    Written 21 Dec 2015 by Karl D. Gordon
    Revised 21 Feb 2017 by Meredith J. Durbin
"""

from __future__ import print_function

import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.colors import LogNorm

# local imports
from beastplotlib import fancify_colname, initialize_parser

def plot_2dhist(xcol, ycol, fig, ax, bins=51):
    '''Updates a given figure/axis object in-place with a new 2D histogram.

    Parameters
    ----------
    xcol : array-like
        Column with x-values
    ycol : array-like
        Column with y-values
    fig : matplotlib figure object
        Figure to plot histogram on
    fig : matplotlib figure object
        Figure to plot histogram on
    bins : int or array-like, optional
        Bins to pass to 2D histogram.

    Returns
    -------
    Nothing.
    '''
    hist = ax.hist2d(xcol, ycol, bins=bins, norm=LogNorm())
    cbar = fig.colorbar(hist[-1], label='Number of stars', ax=ax)
    ax.set_xlabel(fancify_colname(xcol.name))
    ax.set_ylabel(fancify_colname(ycol.name))
    if xcol.name.startswith('logT'):
        ax.invert_xaxis()
    if ycol.name.startswith('logT'):
        ax.invert_yaxis()

def make_diagnostic_plots(statsfile, suffix='Exp'):
    '''Makes a set of 6 diagnostic 2D histograms for BEAST output.

    Parameters
    ----------
    statsfile : str
        Path to file with BEAST output.
    suffix : str, optional
        Column type ('Exp', 'Best', 'p16', 'p50', 'p84')

    Returns
    -------
    fig : matplotlib figure object
        Figure with diagnostic plots
    '''
    stats = Table.read(statsfile)
    base_cnames = ['logT', 'logL', 'Av', 'Rv', 'f_A']
    cnames = ['{}_{}'.format(n, suffix) for n in base_cnames]
    fig, ax = plt.subplots(2,3,figsize=(15,8))
    plot_2dhist(stats[cnames[0]], stats[cnames[1]], fig, ax[0,0])
    plot_2dhist(stats[cnames[0]], stats['chi2min'], fig, ax[0,1])
    plot_2dhist(stats[cnames[2]], stats[cnames[3]], fig, ax[0,2])
    plot_2dhist(stats[cnames[0]], stats[cnames[2]], fig, ax[1,0])
    plot_2dhist(stats[cnames[1]], stats[cnames[2]], fig, ax[1,1])
    plot_2dhist(stats[cnames[2]], stats[cnames[4]], fig, ax[1,2])
    fig.tight_layout()
    return fig

if __name__ == '__main__':
    parser = initialize_parser()
    parser.add_argument('filename', type=str,
                        help='Path to FITS file with output stats')
    suffixes = ['Exp', 'Best', 'p16', 'p50', 'p84']
    parser.add_argument('--suffix', action='store', default='Exp', choices=suffixes,
                        help='Choose column type to plot. \
                              Must be one of: {}'.format(', '.join([s for s in suffixes]))
                        )
    args = parser.parse_args()
    if args.tex:
        plt.rc({'usetex':True})
    basename = args.filename.replace('.fits', '_diagnostics')
    fig = make_diagnostic_plots(args.filename, suffix=args.suffix)
    if args.savefig:
        fig.savefig('{}.{}'.format(basename, args.savefig))
    else:
        plt.show()
