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
from .beastplotlib import fancify_colname, initialize_parser, plot_generic

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
    plot_pairs = [[cnames[0], cnames[1]], # logL vs logT
                  [cnames[0], 'chi2min'], # chi2 vs logT
                  [cnames[2], cnames[3]], # Rv vs Av
                  [cnames[0], cnames[2]], # Av vs logT
                  [cnames[1], cnames[2]], # Av vs logL
                  [cnames[2], cnames[1]]] # f_A vs Av
    fig, axes = plt.subplots(2,3,figsize=(15,8))
    ax = axes.ravel()
    for i, pair in enumerate(plot_pairs):
        plot_generic(stats[pair[0]], stats[pair[1]], fig, ax[i], 
                     plot_kwargs={'norm':LogNorm()})
    fig.tight_layout()
    return fig

if __name__ == '__main__':
    parser = initialize_parser()
    parser.add_argument('filename', type=str,
                        help='Path to FITS file with output stats')
    suffixes = ['Exp', 'Best', 'p16', 'p50', 'p84']
    parser.add_argument('--suffix', action='store', default='Exp', choices=suffixes,
                        help='Choose column type to plot. \
                              Must be one of: "{}"'.format('", "'.join(suffixes))
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
