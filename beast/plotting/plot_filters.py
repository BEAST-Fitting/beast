#!/usr/bin/env python
""" Make a nice plot of the filter response functions

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import argparse
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt

from ..observationmodel import phot
from .beastplotlib import initialize_parser


def plot_filters(args, filter_names, save_name='beast_filters',
                 xlim=[0.19, 2.0], ylim=[1e-8, 2e1]):
    '''Plots transmission curves in log-log space.

    Parameters
    ----------
    args : argparse parser object
        Command line arguments
    filter_names : list of str
        List of full names of filters to plot
    save_name : str, optional
        Filename to save plot as
    xlim : length 2 list
        Values to set plot x-limits to
    ylim : length 2 list
        Values to set plot y-limits to

    Returns
    -------
    Nothing.
    '''
    fig, ax = plt.subplots(1, 1, figsize=(10,6))

    # wavelength grid in angstroms for response functions
    waves = np.logspace(3, np.log10(3e4), 501)

    # read in the filter response functions
    flist = phot.load_filters(filter_names, interp=True, lamb=waves)

    color_indices = np.log10(np.array(sorted([f.lam_eff for f in flist])))
    color_indices -= color_indices.min()
    color_indices /= color_indices.max()
    try:
        cmap = mpl.cm.plasma
    except AttributeError: # for matplotlib version < 1.5
        cmap = mpl.cm.rainbow
    try:
        ax.set_prop_cycle(color=[cmap(i) for i in color_indices])
    except AttributeError: # for matplotlib version < 1.5
        ax.set_color_cycle(cmap(i) for i in color_indices)

    for f in flist:
        ax.plot(f.wavelength[f.nonzero], f.transmit[f.nonzero])
        ax.text(f.lam_eff/1e4, 0.1*ylim[-1], f.name.split('_')[-1], ha='center')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('$\lambda$ [$\mu m$]')
    ax.set_ylabel('$B_i(\lambda)$')

    ax.set_xticks([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    fig.tight_layout()

    if args.tex:
        plt.rc({'usetex':True})
    if args.savefig:
        fig.savefig('{}.{}'.format(save_name, args.savefig))
    else:
        plt.show()

if __name__ == '__main__':
    parser = initialize_parser()
    args = parser.parse_args()

    filter_names = ['HST_WFC3_F225W', 'HST_WFC3_F275W', 'HST_WFC3_F336W', 
                    'HST_ACS_WFC_F475W', 'HST_ACS_WFC_F550M',
                    'HST_ACS_WFC_F814W',
                    'HST_WFC3_F110W', 'HST_WFC3_F160W']

    plot_filters(args, filter_names)

