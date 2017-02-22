#!/usr/bin/env python
"""
Library of general functions for the BEAST plotting scripts

.. history::
    Written 21 Feb 2017 by Meredith J. Durbin
"""
import argparse
import numpy as np
import re

def fancify_colname(name):
    '''someone who knows regex pls rewrite
    '''
    ns = name.split('_')
    if re.match(re.compile('log._*'), name) is not None:
        fancyname = '$\log ({})$ ({})'.format(name[3], ns[-1])
    elif re.match(re.compile('.v_*'), name) is not None:
        fancyname = '${}_V$ ({})'.format(name[0], ns[-1])
    elif re.match(re.compile('._._*'), name) is not None:
        fancyname = '${}_'.format(ns[0]) + '{' + '{}'.format(ns[1]) + '}' +\
                    '$ ({})'.format(ns[-1])
    elif name == 'chi2min':
        fancyname = '$\chi^2$'
    else:
        fancyname = ' '.join(ns)
    return r'{}'.format(fancyname)

def initialize_parser():
    '''For running from command line, initialize argparse with common args
    '''
    ftypes = ['png', 'jpg', 'jpeg', 'pdf', 'ps', 'eps', 'rgba',
              'svg', 'tiff', 'tif', 'pgf', 'svgz', 'raw']
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--savefig', action='store', default=False, choices=ftypes,
                        help='Save figure to a file of specified type rather than show it. \
                              Must be one of: "{}"'.format('", "'.join(ftypes))
                        )
    parser.add_argument('-t', '--tex', action='store_true', 
                        help='Configure Matplotlib to use LaTeX; \
                              sets rcParam "usetex" to True. \
                              (Requires working TeX installation.)\
                              Defaults to false.')
    return parser


def plot_generic(xcol, ycol, fig, ax, plottype='hist', bins=51, chi_col=None, 
                 chi_thresh=0, chi_op='less', plot_kwargs={}):
    '''Plot anything from the BEAST output against anything else.

    Parameters
    ----------
    xcol : astropy table column or pandas series
        Column with x-values
    ycol : astropy table column or pandas series
        Column with y-values
    fig : matplotlib figure object
        Figure to plot histogram on
    ax : matplotlib axes object
        Subplot of figure to plot on
    plottype : str
        String indicating what type of plot to make.
        Must be one of: {'hist', 'scatter'}
    bins : int or array-like, optional
        Number of bins to pass to 2D histogram.
        Only used if plottype == 'hist'.
    chi_col : astropy table column or pandas series
        Column with chi^2 values if chi^2 cutoff is desired
    chi_thresh : scalar
        Threshold of chi^2 values to compare chi_col to
    chi_op : string
        Operator comparison to perform between chi_col and chi_thresh.
        Must be one of Numpy logical operator functions: {'greater', 
        'greater_equal', 'less', 'less_equal', 'equal', 'not_equal'}
    plot_kwargs : dict
        Additional keyword arguments to pass to plotting function.

    Returns
    -------
    Nothing.
    '''
    if chi_col is not None:
        condition = getattr(np,chi_op)(chi_col,chi_thresh)
        xcol = xcol[condition]
        ycol = ycol[condition]
    if plottype == 'hist':
        hist = ax.hist2d(xcol, ycol, bins=bins, **plot_kwargs)
        cbar = fig.colorbar(hist[-1], label='Number of stars', ax=ax)
    elif plottype == 'scatter':
        sc = ax.scatter(xcol, ycol, **plot_kwargs)
    ax.set_xlabel(fancify_colname(xcol.name))
    ax.set_ylabel(fancify_colname(ycol.name))
    if xcol.name.startswith('logT'):
        ax.invert_xaxis()
    if ycol.name.startswith('logT'):
        ax.invert_yaxis()
