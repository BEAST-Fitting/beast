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

def plot_generic(xcol, ycol, fig=None, ax=None, plottype='hist', bins=51,
                 thresh_col=None, thresh=0, thresh_op='less', plot_kwargs={}):
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
    thresh_col : astropy table column or pandas series
        Column with values to filter xcol and ycol by; for example, 
        if you want to look at only stars with chi^2 greater than 20 
        you would set thresh_col = table['chi2min'], thresh=20, and 
        thresh_op = 'greater'
    thresh : scalar
        Threshold of values to compare thresh_col to
    thresh_op : string
        Operator comparison to perform between thresh_col and thresh.
        Must be one of Numpy logical operator functions: {'greater', 
        'greater_equal', 'less', 'less_equal', 'equal', 'not_equal'}
    plot_kwargs : dict
        Additional keyword arguments to pass to plotting function.

    Returns
    -------
    Nothing.
    '''
    if fig is None:
        fig = plt.gcf()
    if ax is None:
        ax = plt.gca()
    if thresh_col is not None:
        condition = getattr(np, thresh_op)(thresh_col, thresh)
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
