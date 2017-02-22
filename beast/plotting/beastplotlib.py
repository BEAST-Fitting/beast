#!/usr/bin/env python
"""
Library of general functions for the BEAST plotting scripts

.. history::
    Written 21 Feb 2017 by Meredith J. Durbin
"""
import argparse
import numpy as np

def fancify_colname(name):
    '''Formats column name for plotting axis labels.

    Parameters
    ----------
    name : str
        Name of column.

    Returns
    -------
    fancyname : str
        Texified version of name.
    '''
    ns = name.split('_')
    if (len(ns) > 1) & (ns[-1] != 'indx'):
        if not name.startswith('HST'):
            newname = []
            if ns[0] in ['Av', 'Rv']:
                newname.append('${}_V$'.format(ns[0][0]))
                if len(ns) == 3:
                    newname.append(ns[1])
            elif name.startswith('logHST'):
                hstname = format_hstname(['HST'] + ns[1:-1])[1:]
                newname.append('$\log (\mathrm{' + '\ '.join(hstname) + '})$')
            elif ns[0].startswith('log'):
                newname.append('$\log ({})$'.format(ns[0][-1]))
            elif len(ns) == 3:
                newname.append('${}'.format(ns[0]) + '_{\mathrm{' + ns[1] + r'}}$')
            elif len(ns[0]) > 2:
                if ns[0] == 'radius':
                    newname.append('radius')
                else:
                    newname.append('${}'.format(ns[0][0]) + '_{\mathrm{' + ns[0][1:] + r'}}$')
            else:
                newname.append('${}$'.format(ns[0]))
            newname.append('({})'.format(ns[-1]))
        elif name.startswith('HST'):
            newname = format_hstname(ns)
    elif name[:4] == 'chi2':
        newname = ['$\chi^2$', name[4:]]
    else:
        newname = ns
    fancyname = ' '.join(newname)
    return r'{}'.format(fancyname)

def format_hstname(ns):
    ''' Formats HST filter names in style of 'HST/WFC3 F160W' or
    'HST/ACS-WFC F475W'

    Parameters
    ----------
    ns : list
        Column name split by underscore, beginning with HST

    Returns
    -------
    newname : list
        Formatted version of HST filter name plus any other items
    '''
    if 'ACS' in ns:
        newname = ['{}/{}-{}'.format(*ns[:3])] + ns[3:]
    else:
        newname = ['{}/{}'.format(*ns[:2])] + ns[2:]
    return newname

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
    '''Plot anything from the BEAST output against anything else as a 2D
    histogram or scatterplot.

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
