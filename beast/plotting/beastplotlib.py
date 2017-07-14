#!/usr/bin/env python
"""
Library of general functions for the BEAST plotting scripts

.. history::
    Written 21 Feb 2017 by Meredith J. Durbin
"""
from __future__ import (absolute_import, division, print_function)

import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

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
    parser.add_argument('-s', '--savefig', action='store',
                        default=False, choices=ftypes,
                        help='Save figure to a file of specified type \
                            rather than show it. \
                            Must be one of: "{}"'.format('", "'.join(ftypes))
                        )
    parser.add_argument('-t', '--tex', action='store_true', 
                        help='Configure Matplotlib to use LaTeX; \
                              sets rcParam "usetex" to True. \
                              (Requires working TeX installation.)\
                              Defaults to false.')
    return parser

def plot_generic(table, xcolname, ycolname, fig=None, ax=None, plottype='hist',
                 bins=51, thresh_col=None, thresh=0, thresh_op='less',
                 xlog=False, ylog=False,
                 plot_kwargs={}):
    '''Plot anything from the BEAST output against anything else as a 2D
    histogram or scatterplot.

    Parameters
    ----------
    table : astropy table or pandas dataframe
        Table of BEAST output
    xcolname : str
        Name of table column with x values
    ycolname : str
        Name of table column with x values
    fig : matplotlib figure object, optional
        Figure to plot histogram on
    ax : matplotlib axes object, optional
        Subplot of figure to plot on
    plottype : str, optional
        String indicating what type of plot to make.
        Must be one of: {'hist', 'scatter'}. Defaults to 'hist'
    bins : int or array-like, optional
        Number of bins to pass to 2D histogram.
        Only used if plottype == 'hist'.
    thresh_col : str, optional
        Name of column with values to filter xcol and ycol by; for example, 
        if you want to look at only stars with chi^2 greater than 20 
        you would set thresh_col = 'chi2min', thresh=20, and 
        thresh_op = 'greater'
    thresh : scalar, optional
        Threshold of values to compare thresh_col to
    thresh_op : str, optional
        Operator comparison to perform between thresh_col and thresh.
        Must be one of Numpy logical operator functions: {'greater', 
        'greater_equal', 'less', 'less_equal', 'equal', 'not_equal'}
    xlog : boolean, optional
        plot the x axis in log units
    ylog : boolean, optional
        plot the y axis in log units
    plot_kwargs : dict, optional
        Additional keyword arguments to pass to plotting function.

    Returns
    -------
    fig : matplotlib figure object
        Figure
    ax : matplotlib subplot object
        Subplot
    '''
    return_vals = []
    if fig is None:
        fig = plt.gcf()
        return_vals.append(fig)
    if ax is None:
        ax = plt.gca()
        return_vals.append(ax)
    if thresh_col is not None:
        condition = getattr(np, thresh_op)(table[thresh_col], thresh)
        table = table[condition]
    xcol = table[xcolname]
    ycol = table[ycolname]

    if len(xcol) <= 0:
        return
    
    if xlog:
        xcol = np.log10(xcol)
        xcolname = 'log10(' + xcolname + ')'
    if ylog:
        ycol = np.log10(ycol)
        ycolname = 'log10(' + ycolname + ')'
    if plottype == 'hist':
        hist = ax.hist2d(xcol, ycol, bins=bins, **plot_kwargs)
        cbar = fig.colorbar(hist[-1], label='Number of stars', ax=ax)
        return_vals.append(cbar)
    elif plottype == 'scatter':
        if ('c' in plot_kwargs):
            if (thresh_col is not None):
                if len(plot_kwargs['c']) == len(condition):
                    plot_kwargs['c'] = plot_kwargs['c'][condition]
        sc = ax.scatter(xcol, ycol, **plot_kwargs)
        if ('c' in plot_kwargs):
            if (len(plot_kwargs['c']) == len(xcol)):
                cbar = fig.colorbar(sc, label=fancify_colname(plot_kwargs['c'].name), ax=ax)
                return_vals.append(cbar)
    ax.set_xlabel(fancify_colname(xcolname))
    ax.set_ylabel(fancify_colname(ycolname))
    xmin, xmax = ax.get_xlim()
    if (xmin < xmax) & xcolname.startswith('logT'):
        ax.invert_xaxis()
    ymin, ymax = ax.get_ylim()
    if (ymin < ymax) & (ycolname.startswith('logT') | ycolname.startswith('mbol')):
        ax.invert_yaxis()
    if len(return_vals) > 0:
        return return_vals

def set_params(lw=1.5, universal_color='#262626', fontsize=12, helvetica=True, usetex=True):
    '''Configure some matplotlib rcParams.

    Parameters
    ----------
    lw : scalar
        Linewidth of plot and axis lines. Default is 1.5.
    universal_color : str, matplotlib named color, rgb tuple
        Color of text and axis spines. Default is #262626, off-black
    fontsize : scalar
        Font size in points. Default is 12
    helvetica : boolean
        Whether to try and change font to Helvetica. Default is true.
    usetex : boolean
        Whether to set rc.text.usetex = True. Default is true.

    Returns
    -------
    fig : matplotlib figure object
        Figure
    ax : matplotlib subplot object
        Subplot
    '''
    rc('font', size=fontsize)
    if usetex:
        try:
            rc('text', usetex=True)
        except:
            print("Couldn't configure usetex. Make sure you have a working TeX installation.")
    if helvetica:
        try:
            rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
            preamble = [r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
                        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
                        r'\usepackage{helvet}',    # set the normal font here
                        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
                        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
                       ]
            rc('text.latex', preamble = preamble)
        except:
            print("Couldn't change font to Helvetica. Falling back to Matplotlib defaults.")
    rc('lines', linewidth=lw, markeredgewidth=lw*0.5)
    rc('patch', linewidth=lw, edgecolor='#FAFAFA')
    rc('axes', linewidth=lw, edgecolor=universal_color, labelcolor=universal_color, 
        axisbelow=True)
    rc('image', origin='lower') # fits images
    rc('xtick.major', width=lw*0.75)
    rc('xtick.minor', width=lw*0.5)
    rc('xtick', color=universal_color)
    rc('ytick.major', width=lw*0.75)
    rc('ytick.minor', width=lw*0.5)
    rc('ytick', color=universal_color)
    rc('grid', linewidth=lw)
    rc('legend', loc='best', numpoints=1, scatterpoints=1, handlelength=1.5,
        fontsize=fontsize, columnspacing=1, handletextpad=0.75)
    
