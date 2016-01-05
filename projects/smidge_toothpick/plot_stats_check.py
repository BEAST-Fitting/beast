#!/usr/bin/env python
"""
Show diagnostic plots from BEAST runs
  Meant as a quick check that the results are reasonable

.. history::
    Written 21 Dec 2015 by Karl D. Gordon
"""

from __future__ import print_function

import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import matplotlib

from matplotlib.colors import LogNorm

from scipy.stats import binned_statistic_2d

from astropy.table import Table

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str,
                        help='FITS file with n band ASTs')
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    args = parser.parse_args()

    # pretty plotting details
    fontsize = 12
    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    # read in the results
    stats = Table.read(args.filename)

    print(np.sort(stats.colnames))

    # make the nice plot
    fig, ax = pyplot.subplots(nrows=2, ncols=3,
                              figsize=(15,10))

    # number of bins in x & y for 2d histrograms
    n_bins = 50
    dummy = 0

    # plot the HR diagram
    count_hr, xedges, yedges, binnum = \
              binned_statistic_2d(stats['logT_Exp'].quantity,
                                  stats['logL_Exp'].quantity,
                                  dummy,
                                  'count',
                                  bins=n_bins)
    cax1 = ax[0,0].imshow(count_hr.T, origin='lower',
                          norm=LogNorm(vmin=0.01),
                          extent=[xedges[0], xedges[-1], yedges[0],
                                  yedges[-1]],
                          aspect='auto', interpolation='nearest')

    ax[0,0].set_xlabel(r'$logT(Exp)$')
    ax[0,0].set_ylabel(r'$logL(Exp)$')

    ax[0,0].set_xlim(xedges[-1], xedges[0])

    # plot the chisqr values versus Teff
    indxs, = np.where(stats['chi2min'].quantity > 0.0)
    if len(indxs) > 0:
        count_chi2, xedges, yedges, binnum = \
                    binned_statistic_2d(stats['logT_Exp'].quantity[indxs],
                                        np.log10(stats['chi2min'].quantity[indxs]),
                                        dummy,
                                        'count',
                                        bins=n_bins)
        cax1 = ax[0,1].imshow(count_chi2.T, origin='lower',
                              norm=LogNorm(vmin=0.1),
                              extent=[xedges[0], xedges[-1], yedges[0],
                                      yedges[-1]],
                              aspect='auto', interpolation='nearest')

        ax[0,1].set_xlabel(r'$logT(Exp)$')
        ax[0,1].set_ylabel(r'$log(\chi^2)$')

        ax[0,1].set_xlim(xedges[-1], xedges[0])

    # plot the A(V) values versus Teff
    count_av, xedges, yedges, binnum = \
              binned_statistic_2d(stats['logT_Exp'],
                                  stats['Av_Exp'],
                                  stats['Av_Exp'],
                                  'count',
                                  bins=n_bins)
    cax1 = ax[1,0].imshow(count_av.T, origin='lower',
                          norm=LogNorm(vmin=0.1),
                          extent=[xedges[0], xedges[-1], yedges[0],
                                  yedges[-1]],
                          aspect='auto', interpolation='nearest')

    ax[1,0].set_xlabel(r'$logT(Exp)$')
    ax[1,0].set_ylabel(r'$A(V)$')

    ax[1,0].set_xlim(xedges[-1], xedges[0])

    # plot the A(V) values versus LogL
    count_av, xedges, yedges, binnum = \
              binned_statistic_2d(stats['logL_Exp'],
                                  stats['Av_Exp'],
                                  stats['Av_Exp'],
                                  'count',
                                  bins=n_bins)
    cax1 = ax[1,1].imshow(count_av.T, origin='lower',
                          norm=LogNorm(vmin=0.1),
                          extent=[xedges[0], xedges[-1], yedges[0],
                                  yedges[-1]],
                          aspect='auto', interpolation='nearest')

    ax[1,1].set_xlabel(r'$logL(Exp)$')
    ax[1,1].set_ylabel(r'$A(V)$')

    # plot the A(V) values versus R(V)
    count_av_rv, xedges, yedges, binnum = \
              binned_statistic_2d(stats['Av_Exp'],
                                  stats['Rv_Exp'],
                                  stats['Av_Exp'],
                                  'count',
                                  bins=n_bins)
    cax1 = ax[0,2].imshow(count_av_rv.T, origin='lower',
                          norm=LogNorm(vmin=0.1),
                          extent=[xedges[0], xedges[-1], yedges[0],
                                  yedges[-1]],
                          aspect='auto', interpolation='nearest')

    ax[0,2].set_xlabel(r'$A(V)$')
    ax[0,2].set_ylabel(r'$R(V)$')

    # plot the A(V) values versus f_A
    count_av_fa, xedges, yedges, binnum = \
              binned_statistic_2d(stats['Av_Exp'],
                                  stats['f_A_Exp'],
                                  stats['Av_Exp'],
                                  'count',
                                  bins=n_bins)
    cax1 = ax[1,2].imshow(count_av_fa.T, origin='lower',
                          norm=LogNorm(vmin=0.1),
                          extent=[xedges[0], xedges[-1], yedges[0],
                                  yedges[-1]],
                          aspect='auto', interpolation='nearest')

    ax[1,2].set_xlabel(r'$A(V)$')
    ax[1,2].set_ylabel(r'$f_A$')

    # colorbar 1
    #fig.colorbar(cax1, cax=(pyplot.subplot(gs[1:n_filters,0])))

    # colorbar 2
    #fig.colorbar(cax2, cax=(pyplot.subplot(gs[0:n_filters-1,n_filters+1])))

    # optimize the figure layout
    pyplot.tight_layout()

    # show or save
    basename = args.filename.replace('.fits','') + '_diagnostics'
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    else:
        pyplot.show()
    
