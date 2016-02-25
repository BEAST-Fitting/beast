#!/usr/bin/env python
"""
Code to plot the nxn band covariance matrix and biases for all
stars in AST file

.. history::
    Written 14 Dec 2015 by Karl D. Gordon
"""

from __future__ import print_function

import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import matplotlib
from matplotlib.ticker import MaxNLocator

from scipy.stats import binned_statistic_2d

from astropy.table import Table

from beast.core.noisemodel import trunchen
import datamodel_small as datamodel

class PHAT_Trunchen_Noisemodel(trunchen.MultiFilterASTs):
    pass

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str,
                        help='FITS file with n band ASTs')
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    # pretty plotting details
    fontsize = 14
    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    # define the filters
    filters = datamodel.basefilters
    n_filters = len(filters)

    # read in the AST file
    asts = PHAT_Trunchen_Noisemodel(args.filename, datamodel.filters)

    # calcuate all the covariance matrices for all the ASTs in the file
    all_covs, all_biases, all_compls, all_corrs, all_ifluxes, ast_minmax = \
              asts._calc_all_ast_cov(filters)

    # make the nice plot
    n_filters = len(filters)

    fig, ax = pyplot.subplots(nrows=n_filters, ncols=n_filters,
                              figsize=(14,14))
    gs = gridspec.GridSpec(n_filters+2, n_filters,
                           height_ratios=[0.25] + n_filters*[1.0] + [0.25])

    ax = []
    for i in range(n_filters):
        x = []
        for j in range(n_filters):
            x.append(pyplot.subplot(gs[i+1,j]))
        ax.append(x)
    ax = np.array(ax)

    # isolate odd stars
    #i = 4
    #indxs, = np.where((np.log10(all_ifluxes[:,i]) > -6.5) &
    #                  ((np.sqrt(all_covs[:,i,i])/all_ifluxes[:,i]) > 0.1))
    #print(indxs)

    # display images showing the corrlelation versus the fluxes in the
    # two bands
    n_bins = 20
    for i in range(n_filters):
        for j in range(n_filters-1,i,-1):
            mean_corr, xedges, yedges, binnum = \
                       binned_statistic_2d(np.log10(all_ifluxes[:,i]),
                                           np.log10(all_ifluxes[:,j]),
                                           all_corrs[:,j,i],'median',
                                           bins=n_bins)
            cax1 = ax[j,i].imshow(mean_corr, origin='lower',vmin=-1.0,vmax=1.0,
                       extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                       aspect='auto', interpolation='nearest')

            if j == n_filters-1:
                ax[j,i].set_xlabel('log[F('+filters[i]+')]')
            if i == 0:
                ax[j,i].set_ylabel('log[F('+filters[j]+')]')

            ax[j,i].xaxis.set_major_locator(MaxNLocator(4,integer=True))
            ax[j,i].yaxis.set_major_locator(MaxNLocator(4,integer=True))

            std_corr, xedges, yedges, binnum = \
                       binned_statistic_2d(np.log10(all_ifluxes[:,j]),
                                           np.log10(all_ifluxes[:,i]),
                                           all_corrs[:,i,j],
                                           statistic=np.std,
                                           bins=n_bins)
            cax2 = ax[i,j].imshow(std_corr, origin='lower',vmin=0.0,vmax=0.5,
                                  extent=[xedges[0], xedges[-1], yedges[0],
                                          yedges[-1]],
                                  aspect='auto', interpolation='nearest',
                                  cmap=pyplot.get_cmap('summer'))

            if i == n_filters-1:
                ax[i,j].set_xlabel('log[F('+filters[j]+')]')
            if j == 0:
                ax[i,j].set_ylabel('log[F('+filters[i]+')]')

            ax[i,j].xaxis.set_major_locator(MaxNLocator(4,integer=True))
            ax[i,j].yaxis.set_major_locator(MaxNLocator(4,integer=True))

    for i in range(n_filters):
        # plot the diagonal terms
        ax[i,i].plot(np.log10(all_ifluxes[:,i]),
                     np.log10(np.sqrt(all_covs[:,i,i])),'ro')
        #np.sqrt(all_covs[:,i,i])/all_ifluxes[:,i],'ro')
        #ax[i,i].set_yscale('log')
        #ax[i,i].set_ylabel(r'$log(\sigma)$')
        if i == n_filters-1:
            ax[i,i].set_xlabel(r'$log[F('+filters[i]+')]$')
        ax[i,i].xaxis.set_major_locator(MaxNLocator(4,integer=True))
        ax[i,i].yaxis.set_major_locator(MaxNLocator(4,integer=True))

        ax[i,i].text(0.1,0.8,r'log($\sigma$)',
                             fontsize=fontsize,
                             transform=ax[i,i].transAxes)

    # colorbar 1
    fig.colorbar(cax1, cax=(pyplot.subplot(gs[n_filters+1,0:n_filters-1])), 
                 orientation='horizontal')

    # colorbar 2
    fig.colorbar(cax2, cax=(pyplot.subplot(gs[0,1:n_filters])), 
                 orientation='horizontal')

    # optimize the figure layout
    pyplot.tight_layout(h_pad=0.1,w_pad=0.1)

    # show or save
    basename = args.filename.replace('.fits','') + '_all_ast_cov'
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    elif args.pdf:
        fig.savefig(basename+'.pdf')
    else:
        pyplot.show()
    
