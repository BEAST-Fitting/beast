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
    args = parser.parse_args()

    # pretty plotting details
    fontsize = 12
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
                              figsize=(20,14))
    gs = gridspec.GridSpec(n_filters, n_filters+2,
                           width_ratios=[0.25] + n_filters*[1.0] + [0.25])

    ax = []
    for i in range(n_filters):
        x = []
        for j in range(n_filters):
            x.append(pyplot.subplot(gs[i,j+1]))
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

            ax[j,i].set_xlabel(r'$log[F('+filters[i]+')]$')
            ax[j,i].set_ylabel(r'$log[F('+filters[j]+')]$')

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

            ax[i,j].set_xlabel(r'$log[F('+filters[j]+')]$')
            ax[i,j].set_ylabel(r'$log[F('+filters[i]+')]$')

    for i in range(n_filters):
        # plot the diagonal terms
        ax[i,i].plot(np.log10(all_ifluxes[:,i]),
                     np.sqrt(all_covs[:,i,i])/all_ifluxes[:,i],'ro')
        ax[i,i].set_yscale('log')
        ax[i,i].set_ylabel(r'$\sigma/F('+filters[i]+')$')
        ax[i,i].set_xlabel(r'$log[F('+filters[i]+')]$')

    # colorbar 1
    fig.colorbar(cax1, cax=(pyplot.subplot(gs[1:n_filters,0])))

    # colorbar 2
    fig.colorbar(cax2, cax=(pyplot.subplot(gs[0:n_filters-1,n_filters+1])))

    # optimize the figure layout
    pyplot.tight_layout(h_pad=0.1,w_pad=0.1)

    # show or save
    basename = args.filename.replace('.fits','') + '_all_ast_cov'
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    else:
        pyplot.show()
    
