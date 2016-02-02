#!/usr/bin/env python
"""
Code to plot nxn band covariance matrix and biases from n band ASTs

.. history::
    Written 7 Dec 2015 by Karl D. Gordon
"""

from __future__ import print_function

import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

from beast.core.noisemodel import trunchen
import datamodel_small as datamodel

class PHAT_Trunchen_Noisemodel(trunchen.MultiFilterASTs):
    pass

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str,
                        help='FITS file with n band ASTs')
    parser.add_argument("-s", "--starnum", type=int, default=0,
                        help="star number in unique list")
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
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

    # read in the AST file
    asts = PHAT_Trunchen_Noisemodel(args.filename, datamodel.filters)

    # find the stars by using unique values of the magnitude vlaues in
    #    filtername
    filters = datamodel.basefilters
    filtername = filters[5] + '_IN'
    uvals, ucounts = np.unique(asts.data[filtername], return_counts=True)
    print('# of models = ', len(uvals))

    # find all the ASTs for this model
    indxs, = np.where(asts.data[filtername] == uvals[args.starnum])
    n_indxs = len(indxs)

    # calculate the covariance matrix and biases for a single model
    cov_matrix, biases, stddevs, corr_matrix, \
                diffs, ifluxes, compl = asts._calc_ast_cov(indxs, filters,
                                                           return_all=True)

    # make the nice plot
    n_filters = len(filters)
    fig, ax = pyplot.subplots(nrows=n_filters, ncols=n_filters,
                              figsize=(10,10))

    # plot the differences normalized by the standard deviations
    for i in range(n_filters):
        for j in range(n_filters-1,i,-1):
            # plot the AST data
            ax[j,i].plot(diffs[i,:]/stddevs[i], diffs[j,:]/stddevs[j], 'o',
                         markerfacecolor='Chartreuse', markersize=4.)
            # plot the biases
            ax[j,i].plot(biases[i]/stddevs[i], biases[j]/stddevs[j], 'bo',
                         markersize=10.)
            # plot "cross hairs"
            ylim = np.array(ax[j,i].get_ylim())
            ylim[0] = min([ylim[0],-3.0])
            ylim[1] = max([ylim[1],3.0])
            ax[j,i].plot([0.0,0.0], ylim, 'k--')
            xlim = np.array(ax[j,i].get_xlim())
            xlim[0] = min([xlim[0],-3.0])
            xlim[1] = max([xlim[1],3.0])
            ax[j,i].plot(xlim, [0.0,0.0], 'k--')

            # plot an ellipse that illustrates the covariance
            theta = 2.*np.pi*np.linspace(0.0,1.0,num=100)
            rot_angle = np.arctan(corr_matrix[i,j])
            if corr_matrix[i,j] < 0:
                rot_angle = -1.0*rot_angle
            a = 1.0/np.cos(rot_angle)
            b = a*(1.0-corr_matrix[i,j])
            r = a*b/np.sqrt(np.square(b*np.cos(theta)) +
                            np.square(a*np.sin(theta)))
            x = r*np.cos(theta)
            y = r*np.sin(theta)

            new_x = x*np.cos(rot_angle) - y*np.sin(rot_angle)
            new_y = x*np.sin(rot_angle) + y*np.cos(rot_angle)

            new_x += biases[i]/stddevs[i]
            new_y += biases[j]/stddevs[j]

            ax[j,i].plot(new_x,new_y,'k-')

    # label the yaxis
    for i in range(1,n_filters):
        ax[i,0].set_ylabel(r'$\mu/\sigma$')

    # label the xaxis
    for i in range(0,n_filters-1):
        ax[n_filters-1,i].set_xlabel(r'$\mu/\sigma$')

    # suppress labeling of x-axis ticks
    for i in range(0,n_filters-1):
        for j in range(n_filters):
            ax[i,j].xaxis.set_ticklabels([])
    
    # suppress labeling of y-axis ticks
    for i in range(1,n_filters):
        for j in range(n_filters):
            ax[j,i].yaxis.set_ticklabels([])

    # suppress axes for the "text" plots
    # add the correlation value
    for i in range(n_filters):
        for j in range(i,n_filters):
            ax[i,j].get_xaxis().set_visible(False)
            ax[i,j].get_yaxis().set_visible(False)

            if i != j:
                ax[i,j].text(0.5,0.5,'% 5.2f' % corr_matrix[i,j],
                             horizontalalignment='center',
                             verticalalignment='center',
                             fontsize=1.8*fontsize,
                             transform=ax[i,j].transAxes)
                redcol = 1.0
                bluecol = 1.0
                if corr_matrix[i,j] > 0:
                    bluecol = 1.0 - corr_matrix[i,j]
                    greencol = bluecol
                else:
                    redcol = corr_matrix[i,j] + 1.0
                    greencol = redcol
                ax[i,j].set_axis_bgcolor((redcol, greencol, bluecol))
            else:
                ax[i,j].text(0.6,0.6,filters[i],rotation=-45.0,
                             horizontalalignment='center',
                             verticalalignment='center',
                             fontsize=1.8*fontsize,
                             transform=ax[i,j].transAxes)
                imag = -2.5*np.log10(ifluxes[i]/asts.vega_flux[i])
                ax[i,j].text(0.35,0.35,'% 6.2f' % imag + ' mag',
                             rotation=-45.0,
                             horizontalalignment='center',
                             verticalalignment='center',
                             fontsize=1.1*fontsize,
                             transform=ax[i,j].transAxes)
                sn = ifluxes[i]/stddevs[i]
                ax[i,j].text(0.46,0.46,'S/N = ' + '% 6.2f' % sn,
                             rotation=-45.0,
                             horizontalalignment='center',
                             verticalalignment='center',
                             fontsize=1.1*fontsize,
                             transform=ax[i,j].transAxes)
                
            

    # optimize the figure layout
    pyplot.tight_layout(h_pad=0.1,w_pad=0.1)

    # show or save
    basename = 'plot_ast_cov'
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    elif args.pdf:
        fig.savefig(basename+'.pdf')
    else:
        pyplot.show()
    
