#!/usr/bin/env python
""" Make a nice plot of the filter response functions

"""
import argparse

import numpy as np

import matplotlib.pyplot as pyplot
import matplotlib 

from ..core import phot
from ..core.vega import Vega

def plot_filters_core(ax, filter_names, out_names, fontsize=15,
                      xrange=[0.2,2.0]):

    # wavelength grid for response functions
    waves_range = np.array([1e3,3e4])
    waves_range_log = np.log10(waves_range)
    n_waves = 501
    waves_delta_log = (waves_range_log[1] - waves_range_log[0])/(n_waves - 1)
    waves_log = np.arange(waves_range_log[0],waves_range_log[1],
                          waves_delta_log)
    waves = np.power(10.0,waves_log)

    # read in the filter response functions
    flist = phot.load_filters(filter_names, interp=True, lamb=waves)

    # change waves to microns for plotting
    waves /= 1e4

    for i, fname in enumerate(out_names):
        indxs, = np.where(flist[i].transmit > 0.0)
        ax.plot(waves[indxs],flist[i].transmit[indxs])
        lam_eff = np.trapz(flist[i].transmit[indxs]*waves[indxs],
                           waves[indxs])/np.trapz(flist[i].transmit[indxs],
                                                  waves[indxs])
        ax.text(lam_eff,2e0,fname,fontsize=1.5*fontsize,
                horizontalalignment='center')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xrange)
    ax.set_ylim(1e-8,2e1)
    ax.set_xlabel('$\lambda$ [$\mu m$]',fontsize=1.2*fontsize)
    ax.set_ylabel('$B_i(\lambda)$',fontsize=1.2*fontsize)

    ax.set_xticks([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

def plot_filters(args, filter_names, out_names,
                 save_name='beast_filters'):

    # make the nice plot
    fontsize = 15
    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    fig, ax = pyplot.subplots(nrows=1,ncols=1, sharex=True, figsize=(15,8))

    plot_filters_core(ax, filter_names, out_names, fontsize=fontsize)

    pyplot.tight_layout()

    if args.png:
        fig.savefig(save_name+'.png')
    elif args.eps:
        fig.savefig(save_name+'.eps')
    elif args.pdf:
        fig.savefig(save_name+'.pdf')
    else:
        pyplot.show()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    filter_names = ['HST_WFC3_F225W', 'HST_WFC3_F275W', 'HST_WFC3_F336W', 
                    'HST_ACS_WFC_F475W','HST_ACS_WFC_F550M',
                    'HST_ACS_WFC_F814W',
                    'HST_WFC3_F110W', 'HST_WFC3_F160W']
    out_names = ['F225W','F275W','F336W','F475W','F550M','F814W',
                 'F110W','F160W']

    plot_filters(args, filter_names, filters)

