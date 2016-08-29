#!/usr/bin/env python
"""
Script to plot CMDs

.. history::
    Written 25mar15 for the BEAST paper and general discussions
"""

from __future__ import print_function

import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec

from beast.core import grid
from beast.core.vega import Vega
from beast.core import extinction

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    args = parser.parse_args()

    # HR space to use
    #logL_range = [3.0, 3.5]
    #logT_range = [4.2, 4.4]
    logL_range = [1.65, 1.85]
    logT_range = [3.6, 3.8]
    
    # Dust extinction definition
    extLaw = extinction.RvFbumpLaw()

    # read in the BEAST generated model grid
    sedgrid_filename = 'smidge_dec15_small/smidge_dec15_small_seds.grid.hd5'
    sedgrid = grid.FileSEDGrid(sedgrid_filename, backend='cache')

    #print(sedgrid.keys())

    fluxes = sedgrid.seds
    filters = sedgrid.filters

    # get the zero mag flux for "vega" for the filters
    _, vega_flux, _ = (Vega().getFlux(filters))

    n_models, n_filters = fluxes.shape

    # setup the plot
    #fig, ax = pyplot.subplots(nrows=2, ncols=4, figsize=(20,10))
    fig, ax = pyplot.subplots(figsize=(15,10))

    # use gridspec to allow for one plot to be larger than the others
    gs = gridspec.GridSpec(3, 4, height_ratios=[1.0,1.0,2.0])
    ax = []
    for i in range(n_filters):
        gs_indxs = np.array(divmod(i,4))
        ax.append(pyplot.subplot(gs[gs_indxs[0],gs_indxs[1]]))
    # cut plots
    ax.append(pyplot.subplot(gs[2,0:4]))
    #ax.append(pyplot.subplot(gs[1,0:2]))
    #ax.append(pyplot.subplot(gs[2,0:2]))
    #ax.append(pyplot.subplot(gs[3,0:2]))

    # set the x axis filters
    xf1 = 3
    #xf2 = 5
    #xtitle = filters[xf1].split('_')[-1].upper() + ' - ' + filters[xf2].split('_')[-1].upper()
    xtitle = filters[xf1].split('_')[-1].upper()

    # plot the unreddened points
    mindxs, = np.where((sedgrid['Av'] == 0.0) &
                       (sedgrid['Rv'] == 3.0) &
                       (sedgrid['f_A'] == 1.0) &
                       (sedgrid['Z'] == sedgrid['Z'][0]) &
                       (sedgrid['logL'] > logL_range[0]) & (sedgrid['logL'] < logL_range[1]) &
                       (sedgrid['logT'] > logT_range[0]) & (sedgrid['logT'] < logT_range[1]))
                       #(sedgrid['logg'] > 2.2) & (sedgrid['logg'] < 2.6) &
                       #(sedgrid['logT'] > 3.65) & (sedgrid['logT'] < 3.7))
    mindxs_nr = mindxs

    #xvals = -2.5*(np.log10(fluxes[mindxs,xf1]/vega_flux[xf1]) - np.log10(fluxes[mindxs,xf2]/vega_flux[xf2]))
    xvals = -2.5*(np.log10(fluxes[mindxs,xf1]/vega_flux[xf1]))
    for i in range(n_filters):
        yvals = -2.5*np.log10(fluxes[mindxs,i]/vega_flux[i])

        ax[i].plot(xvals,yvals,'ko')
        ax[i].invert_yaxis()        
        ax[i].set_ylabel(filters[i].split('_')[-1].upper())

        dm = divmod(i,4)
        if dm[0] == 1:
            ax[i].set_xlabel(xtitle)

    # now plot a reddended version with "standard" MW dust
    #wrange = [1000.,2e5]
    #logwrange = np.log10(wrange)
    #mwaves = np.arange(logwrange[1],logwrange[0],0.1)
    #mwaves = np.power(10,mwaves)

    mwaves = np.logspace(np.log10(1.8e3),np.log10(2e4),100)

    plot_red = True
    if plot_red:
        mindxs, = np.where((sedgrid['Av'] == 1.0) &
                           (sedgrid['Rv'] == 2.0) &
                           (sedgrid['f_A'] == 1.0) &
                           (sedgrid['Z'] == sedgrid['Z'][0]) &
                           (sedgrid['logL'] > logL_range[0]) & (sedgrid['logL'] < logL_range[1]) &
                           (sedgrid['logT'] > logT_range[0]) & (sedgrid['logT'] < logT_range[1]))

        #xvals = -2.5*(np.log10(fluxes[mindxs,xf1]/vega_flux[xf1]) - np.log10(fluxes[mindxs,xf2]/vega_flux[xf2]))
        xvals = -2.5*(np.log10(fluxes[mindxs,xf1]/vega_flux[xf1]))
        xvals_nr = -2.5*(np.log10(fluxes[mindxs_nr,xf1]/vega_flux[xf1]))

        rv20ext = np.empty((n_filters))
        rv20ext_min = np.empty((n_filters))
        rv20ext_max = np.empty((n_filters))
        for i in range(n_filters):
            yvals = -2.5*np.log10(fluxes[mindxs,i]/vega_flux[i])
            yvals_nr = -2.5*np.log10(fluxes[mindxs_nr,i]/vega_flux[i])

            af1 = xvals_nr - xvals
            af2 = yvals_nr - yvals
            
            rv20ext[i] = np.mean(af2/af1)
            rv20ext_min[i] = np.min(af2/af1)
            rv20ext_max[i] = np.max(af2/af1)
            #print(rv20ext[i],rv20ext_min[i],rv20ext_max[i])
            
            dm = divmod(i,4)
            ax[i].plot(xvals,yvals,'ro')
            #ax[dm[0],dm[1]].plot(af1,af2,'ro')
            #ax[dm[0],dm[1]].set_ylabel(filters[i].split('_')[-1].upper())

        # get the A(lambda)/A(V) values for the nominal filter wavelengths
        alav_vals = extLaw.function(np.concatenate(([2000.],sedgrid.lamb)), Av=1., Rv=2.0, f_A=1.0, Alambda=True)
        alav_vals = alav_vals[1:]
        rv20ext_input = alav_vals/alav_vals[xf1]
        ax[n_filters].plot(sedgrid.lamb,rv20ext_input,'r-')

        talav_vals = extLaw.function(mwaves, Av=1., Rv=2.0, f_A=1.0, Alambda=True)
        rv20ext_input = talav_vals/alav_vals[xf1]
        ax[n_filters].plot(mwaves,rv20ext_input,'r--')

        ax[n_filters].errorbar(sedgrid.lamb,rv20ext,yerr=[rv20ext-rv20ext_min,rv20ext_max-rv20ext],fmt='ro',label=r'R(V)=2.0, $f_A$=1.0')
        ax[n_filters].set_xscale('log')
        ax[n_filters].set_xlim(2e3,2e4)
        ax[n_filters].set_xlabel(r'$\lambda$ [$\AA$]')
        ax[n_filters].set_ylabel(r'A($\lambda$)/A('+xtitle+')')

    plot_red = True
    if plot_red:
        mindxs, = np.where((sedgrid['Av'] == 1.0) &
                           (sedgrid['Rv'] == 6.0) &
                           (sedgrid['f_A'] == 1.0) &
                           (sedgrid['Z'] == sedgrid['Z'][0]) &
                           (sedgrid['logL'] > logL_range[0]) & (sedgrid['logL'] < logL_range[1]) &
                           (sedgrid['logT'] > logT_range[0]) & (sedgrid['logT'] < logT_range[1]))

        #xvals = -2.5*(np.log10(fluxes[mindxs,xf1]/vega_flux[xf1]) - np.log10(fluxes[mindxs,xf2]/vega_flux[xf2]))
        xvals = -2.5*(np.log10(fluxes[mindxs,xf1]/vega_flux[xf1]))
        xvals_nr = -2.5*(np.log10(fluxes[mindxs_nr,xf1]/vega_flux[xf1]))

        rv60ext = np.empty((n_filters))
        rv60ext_min = np.empty((n_filters))
        rv60ext_max = np.empty((n_filters))
        for i in range(n_filters):
            yvals = -2.5*np.log10(fluxes[mindxs,i]/vega_flux[i])
            yvals_nr = -2.5*np.log10(fluxes[mindxs_nr,i]/vega_flux[i])
            
            af1 = xvals_nr - xvals
            af2 = yvals_nr - yvals
            
            rv60ext[i] = np.mean(af2/af1)
            rv60ext_min[i] = np.min(af2/af1)
            rv60ext_max[i] = np.max(af2/af1)

            ax[i].plot(xvals,yvals,'go')

        # get the A(lambda)/A(V) values for the nominal filter wavelengths
        alav_vals = extLaw.function(np.concatenate(([2000.],sedgrid.lamb)), Av=1., Rv=6.0, f_A=1.0, Alambda=True)
        alav_vals = alav_vals[1:]
        rv60ext_input = alav_vals/alav_vals[xf1]
        ax[n_filters].plot(sedgrid.lamb,rv60ext_input,'g-')

        talav_vals = extLaw.function(mwaves, Av=1., Rv=6.0, f_A=1.0, Alambda=True)
        rv60ext_input = talav_vals/alav_vals[xf1]
        ax[n_filters].plot(mwaves,rv60ext_input,'g--')

        ax[n_filters].errorbar(sedgrid.lamb,rv60ext,yerr=[rv60ext-rv60ext_min,rv60ext_max-rv60ext],fmt='go',label=r'R(V)=6.0, $f_A$=1.0')

    plot_red = True
    if plot_red:
        mindxs, = np.where((sedgrid['Av'] == 1.0) &
                           (sedgrid['Rv'] == 3.0) &
                           (sedgrid['f_A'] == 1.0) &
                           (sedgrid['Z'] == sedgrid['Z'][0]) &
                           (sedgrid['logL'] > logL_range[0]) & (sedgrid['logL'] < logL_range[1]) &
                           (sedgrid['logT'] > logT_range[0]) & (sedgrid['logT'] < logT_range[1]))

        #xvals = -2.5*(np.log10(fluxes[mindxs,xf1]/vega_flux[xf1]) - np.log10(fluxes[mindxs,xf2]/vega_flux[xf2]))
        xvals = -2.5*(np.log10(fluxes[mindxs,xf1]/vega_flux[xf1]))
        xvals_nr = -2.5*(np.log10(fluxes[mindxs_nr,xf1]/vega_flux[xf1]))

        rv60ext = np.empty((n_filters))
        rv60ext_min = np.empty((n_filters))
        rv60ext_max = np.empty((n_filters))
        for i in range(n_filters):
            yvals = -2.5*np.log10(fluxes[mindxs,i]/vega_flux[i])
            yvals_nr = -2.5*np.log10(fluxes[mindxs_nr,i]/vega_flux[i])
            
            af1 = xvals_nr - xvals
            af2 = yvals_nr - yvals
            
            rv60ext[i] = np.mean(af2/af1)
            rv60ext_min[i] = np.min(af2/af1)
            rv60ext_max[i] = np.max(af2/af1)

            ax[i].plot(xvals,yvals,'co')

        # get the A(lambda)/A(V) values for the nominal filter wavelengths
        alav_vals = extLaw.function(np.concatenate(([2000.],sedgrid.lamb)), Av=1., Rv=3.0, f_A=1.0, Alambda=True)
        alav_vals = alav_vals[1:]
        rv60ext_input = alav_vals/alav_vals[xf1]
        ax[n_filters].plot(sedgrid.lamb,rv60ext_input,'c-')

        talav_vals = extLaw.function(mwaves, Av=1., Rv=3.0, f_A=0.2, Alambda=True)
        rv60ext_input = talav_vals/alav_vals[xf1]
        ax[n_filters].plot(mwaves,rv60ext_input,'c--')

        ax[n_filters].errorbar(sedgrid.lamb,rv60ext,yerr=[rv60ext-rv60ext_min,rv60ext_max-rv60ext],fmt='co',label=r'R(V)=3.0, $f_A$=1.0')

    plot_red = True
    if plot_red:
        mindxs, = np.where((sedgrid['Av'] == 1.0) &
                           (sedgrid['Rv'] == 3.0) &
                           (sedgrid['f_A'] == 0.2) &
                           (sedgrid['Z'] == sedgrid['Z'][0]) &
                           (sedgrid['logL'] > logL_range[0]) & (sedgrid['logL'] < logL_range[1]) &
                           (sedgrid['logT'] > logT_range[0]) & (sedgrid['logT'] < logT_range[1]))

        #xvals = -2.5*(np.log10(fluxes[mindxs,xf1]/vega_flux[xf1]) - np.log10(fluxes[mindxs,xf2]/vega_flux[xf2]))
        xvals = -2.5*(np.log10(fluxes[mindxs,xf1]/vega_flux[xf1]))
        xvals_nr = -2.5*(np.log10(fluxes[mindxs_nr,xf1]/vega_flux[xf1]))

        rv60ext = np.empty((n_filters))
        rv60ext_min = np.empty((n_filters))
        rv60ext_max = np.empty((n_filters))
        for i in range(n_filters):
            yvals = -2.5*np.log10(fluxes[mindxs,i]/vega_flux[i])
            yvals_nr = -2.5*np.log10(fluxes[mindxs_nr,i]/vega_flux[i])
            
            af1 = xvals_nr - xvals
            af2 = yvals_nr - yvals
            
            rv60ext[i] = np.mean(af2/af1)
            rv60ext_min[i] = np.min(af2/af1)
            rv60ext_max[i] = np.max(af2/af1)

            ax[i].plot(xvals,yvals,'bo')

        # get the A(lambda)/A(V) values for the nominal filter wavelengths
        alav_vals = extLaw.function(np.concatenate(([2000.],sedgrid.lamb)), Av=1., Rv=3.0, f_A=0.2, Alambda=True)
        alav_vals = alav_vals[1:]
        rv60ext_input = alav_vals/alav_vals[xf1]
        ax[n_filters].plot(sedgrid.lamb,rv60ext_input,'b-')

        talav_vals = extLaw.function(mwaves, Av=1., Rv=3.0, f_A=0.2, Alambda=True)
        rv60ext_input = talav_vals/alav_vals[xf1]
        ax[n_filters].plot(mwaves,rv60ext_input,'b--')

        ax[n_filters].errorbar(sedgrid.lamb,rv60ext,yerr=[rv60ext-rv60ext_min,rv60ext_max-rv60ext],fmt='bo',label=r'R(V)=3.0, $f_A$=0.2')

    ax[n_filters].legend()

    ax[n_filters].text(5000.,3.0,'%.1f' % logL_range[0] + ' < LogL < ' + '%.1f' % logL_range[1])
    ax[n_filters].text(5000.,2.8,'%.1f' % logT_range[0] + ' < LogT < ' + '%.1f' % logT_range[1])
    
    # optimize the figure layout
    gs.tight_layout(fig, rect=[0, 0.03, 1, 0.96])

    if args.png:
        fig.savefig('smidge_red_clump_ext_sim.png')
    elif args.eps:
        fig.savefig('smidge_red_clump_ext_sim.eps')
    else:
        pyplot.show()
