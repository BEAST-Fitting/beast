#!/usr/bin/env python
"""
Plot the individual fit for a single observed star

.. history::
    Written 12 Jan 2016 by Karl D. Gordon
      based on code written by Heddy Arab for the BEAST techniques paper figure
"""

from __future__ import print_function

import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
import matplotlib

from astropy.table import Table
from astropy.io import fits

from astropy import units as ap_units
from astropy.coordinates import SkyCoord as ap_SkyCoord

def disp_str(stars, k, keyname):

    dvals = [stats[keyname+'_p50'][k],
             stats[keyname+'_p84'][k],
             stats[keyname+'_p16'][k]]
    if keyname == 'M_ini':
        dvals = np.log10(dvals)
    disp_str = '$' + \
               "{0:.2f}".format(dvals[0]) + \
               '^{+' + \
               "{0:.2f}".format(dvals[1] - dvals[0]) + \
               '}_{-' + \
               "{0:.2f}".format(dvals[0] - dvals[2]) + \
               '}$'

    return disp_str


def plot_1dpdf(ax, pdf1d_hdu, tagname, xlabel, starnum,
               stats=None,
               logx=False, fontsize=None):

    pdf = pdf1d_hdu[tagname].data

    n_objects, n_bins = pdf.shape
    n_objects -= 1

    xvals = pdf[n_objects,:]
    if logx:
        xvals = np.log10(xvals)

    if tagname == 'Z':
        gindxs, = np.where(pdf[starnum,:] > 0.)
        ax.plot(xvals[gindxs],pdf[starnum,gindxs]/max(pdf[starnum,gindxs]),
                color='k')
    else:
        ax.plot(xvals,pdf[starnum,:]/max(pdf[starnum,:]),color='k')

    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim = [xvals.min(), xvals.max()]
    xlim_delta = xlim[1] - xlim[0]
    ax.set_xlim(xlim[0]-0.05*xlim_delta, xlim[1]+0.05*xlim_delta)
    #ax.set_ylim(0.0,1.1*pdf[starnum,:].max())
    ax.set_ylim(0.0,1.1)

    ax.text(0.95, 0.95, xlabel, transform=ax.transAxes,
            verticalalignment='top',horizontalalignment='right',
            fontsize=fontsize)

    if stats is not None:
        ylim = ax.get_ylim()

        y1 = ylim[0] + 0.5*(ylim[1]-ylim[0])
        y2 = ylim[0] + 0.7*(ylim[1]-ylim[0])
        pval = stats[tagname+'_Best'][starnum]
        if logx:
            pval = np.log10(pval)
        ax.plot(np.full((2),pval),[y1,y2],
                '-', color='c')

        y1 = ylim[0] + 0.2*(ylim[1]-ylim[0])
        y2 = ylim[0] + 0.4*(ylim[1]-ylim[0])
        y1m = ylim[0] + 0.25*(ylim[1]-ylim[0])
        y2m = ylim[0] + 0.35*(ylim[1]-ylim[0])
        ym = 0.5*(y1 + y2)
        pvals = [stats[tagname+'_p50'][starnum],
                 stats[tagname+'_p16'][starnum],
                 stats[tagname+'_p84'][starnum]]
        if logx:
            pvals = np.log10(pvals)
        ax.plot(np.full((2),pvals[0]),[y1m,y2m],'-', color='m')
        ax.plot(np.full((2),pvals[1]),[y1,y2],'-', color='m')
        ax.plot(np.full((2),pvals[2]),[y1,y2],'-', color='m')
        ax.plot(pvals[1:3],[ym,ym],'-', color='m')

def plot_beast_ifit(filters, waves, stats, pdf1d_hdu, fontsize,
                    surveyname='PHAT'):

    # setup the plot grid
    gs = gridspec.GridSpec(4, 4,
                           height_ratios=[1.0,1.0,1.0,1.0],
                           width_ratios=[1.0,1.0,1.0,1.0])
    ax = []
    # plots for the 1D PDFs
    for j in range(2):
        for i in range(4):
            ax.append(pyplot.subplot(gs[j+2,i]))
    # now for the big SED plot
    ax.append(pyplot.subplot(gs[0:2,0:3]))

    # plot the SED
    #print(np.sort(stats.colnames))

    n_filters = len(filters)
    
    # get the observations
    waves *= 1e-4
    obs_flux = np.empty((n_filters),dtype=np.float)
    mod_flux = np.empty((n_filters,3),dtype=np.float)
    mod_flux_nd = np.empty((n_filters,3),dtype=np.float)
    mod_flux_wbias = np.empty((n_filters,3),dtype=np.float)
    k = starnum

    c = ap_SkyCoord(ra=stats['RA'][k]*ap_units.degree,
                    dec=stats['DEC'][k]*ap_units.degree,
                    frame='icrs')
    corname = (surveyname + ' J' + 
               c.ra.to_string(unit=ap_units.hourangle, sep="",precision=2,
                              alwayssign=False,pad=True) + 
               c.dec.to_string(sep="",precision=2,
                               alwayssign=True,pad=True))
        
    for i, cfilter in enumerate(filters):
        obs_flux[i] = stats[cfilter][k]
        mod_flux[i,0] = np.power(10.0,stats['log'+cfilter+'_wd_p50'][k])
        mod_flux[i,1] = np.power(10.0,stats['log'+cfilter+'_wd_p16'][k])
        mod_flux[i,2] = np.power(10.0,stats['log'+cfilter+'_wd_p84'][k])
        mod_flux_nd[i,0] = np.power(10.0,stats['log'+cfilter+'_nd_p50'][k])
        mod_flux_nd[i,1] = np.power(10.0,stats['log'+cfilter+'_nd_p16'][k])
        mod_flux_nd[i,2] = np.power(10.0,stats['log'+cfilter+'_nd_p84'][k])
        if 'log'+cfilter+'_wd_bias_p50' in stats.colnames:
            mod_flux_wbias[i,0] = np.power(10.0,stats['log'+cfilter+
                                                      '_wd_bias_p50'][k])
            mod_flux_wbias[i,1] = np.power(10.0,stats['log'+cfilter+
                                                      '_wd_bias_p16'][k])
            mod_flux_wbias[i,2] = np.power(10.0,stats['log'+cfilter+
                                                      '_wd_bias_p84'][k])
        
    ax[8].plot(waves, obs_flux, 'ko', label='observed')

    if 'log'+filters[0]+'_wd_bias_p50' in stats.colnames:
        ax[8].plot(waves, mod_flux_wbias[:,0], 'b-',label='stellar+dust+bias')
        ax[8].fill_between(waves, mod_flux_wbias[:,1], mod_flux_wbias[:,2],
                           color='b', alpha = 0.3)

    ax[8].plot(waves, mod_flux[:,0], 'r-',label='stellar+dust')
    ax[8].fill_between(waves, mod_flux[:,1], mod_flux[:,2],
                       color='r', alpha = 0.2)

    ax[8].plot(waves, mod_flux_nd[:,0], 'y-',label='stellar only')
    ax[8].fill_between(waves, mod_flux_nd[:,1], mod_flux_nd[:,2],
                       color='y', alpha = 0.1)

    ax[8].legend(loc='upper right', bbox_to_anchor=(1.25, 1.025))

    ax[8].set_ylabel(r'Flux [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]',
                     fontsize=1.25*fontsize)
    ax[8].set_yscale('log')

    ax[8].set_xscale('log')
    ax[8].text(0.5,-0.01,r'$\lambda$ [$\AA$]', fontsize=1.25*fontsize,
               transform=ax[8].transAxes, verticalalignment='top')
    ax[8].set_xlim(0.2,2.0)
    ax[8].set_xticks([0.2,0.3,0.4,0.5,0.8,0.9,1.0,2.0])
    ax[8].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    ax[8].text(0.05, 0.95, corname, transform=ax[8].transAxes,
               verticalalignment='top',horizontalalignment='left',
               fontsize=1.25*fontsize)

    # add the text results
    keys = ['Av','M_ini','logA','Rv','f_A','Z','logT','logg','logL']
    dispnames = ['A(V)','log(M)','log(t)','R(V)',r'f$_\mathcal{A}$','Z',
                 r'log(T$_\mathrm{eff})$','log(g)','log(L)']
    laby = 0.72
    ty = np.linspace(laby-0.07,0.1,num=len(keys))
    ty[3:] -= 0.025
    ty[6:] -= 0.025
    tx = [1.12, 1.2, 1.3]
    for i in range(len(keys)):
        ax[8].text(tx[0], ty[i], dispnames[i],
                   horizontalalignment='right',
                   transform=ax[8].transAxes, fontsize=1.25*fontsize)
        ax[8].text(tx[1], ty[i], disp_str(stats, starnum, keys[i]),
                   horizontalalignment='center', color='m',
                   transform=ax[8].transAxes, fontsize=1.25*fontsize)
        best_val = stats[keys[i]+'_Best'][k]
        if keys[i] == 'M_ini':
            best_val = np.log10(best_val)
        ax[8].text(tx[2], ty[i],
                   '$' + "{0:.2f}".format(best_val) + '$', 
                   horizontalalignment='center', color='c', 
                   transform=ax[8].transAxes, fontsize=1.25*fontsize)
    ax[8].text(tx[0],laby, 'Param',
               horizontalalignment='right',
               transform=ax[8].transAxes, fontsize=1.25*fontsize)
    ax[8].text(tx[1],laby, '50%$\pm$33%',
               horizontalalignment='center', color='k',
               transform=ax[8].transAxes, fontsize=fontsize)
    ax[8].text(tx[2],laby, 'Best',color='k',
               horizontalalignment='center',
               transform=ax[8].transAxes, fontsize=1.25*fontsize)

    # now draw boxes around the different kinds of parameters
    tax = ax[8]

    # primary
    rec = Rectangle((tx[0]-0.1,ty[2]-0.02),
                    tx[2]-tx[0]+0.15, (ty[0]-ty[2])*1.5,
                    fill=False, lw=2, transform=tax.transAxes,
                    ls='dashed')
    rec = tax.add_patch(rec)
    rec.set_clip_on(False)

    # secondary
    rec = Rectangle((tx[0]-0.1,ty[5]-0.02),
                    tx[2]-tx[0]+0.15, (ty[3]-ty[5])*1.5,
                    fill=False, lw=2, transform=tax.transAxes,
                    ls='dotted')
    rec = tax.add_patch(rec)
    rec.set_clip_on(False)

    # derived
    rec = Rectangle((tx[0]-0.1,ty[8]-0.02),
                    tx[2]-tx[0]+0.15, (ty[6]-ty[8])*1.5,
                    fill=False, lw=2, transform=tax.transAxes,
                    ls='dashdot')
    rec = tax.add_patch(rec)
    rec.set_clip_on(False)

    # padding for rectangles of 1D PDFs
    pad = 0.1

    # plot the primary parameter 1D PDFs
    plot_1dpdf(ax[0], pdf1d_hdu, 'Av', 'A(V)', starnum,
               stats=stats, fontsize=fontsize)
    plot_1dpdf(ax[1], pdf1d_hdu, 'M_ini', 'log(M)', starnum, logx=True,
               stats=stats, fontsize=fontsize)
    plot_1dpdf(ax[2], pdf1d_hdu, 'logA', 'log(t)', starnum,
               stats=stats, fontsize=fontsize)

    # draw a box around them and label
    tax = ax[0]
    rec = Rectangle((-1.75*pad,-pad),
                    3*(1.0+pad)+1.5*pad,
                    1.0+1.5*pad,
                    fill=False, lw=2, transform=tax.transAxes,
                    ls='dashed')
    rec = tax.add_patch(rec)
    rec.set_clip_on(False)

    tax.text(-2.*pad, 0.5, 'Primary', transform=tax.transAxes, 
             rotation='vertical', fontstyle='oblique',
             verticalalignment='center', horizontalalignment='right',
             fontsize=1.5*fontsize)

    tax.text(0.0, 0.5, 'Probability', transform=tax.transAxes, 
             rotation='vertical', 
             verticalalignment='center', horizontalalignment='right',
             fontsize=1.25*fontsize)

    # plot the secondary parameter 1D PDFs
    plot_1dpdf(ax[4], pdf1d_hdu, 'Rv', 'R(V)', starnum,
               stats=stats, fontsize=fontsize)
    plot_1dpdf(ax[5], pdf1d_hdu, 'f_A', r'f$_\mathcal{A}$', starnum,
               stats=stats, fontsize=fontsize)
    plot_1dpdf(ax[6], pdf1d_hdu, 'Z', 'Z', starnum,
               stats=stats, fontsize=fontsize)

    # draw a box around them
    tax = ax[4]
    rec = Rectangle((-1.75*pad,-pad),
                    3*(1.0+pad)+1.5*pad,
                    1.0+1.5*pad,
                    fill=False, lw=2, transform=tax.transAxes,
                    ls='dotted')
    rec = tax.add_patch(rec)
    rec.set_clip_on(False)

    tax.text(-2*pad, 0.5, 'Secondary', transform=tax.transAxes, 
             rotation='vertical', fontstyle='oblique',
             verticalalignment='center', horizontalalignment='right',
             fontsize=1.5*fontsize)

    tax.text(0.0, 0.5, 'Probability', transform=tax.transAxes, 
             rotation='vertical', 
             verticalalignment='center', horizontalalignment='right',
             fontsize=1.25*fontsize)

    # plot the derived parameter 1D PDFs
    plot_1dpdf(ax[3], pdf1d_hdu, 'logT', r'log(T$_\mathrm{eff})$', starnum,
               stats=stats, fontsize=fontsize)
    plot_1dpdf(ax[7], pdf1d_hdu, 'logg', 'log(g)', starnum,
               stats=stats, fontsize=fontsize)

    # draw a box around them
    tax = ax[7]
    rec = Rectangle((-0.25*pad,-pad),
                    1.0+0.5*pad,
                    2*(1.0+2.*pad)-0.125*pad,
                    fill=False, lw=2, transform=tax.transAxes,
                    ls='dashdot')
    rec = tax.add_patch(rec)
    rec.set_clip_on(False)

    tax.text(0.5, 2*(1.0+2*pad)-1.0*pad, 'Derived', transform=tax.transAxes, 
             rotation='horizontal', fontstyle='oblique',
             verticalalignment='bottom', horizontalalignment='center',
             fontsize=1.5*fontsize)

    # optimize the figure layout
    pyplot.tight_layout(h_pad=2.0, w_pad=1.0)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("filebase", type=str,
                        help='base filenamae of output results')
    parser.add_argument("-s", "--starnum", type=int, default=0,
                        help="star number in observed file")
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    parser.add_argument("--surveyname", type=str, default='PHAT',
                        help="name of survey")
    args = parser.parse_args()

    starnum = args.starnum
    surveyname = args.surveyname

    # pretty plotting details
    fontsize = 12
    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    # base filename
    filebase = args.filebase

    # read in the stats
    stats = Table.read(filebase + '_stats.fits')
    #print(stats.colnames)
    #indxs, = np.where((stats['logT_p50'] > 4.0) &
    #                  (stats['Av_p50'] > 1.0))
    #print(indxs)

    # open 1D PDF file
    pdf1d_hdu = fits.open(filebase+'_pdf1d.fits')

    # filters for PHAT
    #filters = ['HST_WFC3_F225W', 'HST_WFC3_F275W', 'HST_WFC3_F336W', 
    #           'HST_ACS_WFC_F475W','HST_ACS_WFC_F550M', 
    #           'HST_ACS_WFC_F658N', 'HST_ACS_WFC_F814W',
    #           'HST_WFC3_F110W', 'HST_WFC3_F160W']
    #waves = np.asarray([2250., 2750.0, 3360.0, 
    #                    4750., 5500., 6580., 8140.,
    #                    11000., 16000.])
    if surveyname == 'METAL': 
        filters = ['HST_WFC3_F225W','HST_WFC3_F275W','HST_WFC3_F336W',
                   'HST_ACS_WFC_F475W','HST_ACS_WFC_F814W',
                   'HST_WFC3_F110W','HST_WFC3_F160W']
        waves = np.asarray([2250.0, 2722.05531502, 3366.00507206,
                            4763.04670013,8087.36760191,
                            11672.35909295,15432.7387546])
    else:
        filters = ['HST_WFC3_F275W','HST_WFC3_F336W','HST_ACS_WFC_F475W',
                   'HST_ACS_WFC_F814W','HST_WFC3_F110W','HST_WFC3_F160W']
        waves = np.asarray([2722.05531502, 3366.00507206,4763.04670013,
                            8087.36760191,11672.35909295,15432.7387546])
    
    fig, ax = pyplot.subplots(figsize=(12,12))

    # make the plot!
    plot_beast_ifit(filters, waves, stats, pdf1d_hdu, fontsize,
                    surveyname=surveyname)

    # show or save
    basename = filebase + '_ifit_starnum_' + str(starnum)
    print(basename)
    if args.png:
        fig.savefig(basename+'.png')
    elif args.eps:
        fig.savefig(basename+'.eps')
    elif args.pdf:
        fig.savefig(basename+'.pdf')
    else:
        pyplot.show()
    
