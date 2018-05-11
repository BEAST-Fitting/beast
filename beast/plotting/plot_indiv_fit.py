#!/usr/bin/env python
"""
Plot the individual fit for a single observed star

.. history::
    Written 12 Jan 2016 by Karl D. Gordon
      based on code written by Heddy Arab for the BEAST techniques paper figure
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
import matplotlib

from astropy.table import Table
from astropy.io import fits

from astropy import units as ap_units
from astropy.coordinates import SkyCoord as ap_SkyCoord

from beastplotlib import initialize_parser


def disp_str(stats, k, keyname):
    dvals = [stats[keyname + '_p50'][k],
             stats[keyname + '_p84'][k],
             stats[keyname + '_p16'][k]]
    if keyname == 'M_ini':
        dvals = np.log10(dvals)
    if keyname == 'distance':
        if dvals[0] > 1000:
            dvals = [v  / 1000. for v in dvals]
    disp_str = '$' + \
               "{0:.2f}".format(dvals[0]) + \
               '^{+' + \
               "{0:.2f}".format(dvals[1] - dvals[0]) + \
               '}_{-' + \
               "{0:.2f}".format(dvals[0] - dvals[2]) + \
               '}$'

    return disp_str


def plot_1dpdf(ax, pdf1d_hdu, tagname, xlabel, starnum,
               stats=None, logx=False):

    pdf = pdf1d_hdu[tagname].data

    n_objects, n_bins = pdf.shape
    n_objects -= 1

    ax.text(0.95, 0.95, xlabel, transform=ax.transAxes,
            va='top', ha='right')

    if n_bins == 1:
        ax.text(0.5, 0.5, 'unused', transform=ax.transAxes, va='center', ha='center')
        return

    xvals = pdf[n_objects, :]
    if logx:
        xvals = np.log10(xvals)

    if tagname == 'Z':
        gindxs, = np.where(pdf[starnum, :] > 0.)
        ax.plot(xvals[gindxs], pdf[starnum, gindxs] / max(pdf[starnum, gindxs]),
                color='k')
    else:
        ax.plot(xvals, pdf[starnum, :] / max(pdf[starnum, :]), color='k')

    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim = [xvals.min(), xvals.max()]
    xlim_delta = xlim[1] - xlim[0]
    ax.set_xlim(xlim[0] - 0.05 * xlim_delta, xlim[1] + 0.05 * xlim_delta)
    # ax.set_ylim(0.0,1.1*pdf[starnum,:].max())
    ax.set_ylim(0.0, 1.1)

    if stats is not None:
        ylim = ax.get_ylim()

        y1 = ylim[0] + 0.5 * (ylim[1] - ylim[0])
        y2 = ylim[0] + 0.7 * (ylim[1] - ylim[0])
        pval = stats[tagname + '_Best'][starnum]
        if logx:
            pval = np.log10(pval)
        ax.plot(np.full((2), pval), [y1, y2],
                '-', color='c')

        y1 = ylim[0] + 0.2 * (ylim[1] - ylim[0])
        y2 = ylim[0] + 0.4 * (ylim[1] - ylim[0])
        y1m = ylim[0] + 0.25 * (ylim[1] - ylim[0])
        y2m = ylim[0] + 0.35 * (ylim[1] - ylim[0])
        ym = 0.5 * (y1 + y2)
        pvals = [stats[tagname + '_p50'][starnum],
                 stats[tagname + '_p16'][starnum],
                 stats[tagname + '_p84'][starnum]]
        if logx:
            pvals = np.log10(pvals)
        ax.plot(np.full((2), pvals[0]), [y1m, y2m], '-', color='m')
        ax.plot(np.full((2), pvals[1]), [y1, y2], '-', color='m')
        ax.plot(np.full((2), pvals[2]), [y1, y2], '-', color='m')
        ax.plot(pvals[1:3], [ym, ym], '-', color='m')


def plot_beast_ifit(filters, waves, stats, pdf1d_hdu):

    # setup the plot grid
    gridNrow, gridNcol = 5, 12
    gs = gridspec.GridSpec(gridNrow, gridNcol,
                           height_ratios=[1.] * gridNrow,
                           width_ratios=[1.] * gridNcol)
    ax = []

    # axes for the big SED plot. Leave empty columns right of the plot to
    # put the legend and values.
    sed_height = 2
    free_cols = 2
    index_sedplot = len(ax)
    ax.append(plt.subplot(gs[0:sed_height, 0:-1 - free_cols]))

    # axes for the 1D PDFs
    nprim = 4
    nsec = 3
    nderiv = 3

    indices_1dpdf = []
    rows = [sed_height + i for i in range(3)]
    widths = [3, 4, 4]
    naxes = [nprim, nsec, nderiv]
    for r, w, n in zip(rows, widths, naxes):
        for i in range(n):
            indices_1dpdf.append(len(ax))
            ax.append(plt.subplot(gs[r, i * w:(i + 1) * w]))

    # plot the SED
    # print(np.sort(stats.colnames))

    n_filters = len(filters)

    # get the observations
    waves *= 1e-4
    obs_flux = np.empty((n_filters), dtype=np.float)
    mod_flux = np.empty((n_filters, 3), dtype=np.float)
    mod_flux_nd = np.empty((n_filters, 3), dtype=np.float)
    mod_flux_wbias = np.empty((n_filters, 3), dtype=np.float)
    k = starnum

    c = ap_SkyCoord(ra=stats['RA'][k] * ap_units.degree,
                    dec=stats['DEC'][k] * ap_units.degree,
                    frame='icrs')
    corname = ('PHAT J' +
               c.ra.to_string(unit=ap_units.hourangle, sep="", precision=2,
                              alwayssign=False, pad=True) +
               c.dec.to_string(sep="", precision=2,
                               alwayssign=True, pad=True))

    for i, cfilter in enumerate(filters):
        obs_flux[i] = stats[cfilter][k]
        mod_flux[i, 0] = np.power(10.0, stats['log' + cfilter + '_wd_p50'][k])
        mod_flux[i, 1] = np.power(10.0, stats['log' + cfilter + '_wd_p16'][k])
        mod_flux[i, 2] = np.power(10.0, stats['log' + cfilter + '_wd_p84'][k])
        mod_flux_nd[i, 0] = np.power(
            10.0, stats['log' + cfilter + '_nd_p50'][k])
        mod_flux_nd[i, 1] = np.power(
            10.0, stats['log' + cfilter + '_nd_p16'][k])
        mod_flux_nd[i, 2] = np.power(
            10.0, stats['log' + cfilter + '_nd_p84'][k])
        if 'log' + cfilter + '_wd_bias_p50' in stats.colnames:
            mod_flux_wbias[i, 0] = np.power(10.0, stats['log' + cfilter +
                                                        '_wd_bias_p50'][k])
            mod_flux_wbias[i, 1] = np.power(10.0, stats['log' + cfilter +
                                                        '_wd_bias_p16'][k])
            mod_flux_wbias[i, 2] = np.power(10.0, stats['log' + cfilter +
                                                        '_wd_bias_p84'][k])
    sed_ax = ax[index_sedplot]
    sed_ax.plot(waves, obs_flux, 'ko', label='observed')

    if 'log' + filters[0] + '_wd_bias_p50' in stats.colnames:
        sed_ax.plot(waves, mod_flux_wbias[:, 0],
                    'b-', label='stellar+dust+bias')
        sed_ax.fill_between(waves, mod_flux_wbias[:, 1], mod_flux_wbias[:, 2],
                            color='b', alpha=0.3)

    sed_ax.plot(waves, mod_flux[:, 0], 'r-', label='stellar+dust')
    sed_ax.fill_between(waves, mod_flux[:, 1], mod_flux[:, 2],
                        color='r', alpha=0.2)

    sed_ax.plot(waves, mod_flux_nd[:, 0], 'y-', label='stellar only')
    sed_ax.fill_between(waves, mod_flux_nd[:, 1], mod_flux_nd[:, 2],
                        color='y', alpha=0.1)

    sed_ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1.025))

    sed_ax.set_ylabel(r'Flux [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')
    sed_ax.set_yscale('log')

    sed_ax.set_xscale('log')
    sed_ax.text(0.5, -0.01, r'$\lambda$ [$\AA$]',
                transform=sed_ax.transAxes, va='top')
    sed_ax.set_xlim(0.2, 2.0)
    sed_ax.set_xticks([0.2, 0.3, 0.4, 0.5, 0.8, 0.9, 1.0, 2.0])
    sed_ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    sed_ax.text(0.05, 0.95, corname, transform=sed_ax.transAxes,
                va='top', ha='left')

    # add the text results
    keys = ['Av', 'M_ini', 'logA', 'distance',
            'Rv', 'f_A', 'Z', 'logT', 'logg', 'logL']
    dispnames = ['A(V)', 'log(M)', 'log(t)', 'distance(pc)', 'R(V)', r'f$_\mathcal{A}$', 'Z',
                 r'log(T$_\mathrm{eff})$', 'log(g)', 'log(L)']
    startprim, stopprim = 0, nprim - 1  # 0 1 2 3
    startsec, stopsec = stopprim + 1, stopprim + nsec  # 4 5 6
    startderiv, stopderiv = stopsec + 1, stopsec + nderiv  # 7 8 9
    laby = 0.72
    ty = np.linspace(laby - 0.07, 0.1, num=len(keys))
    ty[startsec:] -= 0.025
    ty[startderiv:] -= 0.025
    tx = [1.12, 1.2, 1.3]
    for i in range(len(keys)):
        sed_ax.text(tx[0], ty[i], dispnames[i],
                    ha='right',
                    transform=sed_ax.transAxes)
        sed_ax.text(tx[1], ty[i], disp_str(stats, starnum, keys[i]),
                    ha='center', color='m',
                    transform=sed_ax.transAxes)
        best_val = stats[keys[i] + '_Best'][k]
        if keys[i] == 'M_ini':
            best_val = np.log10(best_val)
        if keys[i] == 'distance':
            best_val /= 1000.
            dispnames[i] = dispnames[i].replace('pc','kpc')
        sed_ax.text(tx[2], ty[i],
                    '$' + "{0:.2f}".format(best_val) + '$',
                    ha='center', color='c',
                    transform=sed_ax.transAxes)
    sed_ax.text(tx[0], laby, 'Param',
                ha='right',
                transform=sed_ax.transAxes)
    sed_ax.text(tx[1], laby, '50%$\pm$33%',
                ha='center', color='k',
                transform=sed_ax.transAxes)
    sed_ax.text(tx[2], laby, 'Best', color='k',
                ha='center',
                transform=sed_ax.transAxes)

    # now draw boxes around the different kinds of parameters
    tax = sed_ax
    left, right = tx[0], tx[-1]

    def draw_box_around_values(start, stop, ls):
        deltaline = ty[start] - ty[start + 1]
        top = ty[start] + deltaline  # Draw the top border ABOVE the text
        bottom = ty[stop]
        rec = Rectangle((left - 0.1, bottom - 0.02), right - left + 0.15, top - bottom + 0.01,
                        fill=False, lw=2, transform=tax.transAxes, ls=ls)
        rec = tax.add_patch(rec)
        rec.set_clip_on(False)

    # primary
    draw_box_around_values(startprim, stopprim, ls='dashed')

    # secondary
    draw_box_around_values(startsec, stopsec, ls='dotted')

    # derived
    draw_box_around_values(startderiv, stopderiv, ls='dashdot')

    # padding for rectangles of 1D PDFs
    pad = 0.1

    # Make these plots:

    # A, M, t, dist,
    # R, fA, Z
    # logT, logg, logL

    # This is done by iterating over the axes created at the start of
    # this function, from left to right, line per line.

    # plot the primary parameter 1D PDFs
    ax_iter = (ax[i] for i in indices_1dpdf)
    first_primary_ax = next(ax_iter)
    plot_1dpdf(first_primary_ax, pdf1d_hdu, 'Av', 'A(V)', starnum,
               stats=stats)
    plot_1dpdf(next(ax_iter), pdf1d_hdu, 'M_ini', 'log(M)', starnum, logx=True,
               stats=stats)
    plot_1dpdf(next(ax_iter), pdf1d_hdu, 'logA', 'log(t)', starnum,
               stats=stats)
    last_primary_ax = next(ax_iter)
    plot_1dpdf(last_primary_ax, pdf1d_hdu, 'distance', 'distance(pc)', starnum,
               stats=stats)

    # plot the secondary parameter 1D PDFs
    first_secondary_ax = next(ax_iter)
    plot_1dpdf(first_secondary_ax, pdf1d_hdu, 'Rv', 'R(V)', starnum,
               stats=stats)
    plot_1dpdf(next(ax_iter), pdf1d_hdu, 'f_A', r'f$_\mathcal{A}$', starnum,
               stats=stats)
    last_secondary_ax = next(ax_iter)
    plot_1dpdf(last_secondary_ax, pdf1d_hdu, 'Z', 'Z', starnum,
               stats=stats)

    # plot the derived parameter 1D PDFs
    first_derived_ax = next(ax_iter)
    plot_1dpdf(first_derived_ax, pdf1d_hdu, 'logT', r'log(T$_\mathrm{eff})$', starnum,
               stats=stats)
    plot_1dpdf(next(ax_iter), pdf1d_hdu, 'logg', 'log(g)', starnum,
               stats=stats)
    last_derived_ax = next(ax_iter)
    plot_1dpdf(last_derived_ax, pdf1d_hdu, 'logL', 'log(L)', starnum,
               stats=stats)

    # A more manual version of tight_layout
    plt.subplots_adjust(top=.95, bottom=.05, left=.125, right=.925,
                        wspace=.5, hspace=.5)

    # PLOT ALL THE BOXES AFTER CALLING TIGHT LAYOUT! Tight layout
    # changes the coordinates of the axes a little, but leaves the boxes
    # untouched. Therefore, we plot the boxes here by extracting the
    # coordinates of the axes after they have been modified by
    # tight_layout.

    def rectangle_around_axes(bottomleft_ax, topright_ax, pad, ls, label=None):
        """
        pad: tuple, (left, right, bottom, top)
        """
        left, bottom = bottomleft_ax.get_position().get_points()[0]
        right, top = topright_ax.get_position().get_points()[1]
        left -= pad[0]
        right += pad[1]
        bottom -= pad[2]
        top += pad[3]
        transf = plt.gcf().transFigure
        rec = Rectangle((left, bottom), right - left, top - bottom,
                        transform=transf, fill=False, lw=2, ls=ls)
        rec = bottomleft_ax.add_patch(rec)
        rec.set_clip_on(False)

        if label:
            middle = (top + bottom) / 2.
            moreleft = left  # pad[0]
            bottomleft_ax.text(moreleft, middle, label, transform=transf,
                               rotation='vertical', fontstyle='oblique',
                               va='center', ha='right')

    rectanglePadding = (0.03, 0.01, 0.03, 0.01)

    # Box around primaries
    tax = first_primary_ax
    rectangle_around_axes(first_primary_ax, last_primary_ax,
                          pad=rectanglePadding, ls='dashed',
                          label='Primary')
    tax.text(0.0, 0.5, 'Probability', transform=tax.transAxes,
             rotation='vertical',
             va='center', ha='right')

    # Box around secondaries
    tax = first_secondary_ax
    rectangle_around_axes(first_secondary_ax, last_secondary_ax,
                          pad=rectanglePadding, ls='dotted',
                          label='Secondary')
    tax.text(0.0, 0.5, 'Probability', transform=tax.transAxes,
             rotation='vertical',
             va='center', ha='right')

    # Box around deriveds
    tax = first_derived_ax
    rectangle_around_axes(first_derived_ax, last_derived_ax,
                          pad=rectanglePadding, ls='dashdot',
                          label='Derived')
    tax.text(0.0, 0.5, 'Probability', transform=tax.transAxes,
             rotation='vertical',
             va='center', ha='right')

if __name__ == '__main__':

    parser = initialize_parser()
    parser.add_argument("filebase", type=str,
                        help='base filename of output results')
    parser.add_argument("--starnum", type=int, default=0,
                        help="star number in observed file")
    args = parser.parse_args()

    starnum = args.starnum

    # base filename
    filebase = args.filebase

    # read in the stats
    stats = Table.read(filebase + '_stats.fits')

    # open 1D PDF file
    pdf1d_hdu = fits.open(filebase + '_pdf1d.fits')

    # filters for PHAT
    # filters = ['HST_WFC3_F225W', 'HST_WFC3_F275W', 'HST_WFC3_F336W',
    #           'HST_ACS_WFC_F475W','HST_ACS_WFC_F550M',
    #           'HST_ACS_WFC_F658N', 'HST_ACS_WFC_F814W',
    #           'HST_WFC3_F110W', 'HST_WFC3_F160W']
    # waves = np.asarray([2250., 2750.0, 3360.0,
    #                    4750., 5500., 6580., 8140.,
    #                    11000., 16000.])
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    waves = np.asarray([2722.05531502, 3366.00507206, 4763.04670013,
                        8087.36760191, 11672.35909295, 15432.7387546])

    fig, ax = plt.subplots(figsize=(8, 8))

    # make the plot!
    plot_beast_ifit(filters, waves, stats, pdf1d_hdu)

    # show or save
    basename = filebase + '_ifit_starnum_' + str(starnum)
    print(basename)
    if args.savefig:
        fig.savefig('{}.{}'.format(basename, args.savefig))
    else:
        plt.show()
