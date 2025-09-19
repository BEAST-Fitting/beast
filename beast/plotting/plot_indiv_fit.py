#!/usr/bin/env python
"""
Plot the individual fit for a single observed star
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
import matplotlib

from astropy.table import Table
from astropy.io import fits

from beast.plotting.beastplotlib import initialize_parser
from beast.tools.symlog import inverse_symlog, symlog_linthreshold


__all__ = ["plot_indiv_fit"]


def disp_str(stats, k, keyname):
    dvals = [
        stats[keyname + "_p50"][k],
        stats[keyname + "_p84"][k],
        stats[keyname + "_p16"][k],
    ]
    if keyname in ["M_ini", "Z"]:
        dvals = np.log10(dvals)
    if keyname == "distance":
        if dvals[0] > 1000:
            dvals = [v / 1000.0 for v in dvals]
    disp_str = (
        "$"
        + "{0:.2f}".format(dvals[0])
        + "^{+"
        + "{0:.2f}".format(dvals[1] - dvals[0])
        + "}_{-"
        + "{0:.2f}".format(dvals[0] - dvals[2])
        + "}$"
    )

    return disp_str


def plot_1dpdf(ax, pdf1d_hdu, tagname, xlabel, starnum, stats=None, logx=False):

    pdf_data = pdf1d_hdu[tagname].data

    if pdf_data.ndim == 2:
        pdf = pdf_data[starnum, :]
        xvals = pdf_data[-1, :]
        n_objects, n_bins = pdf_data.shape
        n_objects -= 1
    elif pdf_data.ndim == 3:
        pdf = pdf_data[starnum, :, 0]
        xvals = pdf_data[starnum, :, 1]
        n_bins = np.sum(~np.isnan(xvals))

    ax.text(0.95, 0.95, xlabel, transform=ax.transAxes, va="top", ha="right")

    if (n_bins == 1) or (n_bins == 0):
        ax.text(0.5, 0.5, "unused", transform=ax.transAxes, va="center", ha="center")
        ax.set_yticklabels([])
        return

    if logx:
        xvals = np.log10(xvals)

    # there is a problem when "Z" in the model grid is set to a single value
    #  it shows up as two values here
    #  likely an issue deep in the BEAST fit.py code - distance handled better
    # if tagname == "Z":
    #     print(xvals, pdf)
    #     (gindxs,) = np.where(pdf > 0.0)
    #     ax.plot(xvals[gindxs], pdf[gindxs] / max(pdf[gindxs]), color="k")
    # else:
    # if tagname == "Z":
    #    print(xvals, pdf)
    ax.plot(xvals, pdf / max(pdf), color="k")

    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    xlim = [xvals.min(), xvals.max()]
    xlim_delta = xlim[1] - xlim[0]
    if ~np.isnan(xlim[0]):
        ax.set_xlim(xlim[0] - 0.05 * xlim_delta, xlim[1] + 0.05 * xlim_delta)
    else:
        bestval = stats[tagname + "_Best"][starnum]
        if tagname == "distance":
            bestval /= 1000.0
        ax.set_xlim(0.95 * bestval, 1.05 * bestval)
    ax.set_ylim(0.0, 1.1)
    ax.set_yticklabels([])

    if stats is not None:
        ylim = ax.get_ylim()

        y1 = ylim[0] + 0.5 * (ylim[1] - ylim[0])
        y2 = ylim[0] + 0.7 * (ylim[1] - ylim[0])
        pval = stats[tagname + "_Best"][starnum]
        if tagname == "distance":
            pval /= 1000.0
        if logx:
            pval = np.log10(pval)
        ax.plot(np.full((2), pval), [y1, y2], "-", color="c")

        y1 = ylim[0] + 0.2 * (ylim[1] - ylim[0])
        y2 = ylim[0] + 0.4 * (ylim[1] - ylim[0])
        y1m = ylim[0] + 0.25 * (ylim[1] - ylim[0])
        y2m = ylim[0] + 0.35 * (ylim[1] - ylim[0])
        ym = 0.5 * (y1 + y2)
        pvals = [
            stats[tagname + "_p50"][starnum],
            stats[tagname + "_p16"][starnum],
            stats[tagname + "_p84"][starnum],
        ]
        if logx:
            pvals = np.log10(pvals)
        ax.plot(np.full((2), pvals[0]), [y1m, y2m], "-", color="m")
        ax.plot(np.full((2), pvals[1]), [y1, y2], "-", color="m")
        ax.plot(np.full((2), pvals[2]), [y1, y2], "-", color="m")
        ax.plot(pvals[1:3], [ym, ym], "-", color="m")


def plot_indiv_fit(filebase, starnum=0, savefig=False, plotfig=True):
    """
    Plot the individual fit for a single observed star including best fit &
    percentile parameters and various 1D pPDFs

    Parameters
    ----------
    filebase : str
        base filename of run

    starnum : int
        number of star in the stats file

    savefig : str
        set to the file extension fo the desired plot file (e.g., png, pdf, etc)

    plotfig : boolean
        plot the figure to a file or the screen based on savefig
        otherwise return the fig object
    """

    starnum = int(starnum)

    # determine how the stats/pdf1d filenames are to be set
    if len(np.atleast_1d(filebase)) == 1:
        stats_fname = f"{filebase}_stats.fits"
        pdf1d_fname = f"{filebase}_pdf1d.fits"
    else:
        stats_fname = filebase[0]
        pdf1d_fname = filebase[1]

    # read in the stats
    stats = Table.read(stats_fname, hdu=1)

    # check how many extensions the stats file has
    # determines how to get the filternames and wavelengths
    with fits.open(stats_fname) as hdul:
        nhdu = len(hdul)
    if nhdu > 2:
        filter_info = Table.read(stats_fname, hdu=2)
        bfilters = filter_info["filternames"].data
        waves = filter_info["wavelengths"].data
        filters = [cfilter.decode("utf-8") for cfilter in bfilters]
    else:  # PHAT values as default to support old stats files
        filters = [
            "HST_WFC3_F275W",
            "HST_WFC3_F336W",
            "HST_ACS_WFC_F475W",
            "HST_ACS_WFC_F814W",
            "HST_WFC3_F110W",
            "HST_WFC3_F160W",
        ]
        waves = np.asarray([2722.05, 3366.01, 4763.05, 8087.37, 11672.36, 15432.74])

    # open 1D PDF file
    pdf1d_hdu = fits.open(pdf1d_fname)

    fig, ax = plt.subplots(figsize=(8, 8))

    # setup the plot grid
    gridNrow, gridNcol = 5, 12
    gs = gridspec.GridSpec(
        gridNrow,
        gridNcol,
        height_ratios=[1.0] * gridNrow,
        width_ratios=[1.0] * gridNcol,
    )
    ax = []

    # axes for the big SED plot. Leave empty columns right of the plot to
    # put the legend and values.
    sed_height = 2
    free_cols = 3
    index_sedplot = len(ax)
    ax.append(plt.subplot(gs[0:sed_height, 0 : -1 - free_cols]))

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
            ax.append(plt.subplot(gs[r, i * w : (i + 1) * w]))

    # plot the SED
    n_filters = len(filters)

    # get the observations
    waves *= 1e-4
    obs_flux = np.zeros((n_filters), dtype=float)
    mod_flux = np.zeros((n_filters, 3), dtype=float)
    mod_flux_nd = np.zeros((n_filters, 3), dtype=float)
    mod_flux_wbias = np.zeros((n_filters, 3), dtype=float)
    k = starnum

    corname = stats["Name"][k]

    for i, cfilter in enumerate(filters):
        obs_flux[i] = stats[cfilter][k]
        fluxname = "log" + cfilter
        mod_flux[i, 0] = np.power(10.0, stats[fluxname + "_wd_p50"][k])
        mod_flux[i, 1] = np.power(10.0, stats[fluxname + "_wd_p16"][k])
        mod_flux[i, 2] = np.power(10.0, stats[fluxname + "_wd_p84"][k])
        mod_flux_nd[i, 0] = np.power(10.0, stats[fluxname + "_nd_p50"][k])
        mod_flux_nd[i, 1] = np.power(10.0, stats[fluxname + "_nd_p16"][k])
        mod_flux_nd[i, 2] = np.power(10.0, stats[fluxname + "_nd_p84"][k])
        if "sym" + fluxname + "_wd_bias_p50" in stats.colnames:
            mod_flux_wbias[i, 0] = inverse_symlog(
                stats["sym" + fluxname + "_wd_bias_p50"][k]
            )
            mod_flux_wbias[i, 1] = inverse_symlog(
                stats["sym" + fluxname + "_wd_bias_p16"][k]
            )
            mod_flux_wbias[i, 2] = inverse_symlog(
                stats["sym" + fluxname + "_wd_bias_p84"][k]
            )

    sed_ax = ax[index_sedplot]
    sed_ax.plot(waves, obs_flux, "ko", label="observed")

    if "symlog" + filters[0] + "_wd_bias_p50" in stats.colnames:
        sed_ax.plot(waves, mod_flux_wbias[:, 0], "b-", label="stellar+dust+bias")
        sed_ax.fill_between(
            waves, mod_flux_wbias[:, 1], mod_flux_wbias[:, 2], color="b", alpha=0.3
        )

    sed_ax.plot(waves, mod_flux[:, 0], "r-", label="stellar+dust")
    sed_ax.fill_between(waves, mod_flux[:, 1], mod_flux[:, 2], color="r", alpha=0.2)

    sed_ax.plot(waves, mod_flux_nd[:, 0], "y-", label="stellar only")
    sed_ax.fill_between(
        waves, mod_flux_nd[:, 1], mod_flux_nd[:, 2], color="y", alpha=0.1
    )

    # can introduce a legend loc option if 'best' produces overlap
    sed_ax.legend(loc='best', fontsize=9)

    sed_ax.set_ylabel(r"Flux [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]")
    sed_ax.set_yscale("symlog", linthresh=symlog_linthreshold)
    sed_ax.grid(True)

    sed_ax.text(0.5, -0.07, r"$\lambda$ [$\mu m$]", transform=sed_ax.transAxes, va="top")
    sed_ax.set_xlim(0.2, 2.0)
    sed_ax.set_xscale("log")
    sed_ax.minorticks_off()
    sed_ax.set_xticks([0.2, 0.3, 0.4, 0.5, 0.8, 1.0, 2.0])
    sed_ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    sed_ax.get_xaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())

    sed_ax.text(0.05, 0.95, corname, transform=sed_ax.transAxes, va="top", ha="left")

    # add the text results
    keys = ["Av", "M_ini", "logA", "distance", "Rv", "f_A", "Z", "logT", "logg", "logL"]
    dispnames = [
        "A(V)",
        "log(M)",
        "log(t)",
        "d(kpc)",
        "R(V)",
        r"f$_\mathcal{A}$",
        "log(Z)",
        r"log(T$_\mathrm{eff})$",
        "log(g)",
        "log(L)",
    ]
    startprim, stopprim = 0, nprim - 1  # 0 1 2 3
    startsec, stopsec = stopprim + 1, stopprim + nsec  # 4 5 6
    startderiv, stopderiv = stopsec + 1, stopsec + nderiv  # 7 8 9
    laby = 0.96
    ty = np.linspace(laby - 0.1, 0.1, num=len(keys))
    ty[startsec:] -= 0.04
    ty[startderiv:] -= 0.04
    tx = [1.14, 1.3, 1.47]
    for i in range(len(keys)):
        sed_ax.text(tx[0], ty[i], dispnames[i], ha="center", transform=sed_ax.transAxes)
        sed_ax.text(
            tx[1],
            ty[i],
            disp_str(stats, starnum, keys[i]),
            ha="center",
            color="m",
            transform=sed_ax.transAxes,
        )
        best_val = stats[keys[i] + "_Best"][k]
        if keys[i] in ["M_ini", "Z"]:
            print(best_val)
            best_val = np.log10(best_val)
        if keys[i] == "distance":
            best_val /= 1000.0
            dispnames[i] = dispnames[i].replace("pc", "kpc")
        sed_ax.text(
            tx[2],
            ty[i],
            "$" + "{0:.2f}".format(best_val) + "$",
            ha="center",
            color="c",
            transform=sed_ax.transAxes,
        )
    sed_ax.text(
        tx[0], laby, "Param", ha="center", transform=sed_ax.transAxes, fontsize=10
    )
    sed_ax.text(
        tx[1],
        laby,
        r"50$\pm$33%",
        ha="center",
        color="k",
        transform=sed_ax.transAxes,
        fontsize=10,
    )
    sed_ax.text(
        tx[2],
        laby,
        "Best",
        color="k",
        ha="center",
        transform=sed_ax.transAxes,
        fontsize=10,
    )

    # now draw boxes around the different kinds of parameters
    tax = sed_ax
    left, right = tx[0], tx[-1]

    def draw_box_around_values(start, stop, ls):
        deltaline = ty[start] - ty[start + 1]
        top = ty[start] + deltaline  # Draw the top border ABOVE the text
        bottom = ty[stop]
        rec = Rectangle(
            (left - 0.1, bottom - 0.02),
            right - left + 0.15,
            top - bottom + 0.01,
            fill=False,
            lw=2,
            transform=tax.transAxes,
            ls=ls,
        )
        rec = tax.add_patch(rec)
        rec.set_clip_on(False)

    # primary
    draw_box_around_values(startprim, stopprim, ls="dashed")

    # secondary
    draw_box_around_values(startsec, stopsec, ls="dotted")

    # derived
    draw_box_around_values(startderiv, stopderiv, ls="dashdot")

    # Make these plots:

    # A, M, t, dist,
    # R, fA, Z
    # logT, logg, logL

    # This is done by iterating over the axes created at the start of
    # this function, from left to right, line per line.

    # plot the primary parameter 1D PDFs
    ax_iter = (ax[i] for i in indices_1dpdf)
    first_primary_ax = next(ax_iter)
    plot_1dpdf(first_primary_ax, pdf1d_hdu, "Av", "A(V)", starnum, stats=stats)
    plot_1dpdf(
        next(ax_iter), pdf1d_hdu, "M_ini", "log(M)", starnum, logx=True, stats=stats
    )
    plot_1dpdf(next(ax_iter), pdf1d_hdu, "logA", "log(t)", starnum, stats=stats)
    last_primary_ax = next(ax_iter)
    plot_1dpdf(last_primary_ax, pdf1d_hdu, "distance", "d(kpc)", starnum, stats=stats)

    # plot the secondary parameter 1D PDFs
    first_secondary_ax = next(ax_iter)
    plot_1dpdf(first_secondary_ax, pdf1d_hdu, "Rv", "R(V)", starnum, stats=stats)
    plot_1dpdf(
        next(ax_iter), pdf1d_hdu, "f_A", r"f$_\mathcal{A}$", starnum, stats=stats
    )
    last_secondary_ax = next(ax_iter)
    plot_1dpdf(last_secondary_ax, pdf1d_hdu, "Z", "log(Z)", starnum, logx=True, stats=stats)

    # plot the derived parameter 1D PDFs
    first_derived_ax = next(ax_iter)
    plot_1dpdf(
        first_derived_ax,
        pdf1d_hdu,
        "logT",
        r"log(T$_\mathrm{eff})$",
        starnum,
        stats=stats,
    )
    plot_1dpdf(next(ax_iter), pdf1d_hdu, "logg", "log(g)", starnum, stats=stats)
    last_derived_ax = next(ax_iter)
    plot_1dpdf(last_derived_ax, pdf1d_hdu, "logL", "log(L)", starnum, stats=stats)

    # A more manual version of tight_layout
    plt.subplots_adjust(
        top=0.95, bottom=0.05, left=0.125, right=0.925, wspace=0.5, hspace=0.5
    )

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
        rec = Rectangle(
            (left, bottom),
            right - left,
            top - bottom,
            transform=transf,
            fill=False,
            lw=2,
            ls=ls,
        )
        rec = bottomleft_ax.add_patch(rec)
        rec.set_clip_on(False)

        if label:
            middle = (top + bottom) / 2.0
            moreleft = left  # pad[0]
            bottomleft_ax.text(
                moreleft,
                middle,
                label,
                transform=transf,
                rotation="vertical",
                fontstyle="oblique",
                va="center",
                ha="right",
            )

    rectanglePadding = (0.03, 0.01, 0.03, 0.01)

    # Box around primaries
    tax = first_primary_ax
    rectangle_around_axes(
        first_primary_ax,
        last_primary_ax,
        pad=rectanglePadding,
        ls="dashed",
        label="Primary",
    )
    tax.text(
        0.0,
        0.5,
        "Probability",
        transform=tax.transAxes,
        rotation="vertical",
        va="center",
        ha="right",
    )

    # Box around secondaries
    tax = first_secondary_ax
    rectangle_around_axes(
        first_secondary_ax,
        last_secondary_ax,
        pad=rectanglePadding,
        ls="dotted",
        label="Secondary",
    )
    tax.text(
        0.0,
        0.5,
        "Probability",
        transform=tax.transAxes,
        rotation="vertical",
        va="center",
        ha="right",
    )

    # Box around deriveds
    tax = first_derived_ax
    rectangle_around_axes(
        first_derived_ax,
        last_derived_ax,
        pad=rectanglePadding,
        ls="dashdot",
        label="Derived",
    )
    tax.text(
        0.0,
        0.5,
        "Probability",
        transform=tax.transAxes,
        rotation="vertical",
        va="center",
        ha="right",
    )

    # show or save
    if plotfig:
        basename = filebase + "_ifit_starnum_" + str(starnum)
        if savefig:
            fig.savefig("{}.{}".format(basename, savefig))
        else:
            plt.show()
    else:
        return fig


if __name__ == "__main__":  # pragma: no cover

    parser = initialize_parser()
    parser.add_argument("filebase", type=str, help="base filename of output results")
    parser.add_argument(
        "--starnum", type=int, default=0, help="star number in observed file"
    )
    args = parser.parse_args()

    # make the plot!
    plot_indiv_fit(args.filebase, args.starnum, args.savefig)
