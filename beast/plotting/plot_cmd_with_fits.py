import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import copy

from beast.plotting.beastplotlib import initialize_parser

__all__ = ["plot_cmd_with_fits"]


def plot_cmd_with_fits(
    data_fits_file,
    beast_stats_file,
    mag1_filter="F475W",
    mag2_filter="F814W",
    mag3_filter="F475W",
    param="chi2min",
    log_param=False,
    plot_all=False,
    savefig=False,
    show_plot=True,
):
    """
    Make a CMD with the data, and color-code points by some other fitted quantity

    Magnitudes are calculated as -2.5*log10(filter_RATE), rather than directly
    extracting the magnitude from the catalog.

    Parameters
    ----------
    data_fits_file : str
        Path+file for the stellar photometry. Photometry will be matched to
        sources in beast_stats_file using RA/Dec, so this can contain sources
        that were not fit with the BEAST.

    beast_stats_file : str
        path+file for the BEAST fitting results

    mag1_filter : str (default='F475W')
        1st color filter (color=mag1-mag2)

    mag2_filter : str (default='F814W')
        2nd color filter (color=mag1-mag2)

    mag3_filter : str (default='F475W')
        filter for the magnitude

    param : str (default='chi2min')
        parameter to use for color-coding points

    log_param : boolean (default=False)
        choose whether to take the log of `param` for assigning color

    plot_all : boolean (default=False)
        If True, plot all points by converting the fluxes into magnitudes.
        If False, only plot sources with Vega mags that are <99 in the
        mag1/mag2/mag3 filters

    savefig : str (default=False)
        to save the figure, set this to the file extension (e.g., 'png', 'pdf')

    show_plot : boolean
        True, show the plot (to screen or a file)
        False, return the fig
    """

    # read in data
    with fits.open(data_fits_file) as data_hdu:
        data_table = data_hdu[1].data
    with fits.open(beast_stats_file) as beast_hdu:
        beast_table = beast_hdu[1].data

    # figure out the subset that were modeled
    data_cat = SkyCoord(
        ra=data_table["RA"] * u.degree, dec=data_table["Dec"] * u.degree
    )
    beast_cat = SkyCoord(
        ra=beast_table["RA"] * u.degree, dec=beast_table["Dec"] * u.degree
    )
    ind, sep, _ = beast_cat.match_to_catalog_sky(data_cat)
    data_table = data_table[ind]

    # Read in band_rate
    mag1_flux = data_table[f"{mag1_filter}_rate"]
    mag2_flux = data_table[f"{mag2_filter}_rate"]
    mag_flux = data_table[f"{mag3_filter}_rate"]

    # read in parameter for color-coding
    color_data = beast_table[param]

    # choose whether to plot all or some of the pointss
    if plot_all:
        # exclude negative or 0 fluxes
        good_ind = np.where((mag1_flux > 0.0) & (mag2_flux > 0.0) & (mag_flux > 0.0))[0]
    else:
        # exclude any that have mag = 0 (indicating either no coverage in
        # that band or sub-optimal sharpness+concentration values)
        temp = [
            copy.copy(data_table[filt + "_VEGA"])
            for filt in set([mag1_filter, mag2_filter, mag3_filter])
        ]
        for col in temp:
            col[col > 99] = np.nan
        good_ind = np.where(np.isfinite(np.sum(temp, axis=0)))[0]

    mag1_flux_pos = mag1_flux[good_ind]
    mag2_flux_pos = mag2_flux[good_ind]
    mag_flux_pos = mag_flux[good_ind]
    color_data_pos = color_data[good_ind]

    # take log of param if set
    if log_param:
        color_data_pos = np.log10(color_data_pos)

    # Convert from flux to mags
    mag1 = (-2.5) * np.log10(mag1_flux_pos)
    mag2 = (-2.5) * np.log10(mag2_flux_pos)
    mag = (-2.5) * np.log10(mag_flux_pos)

    col = mag1 - mag2

    # do the plotting

    fig = plt.figure(figsize=(7, 6))

    im = plt.scatter(
        col,
        mag,
        c=color_data_pos,
        marker="o",
        s=2,
        edgecolors="none",
        cmap="viridis_r",
        alpha=0.25,
        vmin=np.percentile(color_data_pos, 1),
        vmax=np.percentile(color_data_pos, 99),
    )
    ax = plt.gca()

    ax.set_xlim((np.percentile(col, 0.01), np.percentile(col, 99.99)))
    ax.set_ylim((np.percentile(mag, 0.01), np.percentile(mag, 99.99)))

    plt.gca().invert_yaxis()
    plt.xlabel(f"{mag1_filter} - {mag2_filter}", fontsize=15)
    plt.ylabel(mag3_filter, fontsize=15)
    ax.tick_params(axis="both", labelsize=13)

    cbar = plt.colorbar(im)
    cbar.solids.set(alpha=1)
    # cbar = ax.figure.colorbar(color_data_pos, ax=ax)
    cbar_label = param
    if log_param:
        cbar_label = "Log " + param
    cbar.ax.set_ylabel(cbar_label, fontsize=13)  # , rotation=-90, va="bottom")

    # save or show fig
    if show_plot:
        if savefig:
            basename = beast_stats_file.replace(".fits", f"_cmd_{param}")
            fig.savefig("{}.{}".format(basename, savefig))
        else:
            plt.show()
    else:
        return fig


if __name__ == "__main__":  # pragma: no cover

    parser = initialize_parser()
    parser.add_argument(
        "data_fits_file", type=str, help="Path to FITS file with stellar photometry"
    )
    parser.add_argument(
        "beast_stats_file", type=str, help="Path to FITS file with BEAST fits"
    )

    parser.add_argument(
        "--mag1",
        action="store",
        default="F475W",
        help="Choose filter for mag1 (color=mag1-mag2)",
    )
    parser.add_argument(
        "--mag2",
        action="store",
        default="F814W",
        help="Choose filter for mag2 (color=mag1-mag2)",
    )
    parser.add_argument(
        "--magy",
        action="store",
        default="F475W",
        help="Choose filter for the magnitude",
    )

    parser.add_argument(
        "--param",
        action="store",
        default="chi2min",
        help="Choose parameter to color-code the CMD",
    )
    parser.add_argument(
        "--log_param",
        action="store_true",
        help="Set this if you would like to take the log of `param` for assigning color",
    )

    args = parser.parse_args()

    # plot the CMD
    fig = plot_cmd_with_fits(
        args.data_fits_file,
        args.beast_stats_file,
        mag1_filter=args.mag1,
        mag2_filter=args.mag2,
        mag3_filter=args.magy,
        param=args.param,
        log_param=args.log_param,
        savefig=args.savefig,
        show_plot=True,
    )
