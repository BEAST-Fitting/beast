# plot_cmd.py
# Plots a generic CMD from real or simulated BEAST fitting data
# Petia YMJ
# Created 9/13/18
# Updated 10/05/18
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from functools import reduce

from beast.plotting.beastplotlib import initialize_parser, set_params


__all__ = ["plot_cmd"]


def plot_cmd(
    fitsfile,
    mag1_filter="F475W",
    mag2_filter="F814W",
    mag3_filter="F475W",
    savefig=False,
    show_plot=True,
):
    """
    Read in flux from real or simulated data in fitsfile and plot a
    color-magnitude diagram based on specified filters.

    Parameters
    ----------
    fitsfile : str
        input fitsfile (includes full path to file); format = .fits
    mag1_filter : str (default='F475W')
        1st color filter
    mag2_filter : str (default='F814W')
        2nd color filter
    mag3_filter : str (default='F475W')
        magnitude
    savefig : str (default=False)
        to save the figure, set this to the file extension (e.g., 'png', 'pdf')
    show_plot : boolean
        True, show the plot (to screen or a file)
        False, return the fig
    """

    fits_data = fits.open(fitsfile)
    table = fits_data[1].data
    fits_data.close()

    # Read in band_rate
    mag1_flux = table["%s" % (mag1_filter + "_rate")]
    mag2_flux = table["%s" % (mag2_filter + "_rate")]
    mag_flux = table["%s" % (mag3_filter + "_rate")]

    # Exclude negative or 0 fluxes
    m1_pos_inds = np.where(mag1_flux > 0.0)
    m2_pos_inds = np.where(mag2_flux > 0.0)
    m_pos_inds = np.where(mag_flux > 0.0)
    pos_inds = reduce(np.intersect1d, (m1_pos_inds, m2_pos_inds, m_pos_inds))
    mag1_flux_pos = mag1_flux[pos_inds]
    mag2_flux_pos = mag2_flux[pos_inds]
    mag_flux_pos = mag_flux[pos_inds]

    # Convert from flux to mags
    mag1 = (-2.5) * np.log10(mag1_flux_pos)
    mag2 = (-2.5) * np.log10(mag2_flux_pos)
    mag = (-2.5) * np.log10(mag_flux_pos)

    col = mag1 - mag2

    fig = plt.figure(figsize=(9, 9))

    # sets up plots to have larger fonts, wider lines, etc.
    set_params()

    plt.plot(col, mag, "k.", ms=0.1)

    plt.gca().invert_yaxis()

    plt.xlabel("%s - %s" % (mag1_filter, mag2_filter))
    plt.ylabel(mag3_filter)
    plt.title(fitsfile)

    fig.tight_layout()

    # figname
    basename = fitsfile.replace(".fits", "_plot")

    # save or show fig
    if show_plot:
        if savefig:
            fig.savefig("{}.{}".format(basename, savefig))
        else:
            plt.show()
    else:
        return fig


if __name__ == "__main__":  # pragma: no cover

    parser = initialize_parser()
    parser.add_argument("filename", type=str, help="Path to FITS file to plot")
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
        "--mag3",
        action="store",
        default="F475W",
        help="Choose filter for the magnitude",
    )

    args = parser.parse_args()

    # plot the CMD
    fig = plot_cmd(
        args.filename,
        mag1_filter=args.mag1,
        mag2_filter=args.mag2,
        mag3_filter=args.mag3,
        savefig=args.savefig,
    )
