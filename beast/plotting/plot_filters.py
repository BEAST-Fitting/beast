#!/usr/bin/env python
"""
Make a nice plot of the filter response functions
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from beast.observationmodel import phot
from beast.plotting.beastplotlib import initialize_parser

__all__ = ["plot_filters"]


def plot_filters(
    filter_names,
    filterLib=None,
    save_name="beast_filters",
    xlim=[1.4e3, 2e4],
    ylim=[1e-4, 2],
    show_plot=True,
):

    """Plots transmission curves in log-log space.

    Parameters
    ----------
    filter_names : list of str
        List of full names of filters to plot
    filterLib : str, optional
        Filter file (None=default)
    save_name : str, optional
        Filename to save plot as
    xlim : length 2 list
        Values to set plot x-limits to
    ylim : length 2 list
        Values to set plot y-limits to
    """
    if not isinstance(filter_names, list):
        filter_names = [filter_names]

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # wavelength grid in angstroms for response functions
    waves = np.logspace(3, np.log10(3e4), 501)

    # read in the filter response functions
    flist = phot.load_filters(
        filter_names, interp=True, lamb=waves, filterLib=filterLib
    )

    color_indices = np.log10(np.array(np.sort([f.norm for f in flist])))
    color_indices -= color_indices.min()
    color_indices /= color_indices.max()

    cmap = mpl.cm.plasma
    # ax.set_prop_cycle(color=[cmap(i) for i in color_indices])
    color = iter(cmap(np.linspace(0.2, 0.8, len(filter_names))))

    for f in flist:
        c = next(color)
        ax.plot(f.wavelength, f.transmit, color=c, lw=2)
        ax.fill_between(f.wavelength, f.transmit, alpha=0.2, color=c)
        ax.text(
            np.nanmean(f.wavelength[f.transmit > 100.0 * ylim[0]]),
            1.3 * np.nanmax(f.transmit[f.transmit > ylim[0]]),
            f.name.split("_")[-1],
            ha="center",
            color=c,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(r"$\lambda$ [$\mu m$]")
    ax.set_ylabel(r"$B_i(\lambda)$")

    # ax.set_xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    fig.tight_layout()

    if show_plot:
        plt.show()
    else:
        return fig


if __name__ == "__main__":  # pragma: no cover
    parser = initialize_parser()
    parser.add_argument(
        "--save_name",
        action="store",
        default="filters.appendVegaFilter",
        help="Save figure to file",
    )
    args = parser.parse_args()

    filter_names = [
        "HST_WFC3_F225W",
        "HST_WFC3_F275W",
        "HST_WFC3_F336W",
        "HST_ACS_WFC_F475W",
        "HST_ACS_WFC_F550M",
        "HST_ACS_WFC_F814W",
        "HST_WFC3_F110W",
        "HST_WFC3_F160W",
    ]

    fig = plot_filters(filter_names, show_plot=False)

    if args.tex:
        plt.rc({"usetex": True})
    if args.savefig:
        fig.savefig("{}.{}".format(args.save_name, args.savefig))
    else:
        plt.show()
