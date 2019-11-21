import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import copy


def plot_triangle(
    beast_stats_file, param_list=["Av", "Rv", "logA", "f_A", "M_ini", "Z", "distance"],
):
    """
    Make a triangle/corner plot of everything vs everything

    Parameters
    ----------
    beast_stats_file : str
        path+file for the BEAST fitting results

    param_list : list of strings
        names of the parameters to plot (default is the params in Fig 16/17
        in Gordon+16)

    """

    n_params = len(param_list)

    # read in data
    with fits.open(beast_stats_file) as hdu:
        stats_table = hdu[1].data

    # figure
    fig = plt.figure(figsize=(4 * n_params, 4 * n_params))

    # iterate through the panels
    for i, pi in enumerate(param_list):
        for j, pj in enumerate(param_list[i:], i):

            # set x/y axes
            plot_x = stats_table[pi + "_p50"]
            plot_y = stats_table[pj + "_p50"]
            x_label = copy.copy(pi)
            y_label = copy.copy(pj)
            # take log if necessary
            if ("M_" in pi) or (pi == "Z"):
                plot_x = np.log10(plot_x)
                x_label = "log " + pi
            if ("M_" in pj) or (pj == "Z"):
                plot_y = np.log10(plot_y)
                y_label = "log " + pj

            # not along diagonal
            if i != j:

                # set up subplot
                plt.subplot(n_params, n_params, i + j * (n_params) + 1)
                ax = plt.gca()

                # plot points
                plt.plot(
                    plot_x,
                    plot_y,
                    marker="o",
                    mew=0,
                    color="black",
                    markersize=2,
                    linestyle="None",
                    alpha=0.2,
                )

                ax.tick_params(
                    axis="both",
                    which="both",
                    direction="in",
                    labelsize=14,
                    bottom=True,
                    top=True,
                    left=True,
                    right=True,
                )

                # axis labels and ticks
                if i == 0:
                    ax.set_ylabel(y_label, fontsize=16)
                    # ax.get_yaxis().set_label_coords(-0.35,0.5)
                else:
                    ax.set_yticklabels([])
                if j == n_params - 1:
                    ax.set_xlabel(x_label, fontsize=16)
                    plt.xticks(rotation=-45)
                else:
                    ax.set_xticklabels([])

            # along diagonal
            if i == j:

                # set up subplot
                plt.subplot(n_params, n_params, i + j * (n_params) + 1)
                ax = plt.gca()

                # make histogram
                plt.hist(
                    plot_x, bins=20, facecolor="grey", linewidth=0.25, edgecolor="grey"
                )

                ax.tick_params(axis="y", which="both", length=0, labelsize=14)
                ax.tick_params(axis="x", which="both", direction="in", labelsize=14)

                # axis labels and ticks
                ax.set_yticklabels([])
                if i < n_params - 1:
                    ax.set_xticklabels([])
                if i == n_params - 1:
                    ax.set_xlabel(x_label, fontsize=16)
                    plt.xticks(rotation=-45)

    # plt.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.tight_layout()

    fig.savefig(beast_stats_file.replace(".fits", "_param_triangle.pdf"))
    plt.close(fig)
