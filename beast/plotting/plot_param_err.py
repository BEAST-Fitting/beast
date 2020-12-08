import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
from scipy.stats import binned_statistic_2d as stat2d
from matplotlib.ticker import MultipleLocator

__all__ = ["plot_param_err"]


def plot_param_err(
    beast_stats_file,
    param_list=["Av", "Rv", "logA", "f_A", "M_ini", "Z", "logT", "logg", "logL"],
    n_bins=200,
    colormap="cubehelix",
):
    """
    Make a plot of each parameter vs the parameter errors

    Parameters with M (indicating mass) will have log10 taken to help with
    axis scaling

    Parameters
    ----------
    beast_stats_file : str
        path+file for the BEAST fitting results

    param_list : list of strings
        names of the parameters to plot
        default is the params in Figs 16-18 in Gordon+16)

    n_bins : int (default=200)
        number of bins to use in each dimension of the 2D histogram

    colormap : str
        name of a colormap to use included with Matplotlib
    """

    # read in data
    with fits.open(beast_stats_file) as hdu:
        stats_table = hdu[1].data

    n_param = len(param_list)
    cmap = plt.get_cmap(colormap)

    # figure
    fig = plt.figure(figsize=(10, 30))

    # make plots
    for p, param in enumerate(param_list):

        # first column subplot
        ax1 = plt.subplot(n_param, 2, p * 2 + 1)

        param_p50 = stats_table[param + "_p50"]
        param_unc = (stats_table[param + "_p84"] - stats_table[param + "_p16"]) / 2

        if "M_" in param:
            param_p50 = np.log10(param_p50)
            param_unc = param_unc / (param_p50 * np.log(10))

        # plot
        plt.hist2d(param_p50, param_unc, bins=n_bins, cmap=cmap, norm=LogNorm())

        # axis labels
        plt.tick_params(axis="both", which="major", labelsize=13)
        ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
        # ax.set_xlim(ax.get_xlim()[::-1])
        param_label = param
        if "M_" in param:
            param_label = "log " + param
        plt.xlabel(param_label, fontsize=14)
        plt.ylabel(r"$\sigma$ " + param_label, fontsize=14)

        # second column subplot
        ax2 = plt.subplot(n_param, 2, p * 2 + 2)

        logT_p50 = stats_table["logT_p50"]
        logL_p50 = stats_table["logL_p50"]

        # plot
        h = stat2d(logT_p50, logL_p50, param_unc, bins=n_bins)
        plt.imshow(
            h[0].T,
            origin="lower",
            cmap=cmap,
            extent=[h[1].min(), h[1].max(), h[2].min(), h[2].max()],
            aspect="auto",
        )
        cbar = plt.colorbar()
        cbar.set_label(r"$\sigma$ " + param_label, fontsize=14)
        ax2.invert_xaxis()

        # axis labels
        plt.tick_params(axis="both", which="major", labelsize=13)
        plt.xlabel(r"log(T$_{\rm eff}$)", fontsize=14)
        plt.ylabel(r"log(L)", fontsize=14)

    plt.tight_layout()

    fig.savefig(beast_stats_file.replace(".fits", "_param_unc.png"))
    plt.close(fig)
