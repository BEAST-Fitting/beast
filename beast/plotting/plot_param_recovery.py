import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits

__all__ = ["plot_param_recovery"]


def plot_param_recovery(
    sim_data_list,
    stats_file_list,
    output_plot_filename,
    file_label_list=None,
    max_nbins=20,
):
    """
    Make plots comparing the physical parameters from simulated data to the
    recovered physical parameters

    If there are multiple files input, it is presumably because they are from
    different noise models.  If that's the case, you may want to assign labels
    for each of them (file_label_list).

    Parameters
    ----------
    sim_data_list : string or list of strings
        File(s) of simulated data from beast.tools.simulate_obs, which have both
        the photometry and physical parameters

    stats_file_list : string or list of strings
        File(s) of the corresponding stats files with the fit statistics

    output_plot_filename : string
        name of the file in which to save the output plot

    file_label_list : string (default=None)
        Labels to use for each of the files (e.g., their source density ranges)

    max_nbins : int (default=10)
        maximum number of bins to use in each dimension of the 2D histogram
        (fewer will be used if there are fewer unique values)

    """

    # parameters to plot
    param_list = ["Av", "logA", "M_ini", "Rv", "f_A", "Z", "distance"]
    n_param = len(param_list)

    # number of files
    n_stat = len(sim_data_list)

    # figure
    fig = plt.figure(figsize=(5 * n_stat, 4 * n_param))

    # iterate through the files
    for i, (sim_stats, recov_stats) in enumerate(
        zip(np.atleast_1d(sim_data_list), np.atleast_1d(stats_file_list))
    ):

        # read in data
        with fits.open(sim_stats) as hdu_sim, fits.open(recov_stats) as hdu_recov:
            sim_table = hdu_sim[1].data
            recov_table = hdu_recov[1].data

        # make plots
        for p, param in enumerate(param_list):

            # subplot region
            ax = plt.subplot(n_param, n_stat, 1 + n_stat * p + i)

            # set things to plot
            plot_x = sim_table[param]
            plot_y = recov_table[param + "_p50"]

            if ("M_" in param) or (param == "Z"):
                plot_x = np.log10(plot_x)
                plot_y = np.log10(plot_y)

            # number of bins
            n_uniq = len(np.unique(plot_x))
            n_bins = [min(n_uniq, max_nbins), min(3 * n_uniq, max_nbins)]

            # plot
            plt.hist2d(
                plot_x,
                plot_y,
                bins=n_bins,
                cmap="magma",
                norm=matplotlib.colors.LogNorm(),
            )

            # axis labels
            ax.tick_params(axis="both", which="major", labelsize=13)
            # ax.set_xlim(ax.get_xlim()[::-1])
            param_label = param
            if ("M_" in param) or (param == "Z"):
                param_label = "log " + param
            plt.xlabel("Simulated " + param_label, fontsize=14)
            plt.ylabel("Recovered " + param_label, fontsize=14)

    plt.tight_layout()

    fig.savefig(output_plot_filename)
    plt.close(fig)
