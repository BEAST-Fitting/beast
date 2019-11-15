import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits


def plot(
    beast_stats_file,
    param_list=["Av", "Rv", "logA", "f_A", "M_ini", "Z"],  # , 'distance'],
    n_bins=200,
):
    """
    Make a plot of each parameter vs the parameter errors

    Parameters with 'M_' (indicating mass) will have log10 taken to help with
    axis scaling

    Parameters
    ----------
    beast_stats_file : str
        path+file for the BEAST fitting results

    param_list : list of strings
        names of the parameters to plot (default is the params in Fig 16/17
        in Gordon+16)

    n_bins : int (default=200)
        number of bins to use in each dimension of the 2D histogram

    """

    # read in data
    with fits.open(beast_stats_file) as hdu:
        stats_table = hdu[1].data

    n_param = len(param_list)

    # figure
    fig = plt.figure(figsize=(10, 4 * n_param / 2))

    # make plots
    for p, param in enumerate(param_list):

        # subplot region
        ax = plt.subplot(n_param / 2, 2, p + 1)

        plot_x = stats_table[param + "_p50"]
        plot_y = (stats_table[param + "_p84"] - stats_table[param + "_p16"]) / 2

        if "M_" in param:
            plot_y = plot_y / (plot_x * np.log(10))
            plot_x = np.log10(plot_x)

        # plot
        plt.hist2d(
            plot_x, plot_y, bins=n_bins, cmap="magma", norm=matplotlib.colors.LogNorm()
        )

        # axis labels
        ax.tick_params(axis="both", which="major", labelsize=13)
        # ax.set_xlim(ax.get_xlim()[::-1])
        param_label = param
        if "M_" in param:
            param_label = "log " + param
        plt.xlabel(param_label, fontsize=14)
        plt.ylabel(r"$\sigma$ " + param_label, fontsize=14)

    plt.tight_layout()

    fig.savefig(beast_stats_file.replace(".fits", "_param_err.png"))
    plt.close(fig)
