import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm

from beast.observationmodel.noisemodel import toothpick
from beast.physicsmodel.grid import SEDGrid

__all__ = ["plot_ast_diagnostics"]


def plot_ast_diagnostics(ast_file_list, seds_filename, interpolate=True, savefig=False):
    """
    Plot the details of the toothpick noisemodel creation for each filter.
    These plots show the bias, error/uncertainity, and completeness for
    each observational filter, and . If a list of ASTs divid

    Parameters
    ----------
    ast_file_list : list
        list with the filenames of all the AST results
        if there is just one filename, pass it as [ast_filename]
    seds_filename : str
        filename with the SED grid (used just for the filter information)
    interpolate : bool
        interpolate between the known flux bins in each source density bin with enough data
    savefig : str (default=False)
        to save the figure, set this to the file extension (e.g., 'png', 'pdf')
    """
    # fetch the number of filters needed
    sedgrid = SEDGrid(seds_filename, backend="cache")
    seds = sedgrid.seds
    nfilters = len(sedgrid.filters)

    # set up plot
    fig1, ax = plt.subplots(nrows=3, ncols=nfilters, figsize=(5 * nfilters + 1, 12))
    cmaps = plt.get_cmap("viridis_r")

    # figure out what the max source density bin is
    # for plotting the color later
    number_max = 1
    # if filename has bin info
    if "bin" in ast_file_list[0]:
        # go through each file
        for n, nfile in enumerate(ast_file_list):
            # figure out the source density number
            number = int(nfile.split(".")[-2].split("bin")[-1])
            # if it's larger than the previous max, replace it
            if number > number_max:
                number_max = number

    for n, nfile in enumerate(ast_file_list):
        asts_filename = nfile

        if len(ast_file_list) > 1:
            number = int(nfile.split(".")[-2].split("bin")[-1])

        # if just one file
        else:
            number = 1

        # read in AST results
        model = toothpick.MultiFilterASTs(asts_filename, sedgrid.filters)

        # set the column mappings as the external file is BAND_VEGA or BAND_IN
        model.set_data_mappings(
            upcase=True, in_pair=("in", "in"), out_pair=("out", "rate")
        )

        # compute binned biases, uncertainties, and completeness as a function of band flux
        ast_nonrecovered_ratio = 2.0
        model.fit_bins(
            nbins=50,
            ast_nonrecovered_ratio=ast_nonrecovered_ratio,
        )

        print(asts_filename)
        print("Source density : %i" % number)
        print("Observations : %i" % len(model.data))
        print(
            "Non-Zero Detections (%s): %i"
            % (
                model.filters[0].split("_")[-1] + "_VEGA",
                np.sum(model.data[model.filters[0].split("_")[-1] + "_VEGA"] != 99.999),
            )
        )

        if interpolate:
            # interpolate all files with more than 100 detections
            if (
                np.sum(model.data[model.filters[0].split("_")[-1] + "_VEGA"] != 99.999)
                > 100
            ):
                bias, sigma, comp = model.interpolate(sedgrid)
                print("Interpolate : Successful")

            else:
                print("Interpolate : NOT Successful \n")

        for i in range(nfilters):

            if interpolate:

                samp = 100

                good_err = np.where(sigma[:, i] > 0)[0]
                plot_sed = seds[good_err, i][::samp]  # only pulls every 100th point
                plot_err = sigma[good_err, i][::samp]
                plot_bias = bias[good_err, i][::samp]
                plot_comp = comp[good_err, i][::samp]

                marker = "."
                linestyle = "none"
                alpha = 0.1

            else:

                good_err = np.where(model._sigmas[:, i] > 0)[0]
                plot_sed = model._fluxes[good_err, i]
                plot_err = model._sigmas[good_err, i]
                plot_bias = model._biases[good_err, i]
                plot_comp = model._compls[good_err, i]

                marker = "o"
                linestyle = "-"
                alpha = 1

            # plot bias

            ax[0, i].set_yscale("symlog")
            ax[0, 0].set_ylabel(r"Abs Bias ($\mu$/F)", fontsize=10)

            ax[0, i].plot(
                np.log10(plot_sed),
                np.abs(plot_bias) / plot_sed,
                marker=marker,
                linestyle=linestyle,
                mew=0,
                ms=2,
                alpha=alpha,
                c=cmaps((number - 1) / (number_max - 1)),
                label=number,
            )

            # plot sigma/error (uncertainty)

            ax[1, i].set_yscale("log")
            ax[1, 0].set_ylabel(r"Error ($\sigma$/F)", fontsize=10)

            ax[1, i].plot(
                np.log10(plot_sed),
                plot_err / plot_sed,
                marker=marker,
                linestyle=linestyle,
                mew=0,
                ms=2,
                c=cmaps((number - 1) / (number_max - 1)),
                label=number,
                alpha=alpha,
            )

            # plot completeness

            ax[2, 0].set_ylabel(r"Completness", fontsize=10)

            ax[2, i].plot(
                np.log10(plot_sed),
                plot_comp,
                marker=marker,
                linestyle=linestyle,
                mew=0,
                ms=2,
                c=cmaps((number - 1) / (number_max - 1)),
                label=number,
                alpha=alpha,
            )

    # if there is more than 1 source density bin
    # aka number_max is different than the default
    if number_max > 1:

        # create colorbar
        cbar = fig1.colorbar(
            cm.ScalarMappable(cmap=cmaps),
            ax=ax.ravel().tolist(),
            ticks=np.arange(0, number_max, 1) / (number_max - 1),
        )
        cbar.ax.set_yticklabels(np.arange(1, number_max + 1, 1))
        cbar.ax.get_yaxis().labelpad = 15
        cbar.set_label("Source Density", rotation=270)

    else:
        plt.tight_layout()

    for n in range(nfilters):
        ax[2, n].set_xlabel(r"$%s$" % sedgrid.filters[n].split("_")[-1])

    # figname

    # if filename has bin in it
    if "bin" in ast_file_list[0]:
        basename = ast_file_list[0].split("_bin")[0] + "_plot"

    else:
        basename = ast_file_list[0].replace(".fits", "_plot")

    # save or show fig
    if savefig:
        fig1.savefig("{}.{}".format(basename, savefig))
    else:
        plt.show(fig1)
