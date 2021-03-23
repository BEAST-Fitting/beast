import numpy as np
import matplotlib.pyplot as plt
import argparse
import re

from beast.physicsmodel.grid import SEDGrid
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel

__all__ = ["plot_noisemodel"]


def plot_noisemodel(
    sed_file,
    noise_file_list,
    plot_file,
    samp=100,
    cmap_name='viridis',
):
    """
    Make a plot of the noise model: for each of the bandsm make plots of bias
    and uncertainty as a function of flux

    If there are multiple files in noise_file_list, each of them will be
    overplotted in each panel.

    Parameters
    ----------
    sed_file : string
        path+name of the SED grid file

    noise_file_list : list of strings
        path+name of the noise model file(s)

    plot_file : string
        name of the file to save the plot

    samp : int (default=100)
        plotting all of the SED points takes a long time for a viewer to load,
        so set this to plot every Nth point

    cmap_name : string (default=plt.cm.viridis)
        name of a color map to use
    """

    # read in the SED grid
    print("* reading SED grid file")
    sed_object = SEDGrid(sed_file)
    if hasattr(sed_object.seds, "read"):
        sed_grid = sed_object.seds.read()
    else:
        sed_grid = sed_object.seds
    filter_list = sed_object.filters
    n_filter = len(filter_list)

    # figure
    fig, ax = plt.subplots(nrows=3, ncols=n_filter, figsize=(25, 15))

    # setup the plots
    fontsize = 12
    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    plt.set_cmap(cmap_name)

    # go through noise files after sorting them according to
    # their SD bin number
    noise_file_list.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    bin_label = [re.findall(r"bin\d+", x)[0] for x in noise_file_list]
    for n, nfile in enumerate(np.atleast_1d(noise_file_list)):

        print("* reading " + nfile)

        # read in the values
        noisemodel_vals = noisemodel.get_noisemodelcat(nfile)

        # extract error and bias
        noise_err = noisemodel_vals["error"]
        noise_bias = noisemodel_vals["bias"]
        noise_compl = noisemodel_vals["completeness"]

        # plot things
        for f, filt in enumerate(filter_list):

            # error is negative where it's been extrapolated -> trim those
            good_err = np.where(noise_err[:, f] > 0)[0]
            plot_sed = sed_grid[good_err, f][::samp]
            plot_err = noise_err[good_err, f][::samp]
            plot_bias = noise_bias[good_err, f][::samp]
            plot_compl = noise_compl[good_err, f][::samp]

            # bias
            bax = ax[0, f]
            bax.plot(
                np.log10(plot_sed),
                plot_bias / plot_sed,
                marker="o",
                linestyle="none",
                mew=0,
                ms=2,
                alpha=0.1,
                label='SD %s' % (bin_label[n]),
            )

            bax.tick_params(axis="both", which="major")
            bax.set_xlabel("log " + filt)
            bax.set_ylabel(r"Bias ($\mu$/F)")
            leg = bax.legend(loc='lower right', markerscale=3)
            for lh in leg.legendHandles:
                lh._legmarker.set_alpha(1)

            # error
            eax = ax[1, f]
            eax.plot(
                np.log10(plot_sed),
                plot_err / plot_sed,
                marker="o",
                linestyle="none",
                mew=0,
                ms=2,
                alpha=0.1,
            )

            eax.tick_params(axis="both", which="major")
            eax.set_xlabel("log " + filt)
            eax.set_ylabel(r"Error ($\sigma$/F)")

            # completeness
            cax = ax[2, f]
            cax.plot(
                np.log10(plot_sed),
                plot_compl,
                marker="o",
                linestyle="none",
                mew=0,
                ms=2,
                alpha=0.1,
            )

            cax.tick_params(axis="both", which="major")
            cax.set_xlabel("log " + filt)
            cax.set_ylabel(r"Completeness")

    plt.tight_layout()

    fig.savefig(plot_file, dpi=300)
    plt.close(fig)


if __name__ == "__main__":  # pragma: no cover

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sed_file", type=str, help="path+name of the sed grid file",
    )
    parser.add_argument(
        "noise_file_list",
        type=str,
        nargs="+",
        help="path+name of the noise model file(s)",
    )
    parser.add_argument(
        "plot_file", type=str, help="name of the file to save the plot",
    )
    parser.add_argument(
        "--samp", type=int, default=100, help="plot every Nth point",
    )
    parser.add_argument(
        "--cmap_name",
        type=str,
        default="viridis",
        help="color map to use when making plots",
    )

    args = parser.parse_args()

    plot_noisemodel(
        args.sed_file,
        args.noise_file_list,
        args.plot_file,
        samp=args.samp,
        cmap_name=args.cmap_name,
    )

    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()
