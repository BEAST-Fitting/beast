import matplotlib.pyplot as plt

from beast.observationmodel.noisemodel import toothpick
from beast.physicsmodel.grid import SEDGrid


def plot_toothpick_details(asts_filename, seds_filename, savefig=False):
    """
    Plot the details of the toothpick noisemodel creation including the
    individal AST results and the binned results

    Parameters
    ----------
    asts_filename : str
        filename with the AST results

    seds_filename : str
        filename with the SED grid (used just for the filter information)

    savefig : str (default=False)
        to save the figure, set this to the file extension (e.g., 'png', 'pdf')
    """
    sedgrid = SEDGrid(seds_filename, backend="cache")

    # read in AST results
    model = toothpick.MultiFilterASTs(asts_filename, sedgrid.filters)

    # set the column mappings as the external file is BAND_VEGA or BAND_IN
    model.set_data_mappings(upcase=True, in_pair=("in", "in"), out_pair=("out", "rate"))

    # compute binned biases, uncertainties, and completeness as a function of band flux
    model.fit_bins(nbins=20, completeness_mag_cut=-10)

    nfilters = len(sedgrid.filters)
    fig, ax = plt.subplots(nrows=nfilters, figsize=(8, 8), sharex=True)

    for i, cfilter in enumerate(sedgrid.filters):
        mag_in = model.data[model.filter_aliases[cfilter + "_in"]]
        flux_out = model.data[model.filter_aliases[cfilter + "_out"]]

        flux_in = (10 ** (-0.4 * mag_in)) * model.vega_flux[i]
        flux_out *= model.vega_flux[i]

        gvals = flux_out != 0.0

        ax[i].plot(
            flux_in[gvals],
            (flux_in[gvals] - flux_out[gvals]) / flux_in[gvals],
            "ko",
            alpha=0.1,
            markersize=2,
        )

        gmods = model._compls[:, i] > 0.0
        ax[i].plot(
            model._fluxes[gmods, i],
            model._biases[gmods, i] / model._fluxes[gmods, i],
            "b-",
        )

        ax[i].set_ylim(-1e2, 1e2)
        ax[i].set_xscale("log")
        ax[i].set_ylabel(r"$(F_i - F_o)/F_i$")

    ax[nfilters - 1].set_xlabel(r"$F_i$")
    ax[0].set_xlim(1e-21, 1e-13)

    # figname
    basename = asts_filename.replace(".fits", "_plot")

    fig.tight_layout()

    # save or show fig
    if savefig:
        fig.savefig("{}.{}".format(basename, savefig))
    else:
        plt.show()
