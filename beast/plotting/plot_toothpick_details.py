# import numpy as np
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
    model.fit_bins(nbins=50, completeness_mag_cut=-10)

    nfilters = len(sedgrid.filters)
    fig, ax = plt.subplots(nrows=nfilters, ncols=2, figsize=(12, 8), sharex=True)

    for i, cfilter in enumerate(sedgrid.filters):
        mag_in = model.data[model.filter_aliases[cfilter + "_in"]]
        flux_out = model.data[model.filter_aliases[cfilter + "_out"]]

        flux_in = (10 ** (-0.4 * mag_in)) * model.vega_flux[i]
        flux_out *= model.vega_flux[i]

        gvals = flux_out != 0.0

        # delt = np.absolute(((flux_in[gvals] - flux_out[gvals]) / flux_in[gvals]) - 1.0)
        # gvals2 = delt < 0.01
        # print(cfilter, flux_in[gvals][gvals2][0], flux_out[gvals][gvals2][0], delt[gvals2][0])

        ax[i, 0].plot(
            flux_in[gvals],
            (flux_in[gvals] - flux_out[gvals]) / flux_in[gvals],
            "ko",
            alpha=0.1,
            markersize=2,
        )

        # not all bins are filled with good data
        ngbins = model._nasts[i]
        ax[i, 0].plot(
            model._fluxes[0:ngbins, i],
            model._biases[0:ngbins, i] / model._fluxes[0:ngbins, i],
            "b-",
        )

        ax[i, 0].errorbar(
            model._fluxes[0:ngbins, i],
            model._biases[0:ngbins, i] / model._fluxes[0:ngbins, i],
            yerr=model._sigmas[0:ngbins, i] / model._fluxes[0:ngbins, i],
            fmt="bo",
            markersize=2,
            alpha=0.5,
        )

        ax[i, 0].set_ylim(-5e0, 5e0)
        ax[i, 0].set_xscale("log")
        ax[i, 0].set_ylabel(r"$(F_i - F_o)/F_i$")

        ax[i, 1].plot(
            model._fluxes[0:ngbins, i],
            model._compls[0:ngbins, i],
            "b-",
        )

        ax[i, 1].yaxis.tick_right()
        ax[i, 1].yaxis.set_label_position("right")
        ax[i, 1].set_ylim(0, 1)
        ax[i, 1].set_xscale("log")
        sfilt = cfilter.split("_")[-1]
        ax[i, 1].set_ylabel(f"C({sfilt})")

    ax[nfilters - 1, 0].set_xlabel(r"$F_i$")
    ax[nfilters - 1, 1].set_xlabel(r"$F_i$")
    # ax[0, 0].set_xlim(1e-25, 1e-13)

    # figname
    basename = asts_filename.replace(".fits", "_plot")

    fig.tight_layout()

    # save or show fig
    if savefig:
        fig.savefig("{}.{}".format(basename, savefig))
    else:
        plt.show()
