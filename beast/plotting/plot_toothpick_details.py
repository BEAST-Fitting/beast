import matplotlib.pyplot as plt
import tables

from astropy.table import Table

from beast.observationmodel.noisemodel import toothpick
from beast.plotting.beastplotlib import set_params
from beast.config import __ROOT__

__all__ = ["plot_toothpick_details"]


def plot_toothpick_details(asts_filename, savefig=False):
    """
    Plot the details of the toothpick noisemodel creation for each filter.
    These plots show the individual AST results as points as
    (flux_in - flux_out)/flux_in.  In addition, the binned values of these
    points are plotted giving the bias term in the observation model.
    Error bars around the binned bias values give the binned sigma term of
    the observation model.  Finally, as a separate column of plots the
    binned completeness in each filter is plotted.

    Parameters
    ----------
    asts_filename : str
        filename with the AST results

    savefig : str (default=False)
        to save the figure, set this to the file extension (e.g., 'png', 'pdf')
    """
    # determine filters
    filters = fetch_filters(asts_filename)

    # read in AST results
    model = toothpick.MultiFilterASTs(asts_filename, filters)

    # set the column mappings as the external file is BAND_VEGA or BAND_IN
    model.set_data_mappings(upcase=True, in_pair=("in", "in"), out_pair=("out", "rate"))

    # compute binned biases, uncertainties, and completeness as a function of band flux
    ast_nonrecovered_ratio = 2.0
    model.fit_bins(
        nbins=50,
        ast_nonrecovered_ratio=ast_nonrecovered_ratio,
    )

    nfilters = len(filters)
    figsize_y = nfilters * 3
    fig, ax = plt.subplots(nrows=nfilters, ncols=2, figsize=(14, figsize_y), sharex=True)
    set_params()

    for i, cfilter in enumerate(filters):
        mag_in = model.data[model.filter_aliases[cfilter + "_in"]]
        flux_out = model.data[model.filter_aliases[cfilter + "_out"]]

        flux_in = (10 ** (-0.4 * mag_in)) * model.vega_flux[i]
        flux_out *= model.vega_flux[i]

        gvals = flux_out != 0.0

        ax[i, 0].plot(
            flux_in[gvals],
            flux_out[gvals] / flux_in[gvals],
            "ko",
            alpha=0.1,
            markersize=2,
        )

        # not all bins are filled with good data
        ngbins = model._nasts[i]

        ax[i, 0].plot(
            model._fluxes[0:ngbins, i],
            1. + model._biases[0:ngbins, i] / model._fluxes[0:ngbins, i],
            "b-",
        )

        ax[i, 0].errorbar(
            model._fluxes[0:ngbins, i],
            1. + model._biases[0:ngbins, i] / model._fluxes[0:ngbins, i],
            yerr=model._sigmas[0:ngbins, i] / model._fluxes[0:ngbins, i],
            fmt="bo",
            markersize=2,
            alpha=0.5,
        )

        if ast_nonrecovered_ratio is not None:
            ax[i, 0].axhline(
                ast_nonrecovered_ratio, linestyle="--", alpha=0.25, color="k"
            )

        ax[i, 0].set_ylim(-10, 2.5)
        ax[i, 0].set_ylabel(r"$F_o/F_i$")

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

    # add in the zero line
    # do after all the data has been plotted to get the full x range
    pxrange = ax[0, 0].get_xlim()
    for i, cfilter in enumerate(filters):
        ax[i, 0].plot(pxrange, [1.0, 1.0], "k--", alpha=0.5)

    # figname
    basename = asts_filename.replace(".fits", "_plot")

    fig.tight_layout()

    # save or show fig
    if savefig:
        fig.savefig("{}.{}".format(basename, savefig))
    else:
        plt.show()


def fetch_filters(filename):
    """
    For any FITS table file (gst, ast, etc.), collect the filters based on the
    columns containing "VEGA". Returns a sorted list of the full BEAST filter
    names, assuming the camera is HST_WFC3.

    Parameters
    ----------
    filename : str
        filename containing VEGA measurements

    Returns
    -------
    filters : list
        list of full BEAST filter names

    """

    data = Table.read(filename)

    # extract every filter mentioned in the table
    # find all columns mentioning VEGA
    filters = [f.split("_")[0] for f in data.columns if "VEGA" in f]

    # sort the list of filters based on wavelength
    # i.e. anything < 200nm is in the IR
    rep = 0
    for n in range(len(filters)):
        if "L" in filters[n - rep]:
            num = (int(filters[n - rep][-5:-2]))

        else:
            num = (int(filters[n - rep][-4:-1]))

        if num < 200:
            filters.append(filters.pop(filters.index('{0}'.format(filters[n - rep]))))
            rep += 1

    # add camera prefix
    filters = ["HST_WFC3_" + f for f in filters]

    # check to see that the filters actually exist
    source = "{0}/vega.hd5".format(__ROOT__)
    hdf = tables.open_file(source)

    beast_filters = hdf.root.sed.cols.FNAME[:]
    beast_filters = [ele.decode() for ele in beast_filters]
    hdf.close()

    missing_filters = set(filters).difference(beast_filters)

    if missing_filters:
        print("Missing HST WFC3 equivalent filter for {0}".format(missing_filters))

    return filters
