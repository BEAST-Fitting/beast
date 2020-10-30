import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy import interpolate
from scipy.stats import gaussian_kde

from beast.tools import beast_settings
from beast.tools.density_map import BinnedDensityMap

__all__ = ["plot_ast_diagnostics"]


def plot_ast_diagnostics(
    beast_settings_info,
    ast_file,
    map_file,
    source_density_image,
    reference_image,
    filter_selections=["HST_WFC3_F475W", "HST_WFC3_F814W"],
    savefig=True,
):
    """
    Make a summary plot of the input ASTs. These include their spatial
    distribution relative to the source density image, and CMDs as a function of
    source density (to verify that the input ASTs are being correctly replicated
    across source density bins, and to see the color-magnitude distributions
    of the inputs).

    Output plot is saved in the same location/name as ast_file, but with a
    "_diagnostic.png" instead of ".txt".

    Parameters
    ----------
    beast_settings_info : string or instance
        if string: file name with beast settings
        if class: beast.tools.beast_settings.beast_settings instance

    ast_file : string (default=None)
        name of AST input file

    map_file : string (default=None)
        background or source density map file, used to created BinnedDensityMap
        object for replicating source density bins for AST inputs

    source_density_image : string
        name of the source density image FITS file (for plotting)

    reference_image : string
        name of the reference image FITS file, used for converting AST pixel
        coordinates to WCS.

    filter_selections : list (default=['HST_WFC3_F475W', 'HST_WFC3_F814W'])
        names of the two filters to be used a the reference color for the CMDs.
        For example, ['HST_WFC3_F475W', 'HST_WFC3_F814W'] will plot
        F475W-F814W colors for the CMDs.

    savefig : bool
        If True, save figure to file. If False, don't.
    """

    # process beast settings info
    if isinstance(beast_settings_info, str):
        settings = beast_settings.beast_settings(beast_settings_info)
    elif isinstance(beast_settings_info, beast_settings.beast_settings):
        settings = beast_settings_info
    else:
        raise TypeError(
            "beast_settings_info must be string or beast.tools.beast_settings.beast_settings instance"
        )

    # Read in the input AST file (the ascii .txt file that is sent to Ben for processing)
    ast_input = Table.read(ast_file, format="ascii")

    # Import filter information from the BEAST settings file
    filters = settings.filters.copy()
    # Count number of filters and decide how many rows to plot
    ncmds = len(filters)
    nrows = int(ncmds / 2) + 1

    # Make a giant summary plot
    fig = plt.figure(0, [16, nrows * 4])
    outer_grid = fig.add_gridspec(1, 2, wspace=0.5, hspace=0.4)

    # Add RA and Dec information to the input AST file (which is just an ascii filewith only X,Y positions)
    hdu_ref = fits.open(reference_image)
    wcs_ref = WCS(hdu_ref[0].header)
    source_astin = wcs_ref.wcs_pix2world(ast_input["X"], ast_input["Y"], 0)

    # Compute source coordinates in SD image frame
    hdu_sd = fits.open(source_density_image)
    wcs_sd = WCS(hdu_sd[0].header)
    source_sdin = wcs_sd.wcs_world2pix(source_astin[0], source_astin[1], 0)

    ### In the first column, plot the reference image and the AST coordinates:
    inner_grid = outer_grid[0].subgridspec(2, 1, wspace=0.4, hspace=0.5)
    ax = fig.add_subplot(inner_grid[0])  # , projection=wcs_sd)
    im = ax.imshow(hdu_sd[0].data, cmap="Greys", origin="lower")
    plt.colorbar(im, ax=ax, label=r"$\rm Source density$")
    ax.scatter(source_sdin[0], source_sdin[1], color="orange", marker=".", alpha=0.05)

    # Build density_map instance to replicate source density binning for the input AST catalog
    bdm = BinnedDensityMap.create(
        map_file, N_bins=settings.sd_Nbins, bin_width=settings.sd_binwidth
    )
    # Figure out what the bins are
    bin_foreach_source = np.zeros(len(ast_input), dtype=int)
    for i in range(len(ast_input)):
        bin_foreach_source[i] = bdm.bin_for_position(
            source_astin[0][i], source_astin[1][i]
        )
    # compute the AST input indices for each bin
    binnrs = np.unique(bin_foreach_source)
    bin_idxs = []
    for b in binnrs:
        sources_for_bin = np.where(bin_foreach_source == b)
        bin_idxs.append([sources_for_bin])

    ### In the second column, plot the CMDs of all the filters as a function of source density
    inner_grid = outer_grid[1].subgridspec(nrows, 2, wspace=0.4, hspace=0.5)

    # Loop through the number of CMDs (= number of filters)
    for j in range(ncmds):
        ax = fig.add_subplot(inner_grid[j])

        # Pick out reasonable plotting ranges (also used in contour code)
        kvals = [
            np.min(ast_input[filter_selections[0]] - ast_input[filter_selections[1]])
            - 2,
            np.max(ast_input[filter_selections[0]] - ast_input[filter_selections[1]])
            + 2,
            np.max(ast_input[filters[j]]) + 4,
            np.min(ast_input[filters[j]]) - 4,
        ]

        # Set axis labels
        ax.set_xlabel(
            filter_selections[0].split("_")[-1]
            + "-"
            + filter_selections[1].split("_")[-1]
        )
        ax.set_ylabel(filters[j].split("_")[-1])

        colors = iter(cm.magma(np.linspace(0.2, 0.8, len(binnrs))))
        for k, binnr in enumerate(binnrs):
            cat = ast_input[bin_idxs[k]]
            plot_cool_contours(
                cat[filter_selections[0]] - cat[filter_selections[1]],
                cat[filters[j]],
                kvals,
                ax,
                contour_color=next(colors),
                contour_lw=1,
            )

        ax.set_xlim(kvals[0], kvals[1])
        ax.set_ylim(kvals[2], kvals[3])

    plt.subplots_adjust(hspace=0.5)
    if savefig:
        plt.savefig(ast_file.replace(".txt", "_diagnostic.png"), format="png")


def plot_cool_contours(
    xx, yy, kvals, ax, contour_color="black", contour_lw=2, contour_alpha=1
):

    """
    Plot 1, 2, and 3 sigma contours from x,y data using simple Gaussian KDE.

    Parameters
    ----------
    xx : list
        x-axis data

    yy : list
        x-axis data

    kvals : list
        list specifying the extent of the contour range, in the form of [x_min, x_max, y_min, y_max]

    ax : axis instance
        axis on which to plot the contours

    contour_color : string (default 'black')
        color for contours

    contour_lw : float
        line width for the contours

    contour_alpha : float
        alpha value for contours

    """

    # evaluate on a regular grid
    data = np.vstack([xx, yy])
    kde = gaussian_kde(data)

    xgrid = np.linspace(kvals[0], kvals[1], 50)
    ygrid = np.linspace(kvals[2], kvals[3], 50)

    Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
    z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))
    z = z.reshape(Xgrid.shape)
    z = z / z.sum()

    n = 1000
    t = np.linspace(0, z.max(), n)
    integral = ((z >= t[:, None, None]) * z).sum(axis=(1, 2))

    f = interpolate.interp1d(integral, t)
    t_contours = f(np.array([0.98, 0.95, 0.68]))

    ax.contour(
        z,
        t_contours,
        extent=kvals,
        colors=contour_color,
        zorder=1,
        linewidths=contour_lw,
        alpha=contour_alpha,
    )


if __name__ == "__main__":  # pragma: no cover

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "beast_settings_info", type=str, help="file name with beast settings",
    )
    parser.add_argument("ast_file", type=str, help="name of AST input file")
    parser.add_argument(
        "map_file", type=str, help="background or source density map file"
    )
    parser.add_argument(
        "source_density_image",
        type=str,
        help="name of the source density image FITS file",
    )
    parser.add_argument(
        "reference_image", type=str, help="name of the reference image FITS file"
    )
    parser.add_argument(
        "filter_selections", type=list, help="filter names for CMD reference color"
    )
    parser.add_argument(
        "savefig", type=bool, help="to save or not to save figure to file"
    )

    args = parser.parse_args()

    plot_ast_histogram(
        args.beast_settings_info,
        args.ast_file,
        args.map_file,
        args.source_density_image,
        args.reference_image,
        args.n_per_file,
        args.min_n_subfile,
        args.sort_col,
    )
