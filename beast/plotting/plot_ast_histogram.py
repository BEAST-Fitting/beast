import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.table import Table

from beast.physicsmodel.grid import SEDGrid
from beast.observationmodel.vega import Vega

__all__ = ["plot_ast_histogram"]


def plot_ast_histogram(ast_file, sed_grid_file=None):
    """
    Make a histogram of the AST fluxes.  If an SED grid is given, also plot
    a comparison histogram of the SED fluxes.

    The histogram bins are set by the bins originally used to create the ASTs
    (using the flux bin method), which are saved in
    ast_file.replace('inputAST','ASTfluxbins')
    and are automatically read in.

    Output plot is saved in the same location/name as ast_file, but with a .png
    instead of .txt.

    Parameters
    ----------
    ast_file : string
        name of AST input file

    sed_grid_file : string (default=None)
        name of SED grid file
    """
    # read in AST info
    ast_table = Table.read(ast_file, format="ascii")
    ast_fluxbins = Table.read(
        ast_file.split("inputAST")[0] + "ASTfluxbins.txt", format="ascii"
    )

    # get filter names (and it's even in order!)
    filter_cols = [col for col in ast_table.colnames if "_" in col]
    filter_list = [col[-5:] for col in filter_cols]
    n_filter = len(filter_list)

    # if chosen, read in model grid
    if sed_grid_file is not None:
        modelsedgrid = SEDGrid(sed_grid_file)
        with Vega() as v:
            _, vega_flux, _ = v.getFlux(filter_cols)
        sedsMags = -2.5 * np.log10(modelsedgrid.seds[:] / vega_flux)

    # make a histogram for each filter
    fig = plt.figure(figsize=(7, 4 * n_filter))

    for f, filt in enumerate(filter_list):

        # names of table columns with bin values
        min_bin_col = [b for b in ast_fluxbins.colnames if ("mins" in b and filt in b)][
            0
        ]
        max_bin_col = [b for b in ast_fluxbins.colnames if ("maxs" in b and filt in b)][
            0
        ]
        # use those to make a bin list
        bin_list = np.append(ast_fluxbins[min_bin_col], ast_fluxbins[max_bin_col][-1])

        # make histograms
        ax = plt.subplot(n_filter, 1, f + 1)

        ast_col = [b for b in ast_table.colnames if filt in b][0]

        plt.hist(
            ast_table[ast_col],
            bins=bin_list,
            density=True,
            facecolor="black",
            edgecolor="none",
            alpha=0.3,
            label="ASTs",
        )

        if sed_grid_file is not None:
            plt.hist(
                sedsMags[:, f],
                bins=bin_list,
                density=True,
                histtype="step",
                facecolor="none",
                edgecolor="black",
                label="Model grid",
            )

        # labels
        ax.tick_params(axis="both", which="major", labelsize=13)
        ax.set_xlim(ax.get_xlim()[::-1])
        plt.xlabel(filt + " (Vega mag)", fontsize=14)
        plt.ylabel("Normalized Histogram", fontsize=14)

        # add legend in first plot
        if f == 0:
            ax.legend(fontsize=14)

    plt.tight_layout()

    fig.savefig(ast_file.replace(".txt", ".png"))
    plt.close(fig)


if __name__ == "__main__":  # pragma: no cover

    parser = argparse.ArgumentParser()
    parser.add_argument("ast_file", type=str, help="name of AST input file")
    parser.add_argument(
        "--sed_grid_file", type=str, default=None, help="name of SED grid file"
    )

    args = parser.parse_args()

    plot_ast_histogram(args.ast_file, sed_grid_file=args.sed_grid_file)
