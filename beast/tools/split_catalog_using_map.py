#!/usr/bin/env python3
"""
Split a catalog and a set of AST results, by source or background
density bin. Uses one of the maps created by
'create_background_density_map'. From the split AST catalog, individual
noise models for different regions can be made, which can then be used
to fit the stars of the observed catalog which also fall in those
regions.

"""
import argparse
import numpy as np
from astropy.table import Table
from beast.tools.density_map import BinnedDensityMap


def main():  # pragma: no cover
    parser = argparse.ArgumentParser()
    parser.add_argument("catfile", type=str, help="catalog FITS file")
    parser.add_argument("astfile", type=str, help="ast results fits file")
    parser.add_argument(
        "mapfile", type=str, help="background or source density map file"
    )

    nbin_or_binwidth = parser.add_mutually_exclusive_group()
    nbin_or_binwidth.add_argument(
        "--n",
        type=int,
        help="Number of regions to split the catalog in, the number of density bins",
    )
    nbin_or_binwidth.add_argument(
        "--bin_width",
        type=float,
        default=None,
        help="Width of the density bins for splitting catalog",
    )

    parser.add_argument(
        "--n_per_file", type=int, default=6250, help="Number of sources per subfile"
    )

    parser.add_argument(
        "--sort_col",
        type=str,
        default="F475W_RATE",
        help="If n_per_file set, sort catalog by this column before splitting",
    )

    args = parser.parse_args()

    split_main(
        args.catfile,
        args.astfile,
        args.mapfile,
        args.n,
        args.bin_width,
        args.n_per_file,
        args.sort_col,
    )


def split_main(
    catfile,
    astfile,
    mapfile,
    n_bin=None,
    bin_width=None,
    n_per_file=6250,
    sort_col="F475W_RATE",
):

    # Create a binned density map, so both the observed and the ast
    # catalog can be split using a consistent grouping (= binning) of
    # the tiles
    bdm = BinnedDensityMap.create(mapfile, N_bins=n_bin, bin_width=bin_width)

    print("Splitting catalog")
    split_catalog_using_map(catfile, bdm, n_per_file=n_per_file, sort_col=sort_col)
    print("")
    print("Splitting ASTs")
    split_catalog_using_map(
        astfile, bdm, ra_colname="RA_J2000", dec_colname="DEC_J2000", n_per_file=None
    )


def split_catalog_using_map(
    catfile,
    binned_density_map,
    ra_colname="RA",
    dec_colname="DEC",
    n_per_file=6250,
    sort_col="F475W_RATE",
):
    """
    Code to do the splitting of a catalog

    Parameters
    ----------
    catfile : string
        name of the photometry catalog file

    binned_density_map : BinnedDensityMap object
        the binned density map for the field

    ra_colname, dec_colname : string
        labels for the RA and DEC columns

    n_per_file : int or None (default=6250)
        If set, divide the split catalog into sub-catalogs with length
        n_per_file.  Good for photometry, not useful for ASTs.

    sort_col : string (default="F475W_RATE")
        If n_per_file is set, the catalog will be sorted by this column
        before splitting into sub-catalogs.


    """
    cat = Table.read(catfile)

    ras = cat[ra_colname]
    decs = cat[dec_colname]

    bin_foreach_source = np.zeros(len(cat), dtype=int)
    for i in range(len(cat)):
        bin_foreach_source[i] = binned_density_map.bin_for_position(ras[i], decs[i])

    binnrs = np.unique(bin_foreach_source)

    for b in binnrs:
        # write out file for this bin
        sources_for_bin = np.where(bin_foreach_source == b)
        print("bin {0}: {1} sources".format(b, len(sources_for_bin[0])))
        subcat = cat[sources_for_bin]
        subcat.write(catfile.replace(".fits", "_bin{}.fits".format(b)), overwrite=True)

        # write out sub-files, if chosen
        if n_per_file is not None:

            tot_subfiles = int(np.ceil(len(sources_for_bin[0]) / n_per_file))
            print(
                "dividing into "
                + str(tot_subfiles)
                + " subfiles for later fitting speed"
            )

            # Sort the stars
            sort_indxs = np.argsort(subcat[sort_col])

            for i in range(tot_subfiles):
                min_k = i * n_per_file
                if i < tot_subfiles:
                    max_k = (i + 1) * n_per_file
                else:
                    max_k = N

                subcat[sort_indxs[min_k:max_k]].write(
                    catfile.replace(".fits", "_bin{0}_sub{1}.fits".format(b, i)),
                    overwrite=True,
                )


if __name__ == "__main__":  # pragma: no cover
    main()
