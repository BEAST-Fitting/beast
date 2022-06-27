#!/usr/bin/env python3
"""
Split a catalog and AST results into grid subbins to improve fitting efficiency for large physics grids.
This is almost identical to split_catalog_using_map except this function does not require spatial
information (e.g. RA, DEC) to sort by source density or background emission density, making it
optimal for simulated catalogs.
"""
import argparse
import numpy as np
from astropy.table import Table


def split_main(
    beast_settings_info,
    catfile,
    astfile,
    n_per_file=1000,
    min_n_subfile=None,
    sort_col="F475W_RATE",
):
    """
    Making the physics model grid takes a while for production runs.  This
    creates scripts to run each subgrid as a separate job.

    Parameters
    ----------
    beast_settings_info : string or instance
        if string: file name with beast settings
        if class: beast.tools.beast_settings.beast_settings instance

    catfile : string
        name of the photometry catalog file

    astfile : string
        name of the ast catalog file

    n_per_file : int or None (default=1000)
        If set, divide the split catalog into sub-catalogs with length
        n_per_file.  Good for photometry, not useful for ASTs.

    min_n_subfile : int or None (default=None)
        If set, each bin in the photometry catalog will be split into at least
        this many subfiles. Useful if a bin has fewer than n_per_file stars but
        you still want flux-sorted subfiles (which means more trimming and
        faster fitting).

    sort_col : string (default="F475W_RATE")
        If n_per_file or min_n_subfile is set, the catalog will be sorted by this
        column before splitting into sub-catalogs.


    """

    print("Splitting catalog")
    split_simulated_catalog(
        catfile,
        n_per_file=n_per_file,
        min_n_subfile=min_n_subfile,
        sort_col=sort_col,
    )


def split_simulated_catalog(
    catfile,
    n_per_file=1000,
    min_n_subfile=None,
    sort_col="F475W_RATE",
):
    """
    Code to do the splitting of a catalog

    Parameters
    ----------
    catfile : string
        name of the photometry catalog file

    n_per_file : int or None (default=6250)
        If set, divide the split catalog into sub-catalogs with length
        n_per_file.  Good for photometry, not useful for ASTs.

    min_n_subfile : int or None (default=None)
        If set, each bin in the photometry catalog will be split into at least
        this many subfiles. Useful if a bin has fewer than n_per_file stars but
        you still want flux-sorted subfiles (which means more trimming and
        faster fitting).

    sort_col : string (default="F475W_RATE")
        If n_per_file or min_n_subfile is set, the catalog will be sorted by this
        column before splitting into sub-catalogs.


    """
    cat = Table.read(catfile)

    # write out sub-files, if chosen
    if (n_per_file is not None) or (min_n_subfile is not None):

        # calculate number of subfiles and number of stars per file
        # - only n_per_file set
        if (n_per_file is not None) and (min_n_subfile is None):
            tot_subfiles = int(np.ceil(len(cat) / n_per_file))
            curr_n_per_file = n_per_file
        # - only min_n_subfile set
        if (n_per_file is None) and (min_n_subfile is not None):
            tot_subfiles = min_n_subfile
            curr_n_per_file = int(np.ceil(len(cat) / tot_subfiles))
        # - both are set: make sure the largest number of subfiles is used
        if (n_per_file is not None) and (min_n_subfile is not None):
            temp_tot_subfiles = int(np.ceil(len(cat) / n_per_file))
            # n_per_file makes at least min_n_subfile -> use value from n_per_file
            if min_n_subfile <= temp_tot_subfiles:
                tot_subfiles = temp_tot_subfiles
                curr_n_per_file = n_per_file
            # n_per_file doesn't make enough subfiles -> use min_n_subfile
            else:
                tot_subfiles = min_n_subfile
                curr_n_per_file = int(
                    np.ceil(len(cat) / tot_subfiles)
                )

        print(
            "dividing into "
            + str(tot_subfiles)
            + " subfiles for later fitting speed"
        )

        # Sort the stars
        sort_indxs = np.argsort(cat[sort_col])

        for i in range(tot_subfiles):
            min_k = i * curr_n_per_file
            if i < tot_subfiles:
                max_k = (i + 1) * curr_n_per_file
            else:
                max_k = len(cat)

            cat[sort_indxs[min_k:max_k]].write(
                catfile.replace(".fits", "_bin0_sub{0}.fits".format(i)),
                overwrite=True,
            )


if __name__ == "__main__":  # pragma: no cover
    parser = argparse.ArgumentParser()
    parser.add_argument("catfile", type=str, help="catalog FITS file")
    parser.add_argument("astfile", type=str, help="ast results fits file")
    parser.add_argument(
        "--n_per_file", type=int, default=6250, help="Number of sources per subfile"
    )

    parser.add_argument(
        "--min_n_subfile",
        type=int,
        default=None,
        help="""Set the minimum number of subfiles to use per bin (relevant if a
        bin has fewer than n_per_file but you still want flux-sorted subfiles)""",
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
        args.n_per_file,
        args.min_n_subfile,
        args.sort_col,
    )
