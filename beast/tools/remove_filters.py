#!/usr/bin/env python
#
# remove filters from photometry catalogs, physicsgrid, and observationgrid
#   used to modify simulated data to make plots for proposals
import argparse

import numpy as np
from astropy.table import Table
import tables

from beast.physicsmodel.grid import FileSEDGrid, SpectralGrid
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel


def remove_filters_from_files(
    catfile,
    physgrid=None,
    obsgrid=None,
    outbase=None,
    physgrid_outfile=None,
    rm_filters=None,
    beast_filt=None,
):
    """
    Remove filters from catalog, physics grid, and/or obsmodel grid.  This has
    two primary use cases:

    1. When making simulated observations, you want to test how your fit quality
       changes with different combinations of filters.  In that case, put in
       files for both `physgrid` and `obsgrid`.  Set `rm_filters` to the
       filter(s) you wish to remove, and they will be removed both from those
       and from the catalog file.  The three new files will be output with the
       name prefix set in `outbase`.

    2. When running the BEAST, you have a master physics model grid with all
       filters present in the survey, but some fields don't have observations in
       all of those filters.  In that case, put the master grid in `physgrid`
       and set `rm_filters` to None.  The catalog will be used to determine the
       filters to remove (if any).  `obsgrid` should be left as None, because in
       this use case, the obsmodel grid has not yet been generated.  The output
       physics model grid will be named using the filename in `physgrid_outfile`
       (if given) or with the prefix in `outbase`.


    Parameters
    ----------
    catfile : string
        file name of photometry catalog

    physgrid : string (default=None)
        If set, remove filters from this physics model grid

    obsgrid : string (default=None)
        If set, remove filters from this obsmodel grid

    outbase : string (default=None)
        Path+file to prepend to all output file names.  Useful for case 1 above.

    physgrid_outfile : string (default=None)
        Path+name of the output physics model grid.  Useful for case 2 above.

    rm_filters : string or list of strings (default=None)
        If set, these are the filters to remove from all of the files.  If not
        set, only the filters present in catfile will be retained in physgrid
        and/or obsgrid.

    beast_filt : list of strings
        Sometimes there is ambiguity in the filter name (e.g., the grid has
        both HST_ACS_WFC_F475W and HST_WFC3_F475W, and the filter name is
        F475W).  Set this to the BEAST filter name to resolve any
        ambiguities.  For example, ['HST_WFC3_F475W', 'HST_WFC3_F814W'] ensures
        that these are the names used for F475W and F814W.

    """

    # read in the photometry catalog
    cat = Table.read(catfile)

    # if rm_filters set, remove the requested filters from the catalog
    if rm_filters is not None:
        for cfilter in np.atleast_1d(rm_filters):
            colname = "{}_rate".format(cfilter)
            if colname.upper() in cat.colnames:
                cat.remove_column(colname.upper())
            elif colname.lower() in cat.colnames:
                cat.remove_column(colname.lower())
            else:
                print("{} not in catalog file".format(colname))
        cat.write("{}_cat.fits".format(outbase), overwrite=True)

    # if rm_filters not set, extract the filter names that are present
    if rm_filters is None:
        cat_filters = [f[:-5].upper() for f in cat.colnames if f[-4:].lower() == "rate"]

    # if beast_filt is set, make a list of the short versions
    if beast_filt is not None:
        beast_filt_short = [(f.split("_"))[-1].upper() for f in beast_filt]

    # if physgrid set, process the SED grid
    if physgrid is not None:

        # read in the sed grid
        g0 = FileSEDGrid(physgrid, backend="cache")

        # extract info
        filters = g0.header["filters"].split(" ")
        shortfilters = [(cfilter.split("_"))[-1].upper() for cfilter in filters]
        rindxs = []
        rgridcols = []

        # loop through filters and determine what needs deleting
        for csfilter, cfilter in zip(shortfilters, filters):

            # --------------------------
            # if the user chose the filters to remove
            if rm_filters is not None:

                # if the current filter is in the list of filters to remove
                if csfilter in np.atleast_1d(rm_filters):

                    # if there's a list of BEAST instrument+filter references
                    if beast_filt is not None:

                        # if the current filter is in the list of BEAST references
                        if csfilter in beast_filt_short:

                            # if it's the same instrument, delete it
                            # (if it's not the same instrument, keep it)
                            if beast_filt[beast_filt_short.index(csfilter)] == cfilter:
                                rindxs.append(filters.index(cfilter))
                                for grid_col in g0.grid.colnames:
                                    if cfilter in grid_col:
                                        rgridcols.append(grid_col)

                        # if the current filter isn't in the BEAST ref list, delete it
                        else:
                            rindxs.append(filters.index(cfilter))
                            for grid_col in g0.grid.colnames:
                                if cfilter in grid_col:
                                    rgridcols.append(grid_col)

                    # if there isn't a list of BEAST refs, delete it
                    else:
                        rindxs.append(filters.index(cfilter))
                        for grid_col in g0.grid.colnames:
                            if cfilter in grid_col:
                                rgridcols.append(grid_col)

            # --------------------------
            # if the removed filters are determined from the catalog file
            if rm_filters is None:

                # if the current filter is present in the catalog filters
                if csfilter in cat_filters:

                    # if there's a list of BEAST instrument+filter references
                    # (if there isn't a list of BEAST refs, keep it)
                    if beast_filt is not None:

                        # if the current filter is in the list of BEAST references
                        # (if the current filter isn't in the BEAST ref list, keep it)
                        if csfilter in beast_filt_short:

                            # if it's not the same instrument, delete it
                            # (if it's the same instrument, keep it)
                            if beast_filt[beast_filt_short.index(csfilter)] != cfilter:
                                rindxs.append(filters.index(cfilter))
                                for grid_col in g0.grid.colnames:
                                    if cfilter in grid_col:
                                        rgridcols.append(grid_col)

                # if the current filter isn't in the catalog filters, delete it
                else:
                    rindxs.append(filters.index(cfilter))
                    for grid_col in g0.grid.colnames:
                        if cfilter in grid_col:
                            rgridcols.append(grid_col)

        # delete column(s)
        nseds = np.delete(g0.seds, rindxs, 1)
        nlamb = np.delete(g0.lamb, rindxs, 0)
        nfilters = np.delete(filters, rindxs, 0)
        for rcol in rgridcols:
            g0.grid.delCol(rcol)

        print("orig filters: {}".format(" ".join(filters)))
        print(" new filters: {}".format(" ".join(nfilters)))

        # save the modified grid
        g = SpectralGrid(np.array(nlamb), seds=nseds, grid=g0.grid, backend="memory")
        g.grid.header["filters"] = " ".join(nfilters)
        if physgrid_outfile is not None:
            g.writeHDF(physgrid_outfile)
        elif outbase is not None:
            g.writeHDF("{}_seds.grid.hd5".format(outbase))
        else:
            raise ValueError("Need to set either outbase or physgrid_outfile")

    # if obsgrid set, process the observation model
    if obsgrid is not None:
        obsgrid = noisemodel.get_noisemodelcat(obsgrid)
        with tables.open_file("{}_noisemodel.grid.hd5".format(outbase), "w") as outfile:
            outfile.create_array(
                outfile.root, "bias", np.delete(obsgrid["bias"], rindxs, 1)
            )
            outfile.create_array(
                outfile.root, "error", np.delete(obsgrid["error"], rindxs, 1)
            )
            outfile.create_array(
                outfile.root,
                "completeness",
                np.delete(obsgrid["completeness"], rindxs, 1),
            )


if __name__ == "__main__":  # pragma: no cover

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("catfile", help="filename of photometry catalog")
    parser.add_argument(
        "--physgrid",
        type=str,
        default=None,
        help="If set, remove filters from this physics model grid file",
    )
    parser.add_argument(
        "--obsgrid",
        type=str,
        default=None,
        help="If set, remove filters from this observation/noisemodel grid file",
    )
    parser.add_argument(
        "--outbase",
        type=str,
        default=None,
        help="Path+file to prepend to all output file names",
    )
    parser.add_argument(
        "--physgrid_outfile",
        type=str,
        default=None,
        help="""Path+name of the output physics model grid. Takes precendence
        over the default file name constructed from outbase.""",
    )
    parser.add_argument(
        "--rm_filters",
        type=str,
        nargs="*",
        default=None,
        help="""If set, these are the filters to remove from all of the files.
        If not set, only the filters present in catfile will be retained in
        physgrid and/or obsgrid.""",
    )
    args = parser.parse_args()

    # do the filter removal
    remove_filters_from_files(
        args.catfile,
        physgrid=args.physgrid,
        obsgrid=args.obsgrid,
        outbase=args.outbase,
        physgrid_outfile=args.physgrid_outfile,
        rm_filters=args.rm_filters,
    )
