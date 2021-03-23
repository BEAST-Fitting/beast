import numpy as np
import argparse
from matplotlib.path import Path
from scipy.spatial import ConvexHull

from astropy.table import Table


def cut_catalogs(
    input_phot_file,
    output_phot_file,
    input_ast_file=None,
    output_ast_file=None,
    partial_overlap=False,
    flagged=False,
    flag_filter=None,
    region_file=False,
    no_write=False,
):
    """
    Remove sources from the input photometry catalog that are
    - in regions without full imaging coverage
    OR
    - flagged as bad in flag_filter

    Optionally also input an AST file to similarly cut.  If partial_overlap is
    True, the overlap boundary found for the photometry catalog will be used to
    trim the AST catalog.

    Parameters
    ----------
    input_phot_file : string
        file name of the photometry catalog

    output_phot_file : string
        file name for the output photometry catalog

    input_ast_file : string (default=None)
        file name of the AST catalog

    output_ast_file : string (default=None)
        file name for the output AST catalog

    partial_overlap : boolean (default=False)
        if True, remove sources in regions without full imaging coverage.  This
        is done by finding sources with RATE=0 in any filter.  If a source has
        RATE=0 in all filters (which can only happen for ASTs), that means the
        source was not recovered, and since we need to keep that information,
        the source will not be removed.

    flagged : boolean (default=False)
        if True, remove sources with flag=99 in flag_filter.  Note that this
        will only be applied to sources with RATE>0 in flag_filter.

    flag_filter : string or list of strings (default=None)
        if flagged is True, set this to the filter(s) to use

    region_file : boolean (default=False)
        if True, create ds9 region file(s) where good sources are green and
        removed sources are magenta

    no_write : boolean (default=False)
        if True, don't write out the catalog file(s).  Instead, just return them.

    """

    # make sure something is chosen
    if (partial_overlap is False) and (flagged is False):
        raise ValueError("must choose a criteria to cut catalogs")

    # run the cutting for the photometry file
    phot_cat, phot_good_stars = make_cuts(
        input_phot_file,
        partial_overlap=partial_overlap,
        flagged=flagged,
        flag_filter=flag_filter,
    )
    print(
        "removing {0} stars from {1}".format(
            int(len(phot_cat) - np.sum(phot_good_stars)), input_phot_file
        )
    )
    new_phot_cat = phot_cat[phot_good_stars == 1]

    # if chosen, run the cutting for the AST file
    if input_ast_file is not None:
        ast_cat, temp_good_stars = make_cuts(
            input_ast_file,
            partial_overlap=partial_overlap,
            flagged=flagged,
            flag_filter=flag_filter,
        )

        # do convex hull
        ra_col = [x for x in phot_cat.colnames if "RA" in x.upper()][0]
        dec_col = [x for x in phot_cat.colnames if "DEC" in x.upper()][0]
        phot_cat_boundary = convexhull_path(new_phot_cat[ra_col], new_phot_cat[dec_col])
        # check if points are inside
        ra_col = [x for x in ast_cat.colnames if "RA" in x.upper()][0]
        dec_col = [x for x in ast_cat.colnames if "DEC" in x.upper()][0]
        within_bounds = phot_cat_boundary.contains_points(
            np.array([ast_cat[ra_col], ast_cat[dec_col]]).T
        )  # N,2 array of AST RA and Dec positions

        # get final list of good stars
        ast_good_stars = (temp_good_stars == 1) & (within_bounds)
        print(
            "removing {0} stars from {1}".format(
                int(len(ast_cat) - np.sum(ast_good_stars)), input_ast_file
            )
        )
        new_ast_cat = ast_cat[ast_good_stars]

    # write out the sources as a ds9 region file
    if region_file is True:
        write_ds9(phot_cat, phot_good_stars, input_phot_file + ".reg")

        if input_ast_file is not None:
            write_ds9(ast_cat, ast_good_stars, input_ast_file + ".reg")

    # either save it or return it
    # - save it
    if not no_write:
        new_phot_cat.write(output_phot_file, format="fits", overwrite=True)
        if input_ast_file is not None:
            new_ast_cat.write(output_ast_file, format="fits", overwrite=True)
    # - return it
    else:
        return new_phot_cat, phot_good_stars


def make_cuts(cat_file, partial_overlap=False, flagged=False, flag_filter=None):
    """
    Wrapper to do overlap and flag cuts

    Parameters
    ----------
    cat_file : string
        name of the catalog file

    partial_overlap : boolean (default=False)
        see above

    flagged : boolean (default=False)
        see above

    flag_filter : string (default=None)
        see above

    Returns
    -------
    cat : astropy table
        the astropy table of the input catalog

    good_stars : array
        array in which 1=good star and 0=bad star
    """

    # read in the catalog
    cat = Table.read(cat_file)
    filters = [c[0:-5] for c in cat.colnames if "RATE" in c and "RATERR" not in c]
    n_stars = len(cat)

    # array to keep track of which ones are good
    # (1=good, 0=bad)
    good_stars = np.ones(n_stars)

    # partial overlap
    if partial_overlap is True:
        # number of RATE=0 for each source
        n_zero_flux = np.sum([cat[filt + "_RATE"] == 0 for filt in filters], axis=0)
        # remove sources with more than 0 and less than n_filter
        good_stars[(n_zero_flux > 0) & (n_zero_flux < len(filters))] = 0

    # flagged sources
    if flagged is True:
        for fl in np.atleast_1d(flag_filter):
            good_stars[(cat[fl + "_FLAG"] >= 99) & (cat[fl + "_RATE"] > 0)] = 0

    # return results
    return cat, good_stars


def write_ds9(cat, good_stars, ds9_file_name):
    """
    Make a ds9 region file
    """

    ra_col = [x for x in cat.colnames if "RA" in x.upper()][0]
    dec_col = [x for x in cat.colnames if "DEC" in x.upper()][0]

    with open(ds9_file_name, "w") as ds9_file:

        # header
        ds9_file.write("# Region file format: DS9 version 4.1\n")
        ds9_file.write(
            'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
        )
        ds9_file.write("fk5\n")

        # stars
        for i in range(len(cat)):
            if good_stars[i] == 0:
                ds9_file.write(
                    "circle("
                    + str(cat[ra_col][i])
                    + ","
                    + str(cat[dec_col][i])
                    + ',0.1") # color=magenta\n'
                )
            if good_stars[i] == 1:
                ds9_file.write(
                    "circle("
                    + str(cat[ra_col][i])
                    + ","
                    + str(cat[dec_col][i])
                    + ',0.1")\n'
                )


def convexhull_path(x_coord, y_coord):
    """
    Find the convex hull for the given coordinates and make a Path object
    from it.

    Parameters
    ----------
    x_coord, y_coord: arrays
        Arrays of coordinates (can be x/y or ra/dec)

    Returns
    -------
    matplotlib Path object

    """

    coords = np.array(
        [x_coord, y_coord]
    ).T  # there's a weird astropy datatype issue that requires numpy coercion
    hull = ConvexHull(coords)
    bounds_x, bounds_y = coords[hull.vertices, 0], coords[hull.vertices, 1]
    path_object = Path(np.array([bounds_x, bounds_y]).T)

    return path_object


if __name__ == "__main__":  # pragma: no cover
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_phot_file",
        type=str,
        help="file name of the input photometry catalog",
    )
    parser.add_argument(
        "output_phot_file",
        type=str,
        help="file name for the output photometry catalog",
    )
    parser.add_argument(
        "--input_ast_file",
        type=str,
        default=None,
        help="file name for the input AST catalog",
    )
    parser.add_argument(
        "--output_ast_file",
        type=str,
        default=None,
        help="file name for the output AST catalog",
    )
    parser.add_argument(
        "--partial_overlap",
        default=False,
        action="store_true",
        help="if True, remove sources in regions without full imaging coverage",
    )
    parser.add_argument(
        "--flagged",
        default=False,
        action="store_true",
        help="if True, remove sources with flag=99 in flag_filter",
    )
    parser.add_argument(
        "--flag_filter",
        default=None,
        type=str,
        help="if flagged is True, set this to the filter(s) to use",
    )
    parser.add_argument(
        "--region_file",
        default=False,
        action="store_true",
        help="""if True, create a ds9 region file where good sources are green
        and removed sources are magenta""",
    )
    parser.add_argument(
        "--no_write",
        default=False,
        action="store_true",
        help="if True, don't write out the catalog file; instead, just return it",
    )

    args = parser.parse_args()

    cut_catalogs(
        input_phot_file=args.input_phot_file,
        output_phot_file=args.output_phot_file,
        input_ast_file=args.input_ast_file,
        output_ast_file=args.output_ast_file,
        partial_overlap=args.partial_overlap,
        flagged=args.flagged,
        flag_filter=args.flag_filter,
        region_file=args.region_file,
        no_write=args.no_write,
    )
