import numpy as np
import argparse

from astropy.table import Table


def cut_catalogs(
    input_file,
    output_file,
    partial_overlap=False,
    flagged=False,
    flag_filter=None,
    region_file=False,
    no_write=False,
):
    """
    Remove sources from the input catalog that are
    - in regions without full imaging coverage
    OR
    - flagged as bad in flag_filter

    The input catalog can either be photometry or ASTs

    Parameters
    ----------
    input_file : string
        file name of the input catalog (photometry or AST), where data is in
        extension 1 of the fits file

    output_file : string
        file name for the output catalog

    partial_overlap : boolean (default=False)
        if True, remove sources in regions without full imaging coverage

    flagged : boolean (default=False)
        if True, remove sources with flag=99 in flag_filter

    flag_filter : string or list of strings (default=None)
        if flagged is True, set this to the filter(s) to use

    region_file : boolean (default=False)
        if True, create a ds9 region file where good sources are green and
        removed sources are magenta

    no_write : boolean (default=False)
        if True, don't write out the catalog file.  Instead, just return it.

    """

    # make sure something is chosen
    if (partial_overlap is False) and (flagged is False):
        print("must choose a criteria to cut catalogs")
        return

    # read in the catalog
    cat = Table.read(input_file)
    filters = [c[0:-5] for c in cat.colnames if "RATE" in c and "RATERR" not in c]
    n_stars = len(cat)
    if "RA" in cat.colnames:
        ra_col = "RA"
        dec_col = "DEC"
    else:
        ra_col = "RA_J2000"
        dec_col = "DEC_J2000"

    # array to keep track of which ones are good
    # (1=good, 0=bad)
    good_stars = np.ones(n_stars)

    # partial overlap
    if partial_overlap is True:
        for filt in filters:
            good_stars[cat[filt + "_RATE"] == 0] = 0

    # flagged sources
    if flagged is True:
        if type(flag_filter) == str:
            good_stars[cat[flag_filter + "_FLAG"] >= 99] = 0
        else:
            for fl in flag_filter:
                good_stars[cat[fl + "_FLAG"] >= 99] = 0

    # write out the sources as a ds9 region file
    if region_file is True:
        with open(input_file + ".reg", "w") as ds9_file:

            # header
            ds9_file.write("# Region file format: DS9 version 4.1\n")
            ds9_file.write(
                'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
            )
            ds9_file.write("fk5\n")

            # stars
            for i in range(n_stars):
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

    # make a new file with the bad stars removed
    new_cat = cat[good_stars == 1]

    # either save it or return in
    # - save it
    if no_write == False:
        new_cat.write(output_file, format="fits", overwrite=True)
    # - return in
    else:
        return new_cat, good_stars


if __name__ == "__main__":
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file",
        type=str,
        help="file name of the input catalog (photometry or AST)",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="file name for the output catalog",
    )
    parser.add_argument(
        "--partial_overlap",
        default=False,
        action='store_true'
        help="if True, remove sources in regions without full imaging coverage",
    )
    parser.add_argument(
        "--flagged",
        default=False,
        action='store_true'
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
        action='store_true'
        help="""if True, create a ds9 region file where good sources are green
        and removed sources are magenta""",
    )
    parser.add_argument(
        "--no_write",
        default=False,
        action='store_true'
        help="if True, don't write out the catalog file; instead, just return it",
    )

    args = parser.parse_args()

    cut_catalogs(
        input_file=args.input_file,
        output_file=args.output_file,
        partial_overlap=args.partial_overlap,
        flagged=args.flagged,
        flag_filter=args.flag_filter,
        region_file=args.region_file,
        no_write=args.no_write,
    )

    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()
