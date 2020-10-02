#!/usr/bin/env python
#
# condense separate BEAST stats files into a single stats file
import glob

import argparse

from astropy.table import Table, Column, vstack
from astropy.io import fits


def merge_stats_files(stats_files, out_stats_filebase, reorder_tag_list=None):
    """
    Merge stats files.

    Parameters
    ----------
    stats_files : list of strings
        names of the stats files to be merged

    out_stats_filebase : string
        base name of the output file ('_stats.fits' will be appended)

    reorder_tag_list : list of strings (default=None)
        If set, these tags will be used in the output stats file.  If not
        set, the tag names will be derived from the stats file names.
    """

    # grab filter table if it exists
    try:
        filters_tab = Table.read(stats_files[0], hdu=2)
    except ValueError:
        filters_tab = None

    # loop through the stats files, building up the output table
    cats_list = []
    for i, cur_stat in enumerate(stats_files):

        # read in current catalog
        cur_cat = Table.read(cur_stat, hdu=1)

        if reorder_tag_list is None:
            # get the source density and subregion name
            #  in other words, the reordering tag
            #  bpos is the location after the 2nd underscore of pix coords
            #  epos is before the _stats.fits ending string
            bpos = cur_stat.find("_sd") + 1
            epos = cur_stat.find("_stats")
            reorder_tag = cur_stat[bpos:epos]
        else:
            reorder_tag = reorder_tag_list[i]

        # add the reorder tag to each entry in the current catalog
        n_entries = len(cur_cat)
        cur_cat.add_column(Column([reorder_tag] * n_entries, name="reorder_tag"))

        # append to list
        cats_list.append(cur_cat)

    # concatenate all the small catalogs together
    full_cat = vstack(cats_list)

    # output the full pixel catalog
    ohdu = fits.HDUList()
    ohdu.append(fits.table_to_hdu(full_cat))
    if filters_tab is not None:
        ohdu.append(fits.table_to_hdu(filters_tab))
    ohdu.writeto(out_stats_filebase + "_stats.fits", overwrite=True)

    # return the number of sources in the catalog for later use
    return len(full_cat)


if __name__ == "__main__":  # pragma: no cover

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filebase", help="filebase to use (e.g., xxx_*_stats.fits)")
    args = parser.parse_args()

    # get the files to merge
    stats_files = glob.glob(args.filebase + "*_stats.fits")

    # do the merge
    n_objs = merge_stats_files(stats_files, args.filebase)
