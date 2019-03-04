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
from density_map import BinnedDensityMap


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('catfile', type=str, help='catalog FITS file')
    parser.add_argument('astfile', type=str, help='ast results fits file')
    parser.add_argument('mapfile', type=str,
                        help='background or source density map file')
    parser.add_argument('n', type=int,
                        help='Number of regions to split the catalog in, the number of density bins')
    args = parser.parse_args()

    # Create a binned density map, so both the observed and the ast
    # catalog can be split using a consistent grouping (= binning) of
    # the tiles
    bdm = BinnedDensityMap.create(args.mapfile, args.n)

    split_catalog_using_map(args.catfile, bdm, args.n)
    split_catalog_using_map(args.astfile, bdm, args.n,
                            ra_colname='RA_J2000', dec_colname='DEC_J2000')

def split_catalog_using_map(catfile, binned_density_map, n, ra_colname='RA', dec_colname='DEC'):
    cat = Table.read(catfile)

    ras = cat[ra_colname]
    decs = cat[dec_colname]

    bin_foreach_source = np.zeros(len(cat), dtype=int)
    for i in range(len(cat)):
        bin_foreach_source[i] = binned_density_map.bin_for_position(ras[i], decs[i])

    binnrs = np.unique(bin_foreach_source)
        
    for b in binnrs:
        sources_for_bin = np.where(bin_foreach_source == b)
        subcat = cat[sources_for_bin]
        subcat.write(catfile.replace('.fits', '_bin{}.fits'.format(b)))
    

if __name__ == '__main__':
    main()
