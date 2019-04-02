#!/usr/bin/env python
#
# remove filters from photometry catalogs, physicsgrid, and observationgrid
#   used to modify simulated data to make plots for proposals

import argparse

import numpy as np
from astropy.table import Table
import tables

from beast.physicsmodel.grid import (FileSEDGrid, SpectralGrid)
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel


def remove_filters_from_files(catfile,
                              physgrid,
                              obsgrid,
                              outbase,
                              rm_filters):

    # remove the requested filters from the catalog file
    cat = Table.read(catfile)
    for cfilter in rm_filters:
        colname = '{}_rate'.format(cfilter)
        if colname in cat.colnames:
            cat.remove_column(colname)
        else:
            print('{} not in catalog file'.format(colname))
    cat.write('{}_cat.fits'.format(outbase), overwrite=True)

    # get the sed grid and process
    g0 = FileSEDGrid(physgrid, backend='cache')
    filters = g0.header['filters'].split(' ')
    shortfilters = [(cfilter.split('_'))[-1].lower() for cfilter in filters]
    nlamb = []
    nfilters = []
    rindxs = []
    for csfilter, clamb, cfilter in zip(shortfilters, g0.lamb, filters):
        if csfilter not in rm_filters:
            nlamb.append(clamb)
            nfilters.append(cfilter)
        else:
            rindxs.append(shortfilters.index(csfilter))
    nseds = np.delete(g0.seds, rindxs, 1)

    print('orig filters: {}'.format(' '.join(filters)))
    print(' new filters: {}'.format(' '.join(nfilters)))

    g = SpectralGrid(np.array(nlamb), seds=nseds,
                     grid=g0.grid, backend='memory')
    g.grid.header['filters'] = ' '.join(nfilters)
    g.writeHDF('{}_sed.grid.hd5'.format(outbase))

    # get and process the observation model
    obsgrid = noisemodel.get_noisemodelcat(obsgrid)
    with tables.open_file('{}_noisemodel.grid.hd5'.format(outbase), 'w') \
            as outfile:
        outfile.create_array(outfile.root, 'bias',
                             np.delete(obsgrid.root.bias, rindxs, 1))
        outfile.create_array(outfile.root, 'error',
                             np.delete(obsgrid.root.error, rindxs, 1))
        outfile.create_array(outfile.root, 'completeness',
                             np.delete(obsgrid.root.completeness, rindxs, 1))


if __name__ == '__main__':

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("catfile",
                        help='filename of photometry catalog')
    parser.add_argument("physgrid",
                        help='filename of physics grid')
    parser.add_argument("obsgrid",
                        help='filename of observation/nosie grid')
    parser.add_argument("outbase", default='lessfilters',
                        help='filename for simulated observations')
    parser.add_argument('--rm_filters', type=str, nargs='*',
                        help='filters to remove')
    args = parser.parse_args()

    # do the merge
    remove_filters_from_files(args.catfile,
                              args.physgrid,
                              args.obsgrid,
                              args.outbase,
                              args.rm_filters)
