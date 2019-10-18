"""
Make a plot of the evolutionary tracks
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib.pyplot as plt

from beast.physicsmodel.stars.evoltracks import (ETParsec, ETMist)
from beast.plotting.beastplotlib import initialize_parser


if __name__ == '__main__':
    parser = initialize_parser()
    parser.add_argument('filename', type=str, nargs='*',
                        help="name of file(s) with evolutionary track")
    parser.add_argument('--type', default='mist', choices=['mist', 'parsec'],
                        help="source of tracks")
    args = parser.parse_args()

    fig, ax = plt.subplots(2, 2, figsize=(10, 10))

    # switch between the types of tracks
    if args.type == 'parsec':
        et = ETParsec()
    else:
        et = ETMist()

    # read in the tracks
    et.load_orig_tables(args.filename)

    orig_metrics = et.grid_metrics()

    # plot, plot, plot
    nls = ':'
    et.plot_tracks(ax[0, 0], xval='logT', yval='logL', linestyle=nls)
    et.plot_tracks(ax[1, 1], xval='logA', yval='M_ini', linestyle=nls)
    et.plot_tracks(ax[0, 1], xval='phase', yval='M_ini', linestyle=nls)
    et.plot_tracks(ax[1, 0], xval='logT', yval='logg', linestyle=nls)

    # regrid the evolutionary tracks to uniform log(mass) and variable age
    print('size orig = ', len(et.data['M_ini']))

    et.regrid(logmass_range=[-1., 3.],
              logmass_delta=0.05,
              logT_delta=0.05,
              logL_delta=0.05)

    print('size regrid = ', len(et.data['M_ini']))

    regrid_metrics = et.grid_metrics()
    for ckey in regrid_metrics.keys():
        print(ckey)
        print(orig_metrics[ckey])
        print(regrid_metrics[ckey])

    # get the new grid metrics
    # et.grid_metrics()

    # plot, plot, plot
    nls = '-'
    et.plot_tracks(ax[0, 0], xval='logT', yval='logL', linestyle=nls)
    et.plot_tracks(ax[1, 1], xval='logA', yval='M_ini', linestyle=nls)
    et.plot_tracks(ax[0, 1], xval='phase', yval='M_ini', linestyle=nls)
    et.plot_tracks(ax[1, 0], xval='logT', yval='logg', linestyle=nls)

    fig.tight_layout()

    save_name = 'evoltracks'
    if args.tex:
        plt.rc({'usetex': True})
    if args.savefig:
        fig.savefig('{}.{}'.format(save_name, args.savefig))
    else:
        plt.show()
