#!/usr/bin/env python
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
    parser.add_argument('filename', type=str,
                        help="name of file with evolutionary track")
    parser.add_argument('--type', default='mist', choices=['mist', 'parsec'],
                        help="source of tracks")
    args = parser.parse_args()

    fig, ax = plt.subplots(2, 2, figsize=(10, 10))

    if args.type == 'parsec':
        et = ETParsec()
    else:
        et = ETMist()
    et.load_orig_tables(args.filename)

    et.plot_tracks(ax[0, 0], xval='logT', yval='logL')
    et.plot_tracks(ax[1, 0], xval='logA', yval='phase')
    et.plot_tracks(ax[0, 1], xval='logA', yval='eep')
    et.plot_tracks(ax[1, 1], xval='logT', yval='logg')

    fig.tight_layout()

    save_name = ''
    if args.tex:
        plt.rc({'usetex': True})
    if args.savefig:
        fig.savefig('{}.{}'.format(save_name, args.savefig))
    else:
        plt.show()
