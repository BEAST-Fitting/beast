#!/usr/bin/env python
""" Make a nice plot of the PHAT filter response functions

"""
import argparse

from beast.plotting.plot_filters import plot_filters

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    args = parser.parse_args()

    # PHAT filter names
    filter_names = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
                    'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    out_names = ['F275W','F336W','F475W','F814W','F110W','F160W']

    plot_filters(args, filter_names, out_names)
