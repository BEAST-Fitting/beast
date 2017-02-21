#!/usr/bin/env python
"""
Library of general plotting functions for the BEAST plotting scripts

.. history::
    Written 21 Feb 2017 by Meredith J. Durbin
"""
import argparse
import re

def fancify_colname(name):
    ns = name.split('_')
    if re.match(re.compile('log._*'), name) is not None:
        fancyname = '$\log ({})$ ({})'.format(name[3], ns[-1])
    elif re.match(re.compile('.v_*'), name) is not None:
        fancyname = '${}_V$ ({})'.format(name[0], ns[-1])
    elif re.match(re.compile('._._*'), name) is not None:
        fancyname = '${}_{}$ ({})'.format(ns[0], ns[1], ns[-1])
    elif name == 'chi2min':
        fancyname = '$\chi^2$'
    else:
        fancyname = ' '.join(namesplit)
    return r'{}'.format(fancyname)

def initialize_parser():
    ftypes = ['png', 'jpg', 'jpeg', 'pdf', 'ps', 'eps', 'rgba',
              'svg', 'tiff', 'tif', 'pgf', 'svgz', 'raw']
    ftypestr = ', '.join(['"{}"'.format(f) for f in ftypes])
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--savefig', action='store', default=False, choices=ftypes,
                        help='Save figure to a file of specified type rather than show it. \
                              Must be one of: {}'.format(ftypestr)
                        )
    parser.add_argument('-t', '--tex', action='store_true', 
                        help='Configure Matplotlib to use LaTeX; \
                              sets rcParam "usetex" to True. \
                              (Requires working TeX installation.)\
                              Defaults to false.')
    return parser