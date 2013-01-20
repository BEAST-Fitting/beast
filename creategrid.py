""" Create extinguished grid """
__version__ = '0.1dev'

import numpy
from numpy import exp
import inspect
import itertools
import sys

import mytables
from anased import computeLogLikelihood
from decorators import timeit
import extinction
import grid
import phot
import observations


def make_extinguished_grid(stellar_filename, filter_names, avs, rvs, fbumps):

    g0 = grid.FileSpectralGrid(stellar_filename)

    indx = numpy.arange(0,g0.grid.nrows, 100)
    g0.grid = g0.grid.extract(indx)
    g0.seds = g0.seds[indx]

    g0.write('small_grid.fits',clobber=True,append=False)

    return g0.getExtinguishedSEDs(filter_names, toPhotLamb=1.0, toFbumpLawLamb=1e-4, Av_vals=avs, Rv_vals=rvs, f_bump_vals=fbumps)

if __name__ == '__main__':

    # define the filename with the stellar grid
    stellar_filename = 'libs/stellib_kurucz2004_padovaiso.spectralgrid.fits'

    # define filters for the grid
    filter_names  = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    # grid spacing for stars
    #  TBD

    # grid spacing for dust
    avs = numpy.arange(0.0,2.0,2.)
    rvs = numpy.arange(1.0,6.0,5.0)
    fbumps = numpy.arange(0.0,1.0,1.0)

    # make the grid 
    extgrid = make_extinguished_grid(stellar_filename,filter_names, avs, rvs, fbumps)

    print extgrid.seds.shape

    # save grid to file
    extgrid.write(stellar_filename.replace('.fits','_extinguished.fits'),clobber=True,append=False)
