"""
Create extinguished grid

Dec 2012: Originally written by Kirill T.
Jan 2013: Modified by Karl G. to separate the extinguished grid creation into a separate program.

"""

__version__ = '0.2dev'

import numpy
from numpy import exp
import inspect
import itertools
import sys

import mytables
from decorators import timeit
import extinction
import grid
import phot
#import observations
from time import time
import progressbar

def make_extinguished_grid(stellar_filename, filter_names, avs, rvs, fbumps):

    """
    Extinguish and extract fluxes through filters
       (all wavelengths in stellar SEDs and filter response functions assumed to be in Angstroms)
    INPUTS:
        stellar_filename string    FITS file with stellar SEDs (luminosities)
        filter_names   list        list of filter names according to the filter lib
    KEYWORDS:
        Av_vals     numpy array    Av values to iterate over
        Rv_vals     numpy array    Rv values to iterate over
        f_bump_vals numpy array    f_bump values to iterate over
    """

    # get the stellar grid (no dust)
    g0 = grid.FileSpectralGrid(stellar_filename)

    # reduce size of grid for debugging (TBR)
    indx = numpy.arange(0,g0.grid.nrows, 100)
    g0.grid = g0.grid.extract(indx)
    g0.seds = g0.seds[indx]
    g0.write('small_grid.fits',clobber=True,append=False)

    # define the extinction law be be used (fixed for now)
    extLaw = extinction.RvFbumpLaw()

    # ensure that dust parameter verctors are numpy arrays
    Av_vals, Rv_vals, f_bump_vals = numpy.ix_(avs, rvs, fbumps)

    # setup interation over the full dust parameter grid
    it = numpy.nditer([Av_vals,Rv_vals,f_bump_vals])
    niter = Av_vals.size*Rv_vals.size*f_bump_vals.size

    # setup of output
    cols = {'Av':numpy.empty(g0.grid.nrows*niter),'Rv':numpy.empty(g0.grid.nrows*niter),'f_bump':numpy.empty(g0.grid.nrows*niter)}

    for key in g0.keys():
        cols[key] = numpy.empty(g0.grid.nrows*niter)

    itime = time()
    count = 0
    with progressbar.PBar(niter, txt='creategrid') as Pbar:
        for Av,Rv,f_bump in it:
            Pbar.update(count)
            # info showing program is running
            #if ((count+1) % 100) == 0:
            #    print 100.*count/(Av_vals.size*Rv_vals.size*f_bump_vals.size),time()-itime
            #    itime = time()

            # apply extinction and integrate over band response functions
            temp_results = g0.getSEDs(filter_names, extLaw=extLaw, Av=Av, Rv=Rv, f_bump=f_bump)

            # setup the output object
            #  must be done after the 1st extraction to get the wavelength vector
            if count is 0:
                #print "Memory allocated"
                results = grid.MemoryGrid(temp_results.lamb,seds= numpy.empty((g0.grid.nrows*niter,filter_names.__len__())),grid=mytables.Table(iterable=cols,header=g0.grid.header))

            # assign the extinguished SEDs to the output object
            results.seds[g0.grid.nrows*count:g0.grid.nrows*(count+1)] = temp_results.seds

            # free the memory of temp_results
            del temp_results

            # adding the dust parameters to the models
            results.grid['Av'][g0.grid.nrows*count:g0.grid.nrows*(count+1)] = Av
            results.grid['Rv'][g0.grid.nrows*count:g0.grid.nrows*(count+1)] = Rv
            results.grid['f_bump'][g0.grid.nrows*count:g0.grid.nrows*(count+1)] = f_bump

            # the rest of the parameters
            for key in g0.keys():
                results.grid[key][g0.grid.nrows*count:g0.grid.nrows*(count+1)] = g0.grid[key]
            count += 1

    return results

if __name__ == '__main__':

    # define the filename with the stellar grid
    stellar_filename = 'libs/stellib_kurucz2004_padovaiso.spectralgrid.fits'

    # define filters for the grid
    filter_names  = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    # grid spacing for stars
    #  TBD

    # grid spacing for dust
    avs = numpy.arange(0.0,2.0,0.5)
    rvs = numpy.arange(1.0,6.0,1.0)
    fbumps = numpy.arange(0.0,1.0,0.1)

    # make the grid 
    extgrid = make_extinguished_grid(stellar_filename,filter_names, avs, rvs, fbumps)

    # save grid to file
    extgrid.write(stellar_filename.replace('.fits','_extinguished.fits'),clobber=True,append=False)
