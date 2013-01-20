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


    def getExtinguishedSEDs(self, filter_names, toPhotLamb=1.0e4,toFbumpLawLamb=1.0,absFlux=True, Av_vals=numpy.array([3.0]), Rv_vals = numpy.array([2.7]), f_bump_vals = numpy.array([0.0])):
        """
        Extinguish and extract fluxes through filters
        INPUTS:
              filter_names   list    list of filter names according to the filter lib
        KEYWORDS:
              toPhotLamb     float   conversion factor from grid wavelength units to angstroms
              toFbumpLawLamb float   conversion factor from grid wavelength units to microns
              Av_vals   numpy array  Av values to iterate over
              Rv_vals   numpy array  Rv values to iterate over
              f_bump_vals numpy array f_bump values to iterate over
        """
        extLaw = extinction.RvFbumpLaw()
        
        if type(filter_names[0]) == str:
            flist = phot.load_filters(filter_names, interp=True, lamb=self.lamb*toPhotLamb)
        else:
            flist = filter_names

        lamb  = numpy.empty( len(flist), dtype=float)
        for e,k in enumerate(flist):
            lamb[e] = k.cl
        
        Av_vals, Rv_vals, f_bump_vals = numpy.ix_(Av_vals, Rv_vals, f_bump_vals)
        it = numpy.nditer([Av_vals,Rv_vals,f_bump_vals])
        ext_lamb = self.lamb*toFbumpLawLamb

        niter = Av_vals.size*Rv_vals.size*f_bump_vals.size
        cols = {'Av':numpy.empty(self.grid.nrows*niter),'Rv':numpy.empty(self.grid.nrows*niter),'f_bump':numpy.empty(self.grid.nrows*niter)}
            
        for key in self.keys():
            cols[key] = numpy.empty(self.grid.nrows*niter)

        results = MemoryGrid(lamb,seds= numpy.empty((self.grid.nrows*niter,filter_names.__len__())),grid=mytables.Table(iterable=cols,header=self.grid.header))

        itime = time()
        count = 0
        print "Memory allocated, starting on " + str(niter) + " iterations"
        for Av,Rv,f_bump in it:
            if ((count+1) % 100) == 0:
                print 100.*count/(Av_vals.size*Rv_vals.size*f_bump_vals.size),time()-itime
                itime = time()
            ext = numpy.exp(-1*extLaw.function(ext_lamb,Av=Av,Rv=Rv,f_bump=f_bump))
            results.grid['Av'][self.grid.nrows*count:self.grid.nrows*(count+1)] = Av
            results.grid['Rv'][self.grid.nrows*count:self.grid.nrows*(count+1)] = Rv
            results.grid['f_bump'][self.grid.nrows*count:self.grid.nrows*(count+1)] = f_bump
            results.seds[self.grid.nrows*count:self.grid.nrows*(count+1)] = phot.extractExtinguishedSEDs(self,flist,ext,absFlux=absFlux)
            for key in self.keys():
                results.grid[key][self.grid.nrows*count:self.grid.nrows*(count+1)] = self.grid[key]
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
    avs = numpy.arange(0.0,2.0,2.)
    rvs = numpy.arange(1.0,6.0,5.0)
    fbumps = numpy.arange(0.0,1.0,1.0)

    # make the grid 
    extgrid = make_extinguished_grid(stellar_filename,filter_names, avs, rvs, fbumps)

    # save grid to file
    extgrid.write(stellar_filename.replace('.fits','_extinguished.fits'),clobber=True,append=False)
