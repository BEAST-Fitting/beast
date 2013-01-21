"""
Perform a sensitivity test on a model grid

Jan 2013: written by Karl G.

"""

__version__ = '0.1dev'

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
import anased

def getFakeStar(g, idx, err=0.1):
    """ Generate a fake sed from a model grid
    INPUTS:
        idx int                         index number on the grid
    OUTPUTS:
        fakein  int                     the index of the model on the grid
        fakesed ndarray[float, ndim=1]  resulting SED

    KEYWORDS:
        err float                       proportional error to consider on the fluxes
    """
    fakein   = idx
    fakesed  = numpy.copy(g.seds[fakein, :])
    # sampled noise (with a mean of 1 and standard deviation of err (input value)
    unc_sample = numpy.random.normal(1.0,err,len(fakesed))
    fakesed *= unc_sample
    fakeerr = err * fakesed

    return fakein, fakesed, fakeerr


def fit_model_seds(N, n_test, stellar_filename, err=0.1, outdir='Tests/fake_many_0'):

    """
    Fit model seds with noise for sensitivity tests
    INPUTS:
        n                 int       number of models to test (randomly picked)
        n_test            int       number of noise realizations per model
        stellar_filename  string    FITS file with stellar SEDs (luminosities)
        err               float     fractional noise to assume
        outdir            string    output file directory for results
    """

    #Load the recomputed stellar+dust model grid
    ext_grid = grid.FileSEDGrid(stellar_filename)

    ## Initial SEDs are randomly drawn from the model space
    fakein = numpy.random.randint(0, ext_grid.grid.nrows, N)

    # mask for non-detectors (currently none)
    mask = numpy.zeros(len(ext_grid.lamb), dtype=bool)

    for tn in range(N):
        
        for tt in range(n_test):
            #fake DATA
            idx, fakesed, fakeerr = getFakeStar(ext_grid, fakein[tn], err=err)

            with timeit('Likelihood Object %d, Run %d' % (tn, tt) ):
                 #define output filename
                 outname = outdir + '/fake_star_%d_%d.fits' % ( tn, tt )

                 with timeit('\t * Computing Lnp'):
                     with timeit('\t * Direct Lnp computation'):
                         lnp = anased.computeLogLikelihood(fakesed, fakeerr, ext_grid.seds, normed=False, mask=mask)
                     with timeit('\t * Writing outputs'):
                         t = mytables.Table(name='LNP')
                         t.addCol(numpy.arange(ext_grid.grid.nrows, dtype=int), name='idx')
                         t.addCol(lnp, name='lnp')
                         t.header['GFNAME'] = stellar_filename
                         t.header['FAKE_IDX'] = idx
                         t.write( outname, clobber=True, append=False )

                         del lnp, t

                         # save the FAKE SED as another extension
                         t = mytables.Table(name='FAKESED')
                         t.addCol(ext_grid.lamb, name='waves')
                         t.addCol(fakesed, name='fluxessed')
                         t.addCol(fakeerr, name='fluxesunc')
                         t.addCol(ext_grid.seds[fakein[tn]], name='corfluxes')
                         t.write( outname, clobber=False, append=True )

                         del t

if __name__ == '__main__':

    # define the filename with the stellar grid
    stellar_filename = 'libs/stellib_kurucz2004_padovaiso.spectralgrid_sed_extinguished.grid.fits'

    import gc
    gc.enable()
    # generate the likelihoods for the fake stars (models w/ noise)
    fit_model_seds(1,5,stellar_filename,0.1)
    gc.disable()
