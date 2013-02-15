"""
Perform a sensitivity test on a model grid

Jan 2013: written by Karl G.

"""

__version__ = '0.1dev'

import numpy
import progressbar
import eztables
import tables

#from decorators import timeit
#import extinction
import grid
#import phot
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
    unc_sample = numpy.random.normal(1.0, err, len(fakesed))
    fakesed *= unc_sample
    fakeerr = err * fakesed
    return fakein, fakesed, fakeerr


def getFakeInds(ext_grid, logM, logT, delt_logM, delt_logT):
    F814_i = 3
    flux_sel = numpy.where(ext_grid.seds[:,F814_i] > 3.e-18)[0]
    sel = numpy.where((numpy.abs(ext_grid.logM[flux_sel] - logM) < delt_logM)*(numpy.abs(ext_grid.logT[flux_sel] - logT) < delt_logT))[0]
    return flux_sel[sel]


def fit_model_seds(N, n_test, stellar_filename, err=0.1, outdir='Tests/fake_many_0'):

    """
    Fit model seds with noise for sensitivity tests
    INPUTS:
        N                 int       number of models to test (randomly picked)
        n_test            int       number of noise realizations per model
        stellar_filename  string    FITS file with stellar SEDs (luminosities)
        err               float     fractional noise to assume
        outdir            string    output file directory for results
    """


    #Load the recomputed stellar+dust model grid
    ext_grid = grid.FileSEDGrid(stellar_filename)

    # Initial SEDs to be uniformaly spaced throughout the grid
    #fakein = numpy.round(numpy.linspace(0, ext_grid.grid.nrows - 1, num = N)).astype(int)
    #Randomly chosen from main sequence
    MS = where((ext_grid.logT>4.0) * (ext_grid.logL*4.5 - 15 < ext_grid.logT))[0]
    fakein = MS[numpy.random.randint(low=0,high=MS.size,size=N)]

    # mask for non-detectors (currently none)
    mask = numpy.zeros(len(ext_grid.lamb), dtype=bool)
    with progressbar.PBar(N,txt="Calculating lnp") as pbar:
        for tn in range(N):
            for tt in range(n_test):
            #fake DATA
                idx, fakesed, fakeerr = getFakeStar(ext_grid, fakein[tn], err=err)

                #define output filename
                outname = outdir + '/fake_star_%d_%d.fits' % ( tn, tt )

                lnp = anased.computeLogLikelihood(fakesed, fakeerr, ext_grid.seds, normed=False, mask=mask)

                indx = numpy.where((lnp - max(lnp)) > -40.)

                t = eztables.Table(name='LNP')
                t.addCol('idx',numpy.array(indx[0],dtype=numpy.int32))
                t.addCol('lnp',numpy.array(lnp[indx],dtype=numpy.float))
                t.header['GFNAME'] = stellar_filename
                t.header['FAKE_IDX'] = idx
                t.write( outname,type='fits',clobber=True, append=False,silent=True)
                del lnp, t

                # save the FAKE SED as another extension
                t = eztables.Table(name='FAKESED')
                t.addCol('waves',ext_grid.lamb)
                t.addCol('fluxessed',fakesed)
                t.addCol('fluxesunc',fakeerr)
                t.addCol('corfluxes',ext_grid.seds[fakein[tn]])
                t.write( outname, clobber=False, append=True,silent=True)
                del t
            pbar.update(tn)

def fit_model_seds_pytables(fakein, N, n_test, stellar_filename, filters = numpy.arange(6), err=0.1, outname='Tests/fake_many_0/test1.hf5'):

    """
    Fit model seds with noise for sensitivity tests
    INPUTS:
        N                 int       number of models to test (randomly picked)
        n_test            int       number of noise realizations per model
        stellar_filename  string    FITS file with stellar SEDs (luminosities)
        err               float     fractional noise to assume
        outdir            string    output file directory for results
    """
    #Load the recomputed stellar+dust model grid
    ext_grid = grid.FileSEDGrid(stellar_filename)

    # Initial SEDs to be uniformaly spaced throughout the grid
    #fakein = numpy.round(numpy.linspace(0, ext_grid.grid.nrows - 1, num = N)).astype(int)
    # Pick randomly spaced SEDs
    #fakein = numpy.random.randint(ext_grid.grid.nrows,size=N)

    ext_grid.seds = ext_grid.seds[:,filters]
    ext_grid.lamb = ext_grid.lamb[:,filters]

    # mask for non-detectors (currently none)
    mask = numpy.zeros(len(ext_grid.lamb), dtype=bool)
    outfile = tables.openFile(outname, 'w')
    outfile.createArray(outfile.root,'waves',ext_grid.lamb) #Save wavelengths in root, remember
                                                            #n_stars = root._v_nchildren -1
    with progressbar.PBar(N,txt="Calculating lnp") as pbar:
        for tn in range(N):
            star_group = outfile.createGroup('/','fakeStar_%d' %tn, title="Fake star %d" %tn)
            star_group._v_attrs.fakein =  fakein[tn]
            for tt in range(n_test):
                #fake DATA
                idx, fakesed, fakeerr = getFakeStar(ext_grid, fakein[tn], err=err)

                lnp = anased.computeLogLikelihood(fakesed, fakeerr, ext_grid.seds, normed=False, mask=mask)

                fake_group = outfile.createGroup(star_group,'fake_%d' %tt)
                fake_group._v_attrs.test_n = tt
                #Need ragged arrays rather than uniform table
                outfile.createArray(fake_group,'fakesed',fakesed)
                outfile.createArray(fake_group,'fakerr',fakeerr)
                indx = numpy.where((lnp - max(lnp)) > -40.)

                outfile.createArray(fake_group,'idx',numpy.array(indx[0],dtype=numpy.int32))
                outfile.createArray(fake_group,'lnp',numpy.array(lnp[indx[0]],dtype=numpy.float32))

            outfile.flush()
            pbar.update(tn)
    outfile.close()

if __name__ == '__main__':

    # define the filename with the stellar grid
    stellar_filename = 'no_smc.fits'

    # generate the likelihoods for the fake stars (models w/ noise)
    fit_model_seds_pytables(10,2,stellar_filename,0.1,outname='Tests/all_phat.hf5')
    cut = 3.e-18

