"""
Create extinguished grid

Jan 2013: Create by Karl G. to display the SED fitting results for a single SED

"""

__version__ = '0.1dev'

import numpy as np
from numpy import exp
import inspect
import itertools
import sys

import pylab
import matplotlib as mpl
import mytables
import grid

def plot_priors (stellar_filename):
    """
    plot the priors (just histogram the grid)
    INPUTS:
        stellar_filename string    FITS file with stellar SEDs (luminosities)
    KEYWORDS:
    """

    def OneDPlot(xvals, fig, subplot_val, xlabeltext, nbins=20):

        onedplot = fig.add_subplot(subplot_val)
        n, b = np.histogram(xvals, bins=nbins)
        n = n.astype(float)/n.sum()
        onedplot.step(b[:-1], n, where='pre')
        onedplot.set_ylim(0,n.max()*1.2)

        onedplot.set_xlabel(xlabeltext)
        
    # get the nD stellar/dust SED grid
    ext_grid = grid.FileSEDGrid(stellar_filename)

    # setup the plot
    fig = mpl.pylab.figure()

    # customize the font size
    pylab.rc('font', family='serif', size=10)

    # plot 1D histograms
    OneDPlot(ext_grid.logA, fig, 331, 'log(age)')
    OneDPlot(ext_grid.logM, fig, 332, 'log(mass)')
    OneDPlot(ext_grid.Z, fig, 333, 'Z')

    OneDPlot(ext_grid.Av, fig, 334, 'A(V)')
    OneDPlot(ext_grid.Rv, fig, 335, 'R(V)',nbins=10)
    OneDPlot(ext_grid.f_bump, fig, 336, 'f(bump)',nbins=10)


    OneDPlot(ext_grid.logT, fig, 338, 'log(Teff)')
    OneDPlot(ext_grid.logg, fig, 339, 'log(g)')

    # display
    mpl.pylab.show()

def plot_single_1d_fake (sed_results_filename, stellar_filename):
    """
    plot the results of fitting a single SED: 1D PDFs
    INPUTS:
        stellar_filename string    FITS file with stellar SEDs (luminosities)
    KEYWORDS:
    """

    def OneDPlot(xvals, weights, cor_val, fig, subplot_val, xlabeltext, nbins=20):

        onedplot = fig.add_subplot(subplot_val)
        n, b = np.histogram(xvals, weights = weights, bins=nbins)
        n = n.astype(float)/n.sum()
        onedplot.step(0.5*(b[:-1] + b[1:]), n, where='mid')
        onedplot.vlines(cor_val, 0, n.max()*1.5)
        onedplot.set_ylim(0,n.max()*1.2)

        onedplot.set_xlabel(xlabeltext)
        
    # get the nD likelihood function
    lnp = mytables.load(sed_results_filename, extension='LNP')

    # get the fake SED
    fakesed = mytables.load(sed_results_filename, extension='FAKESED')

    # get the nD stellar/dust SED grid
    ext_grid = grid.FileSEDGrid(stellar_filename)

    # get the index of the fake star
    fake_indx = lnp.header['FAKE_IDX']

    # setup the plot
    fig = mpl.pylab.figure()

    # customize the font size
    pylab.rc('font', family='serif', size=10)

    # plot 1D histograms
    exp_lnp = exp(lnp['lnp'])
    #indx = np.where(np.isfinite(exp_lnp))
    indx = lnp['idx'].astype(int)

    OneDPlot(ext_grid.logA[indx], exp_lnp, ext_grid.logA[fake_indx], fig, 331, 'log(age)')
    OneDPlot(ext_grid.logM[indx], exp_lnp, ext_grid.logM[fake_indx], fig, 332, 'log(mass)')
    OneDPlot(ext_grid.Z[indx], exp_lnp, ext_grid.Z[fake_indx], fig, 333, 'Z')

    OneDPlot(ext_grid.Av[indx], exp_lnp, ext_grid.Av[fake_indx], fig, 334, 'A(V)')
    OneDPlot(ext_grid.Rv[indx], exp_lnp, ext_grid.Rv[fake_indx], fig, 335, 'R(V)',nbins=10)
    OneDPlot(ext_grid.f_bump[indx], exp_lnp, ext_grid.f_bump[fake_indx], fig, 336, 'f(bump)',nbins=10)

    OneDPlot(ext_grid.logT[indx], exp_lnp, ext_grid.logT[fake_indx], fig, 338, 'log(Teff)')
    OneDPlot(ext_grid.logg[indx], exp_lnp, ext_grid.logg[fake_indx], fig, 339, 'log(g)')

    # plot the fake SED w/ and w/o noise
    sedplot = fig.add_subplot(337)
    plot_waves = fakesed['waves']/1e4
    sedplot.set_xlabel('$\lambda$ [$\mu$m]')
    sedplot.set_yscale('log')

    # plot the top 100 models w/ poorer fitting models becoming more transparent
    sindx = np.argsort(exp_lnp)
    topsindx = sindx[-10000:]

    max_lnp = max(exp_lnp[topsindx])
    for z in topsindx[::100]:
        sedplot.plot(plot_waves,ext_grid.seds[indx[z]],color='b',alpha=exp_lnp[z]/max_lnp)

    sedplot.errorbar(plot_waves,fakesed['fluxessed'],yerr=fakesed['fluxesunc'], fmt='o', color='r')
    sedplot.plot(plot_waves,fakesed['corfluxes'], color='r')

    # display
    mpl.pylab.show()

if __name__ == '__main__':

    stellar_filename='libs/stellib_kurucz2004_padovaiso.spectralgrid_sed_extinguished.grid.fits'

    plot_priors(stellar_filename)
    
    #plot_single_1d_fake(sys.argv[-1], stellar_filename)

