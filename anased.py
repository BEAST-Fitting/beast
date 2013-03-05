__version__ = '0.1dev'

from .config import __NTHREADS__
from .config import __USE_NUMEXPR__
if __USE_NUMEXPR__:
    import numexpr
    numexpr.set_num_threads(__NTHREADS__)
else:
    #import numpy ufunc (c-coded for speed-up)
    from numpy import subtract, divide, power

import numpy
from numpy import log, log10, exp

from . import stellib
from . import extinction
from . import photometry
from tools.decorators import timeit


def fluxToMag(flux):
    """ Return the magnitudes from flux values
    INPUTS:
        flux    np.ndarray[float, ndim=N]   array of fluxes
    OUTPUTS:
        mag np.ndarray[float, ndim=N]   array of magnitudes
    """
    return -2.5 * log10(flux)


def fluxErrTomag(flux, fluxerr):
    """ Return the magnitudes and associated errors from fluxes and flux error values
    INPUTS:
        flux    np.ndarray[float, ndim=1]   array of fluxes
        fluxerr np.ndarray[float, ndim=1]   array of flux errors
    OUTPUTS:
        mag np.ndarray[float, ndim=1]   array of magnitudes
        err np.ndarray[float, ndim=1]   array of magnitude errors
    """
    mag = fluxToMag(flux)
    return mag, -2.5 * log10( 1. - fluxerr / flux )


def magToFlux(mag):
    """ Return the flux from magnitude values
    INPUTS:
        mag np.ndarray[float, ndim=N]   array of magnitudes
    OUTPUTS:
        flux    np.ndarray[float, ndim=N]   array of fluxes
    """
    return power(10, -0.4 * mag)


def magErrToFlux(mag, err):
    """ Return the flux and associated errors from magnitude and mag error values
    INPUTS:
        mag np.ndarray[float, ndim=1]   array of magnitudes
        err np.ndarray[float, ndim=1]   array of magnitude errors
    OUTPUTS:
        flux    np.ndarray[float, ndim=1]   array of fluxes
        fluxerr np.ndarray[float, ndim=1]   array of flux errors
    """
    flux = magToFlux(mag)
    return flux, flux * ( 1. - magToFlux(err) )


def getFluxAttenuation(law, lamb, **kwargs):
    """ Get flux attenuations from a given extinciton law
    INPUTS:
        law extinction.ExtinctionLaw    instance of extinction law
        lamb    np.ndarray[float, ndim=1]   array of wavelengths in AA
    KEYWORDS:
        **kwargs is forwarded to the call of the law instance: law(lamb, **kwargs)
    OUTPUTS:
        tau     np.ndarray[float, ndim=1]   tau as in redflux = flux*exp(-tau)
    """
    assert(isinstance(law, extinction.ExtinctionLaw))
    r = law.function( lamb * 1e-4, Alambda=False, **kwargs)
    return r


def computeChi2(flux, fluxerr, fluxmod):
    """ Compute the non-reduced chi2 between data with uncertainties and
        perfectly known models
    INPUTS:
        flux    np.ndarray[float, ndim=1]   array of fluxes
        fluxerr np.ndarray[float, ndim=1]   array of flux errors
        fluxmod np.ndarray[float, ndim=2]   array of modeled fluxes (Nfilters , Nmodels)
    OUTPUTS:
        chi2    np.ndarray[float, ndim=1]   array of chi2 values (Nmodels)

    Note: using ufunc because ufuncs are written in C (for speed) and linked into NumPy.
    """
    # make sure errors are not null
    fluxerr[fluxerr == 0.] = 1.
    if __USE_NUMEXPR__:
        return numexpr.evaluate('sum(((flux-fluxmod)/fluxerr)**2,axis=1)',
                                local_dict={'flux': flux, 'fluxmod': fluxmod, 'fluxerr': fluxerr})
    else:
        return power( divide( subtract(flux[None, :], fluxmod), fluxerr[None, :]), 2.).sum(axis=1)


def computeLogLikelihood(flux, fluxerr, fluxmod, normed=True, mask=None, lnp_threshold=50.):
    """ Compute the log of the chi2 likelihood between data with uncertainties and
        perfectly known models
    INPUTS:
        flux    np.ndarray[float, ndim=1]   array of fluxes
        fluxerr np.ndarray[float, ndim=1]   array of flux errors
        fluxmod np.ndarray[float, ndim=2]   array of modeled fluxes (Nfilters , Nmodels)
    KEYWORDS:
        normed  bool                if set normalize the result
        mask    np.ndarray[bool, ndim=1]    mask array to apply during the calculations
                                mask.shape = flux.shape
        lnp_threshold               cut the values outside -x, x in lnp
    OUTPUTS:
        ln(L)   np.ndarray[float, ndim=1]   array of ln(L) values (Nmodels)

        with L = 1/[sqrt(2pi) * sig**2 ] * exp ( - 0.5 * chi2 )
             L propto  exp ( - 0.5 * chi2 )
    Note: using ufunc because ufuncs are written in C (for speed) and linked into NumPy.
    """
    if __USE_NUMEXPR__:
        fluxerr = numexpr.evaluate('where((fluxerr==0.), 1., fluxerr)', local_dict={'fluxerr': fluxerr})
        flux = numexpr.evaluate('where((flux==0.), 1e-5, flux)', local_dict={'flux': flux})
    else:
        fluxerr[fluxerr == 0.] = 1.
        flux[flux == 0.] = 1e-5

    if not mask is None:
        _m   = ~mask
        dof = _m.sum()
        if __USE_NUMEXPR__:
            chi2 = numexpr.evaluate('- 0.5 / b * a', local_dict={'a': computeChi2( flux[_m], fluxerr[_m], fluxmod[:, _m] ), 'b': dof})
        else:
            chi2 = - 0.5 / dof * computeChi2( flux[_m], fluxerr[_m], fluxmod[:, _m] )

    else:
        dof = len(flux)
        if __USE_NUMEXPR__:
            chi2 = numexpr.evaluate('- 0.5 / b * a', local_dict={'a': computeChi2( flux, fluxerr, fluxmod ), 'b': dof})
        else:
            chi2 = - 0.5 / dof * computeChi2( flux, fluxerr, fluxmod )

    #taking care of possible overflows...
    #if __USE_NUMEXPR__:
    #    chi2 = numexpr.evaluate('where( (chi2 > lim), lim,chi2)', local_dict={'chi2': chi2, 'lim': lnp_threshold})
    #    chi2 = numexpr.evaluate('where( (chi2 < -lim), -lim,chi2)', local_dict={'chi2': chi2, 'lim': lnp_threshold})
    #else:
    #    chi2 = numpy.clip( chi2, -lnp_threshold, lnp_threshold )

    if normed is True:
        if __USE_NUMEXPR__:
            expchi2 = numexpr.evaluate('exp(chi2)', local_dict={'chi2': chi2})
            # not really sure take works as expected with the inf values...
            # if it does, then if I follow the documentation we only need:
            # expchi2 = expchi2.take(numpy.isfinite(expchi2), mode='clip')
            expchi2[numpy.isinf(expchi2)] = expchi2.take(numpy.isfinite(expchi2), mode='clip').max()
        else:
            expchi2 = exp(chi2)
            expchi2[numpy.isinf(expchi2)] = expchi2[ numpy.isfinite(expchi2) ].max()
        psum = expchi2.sum()
        if __USE_NUMEXPR__:
            lnp = numexpr.evaluate('chi2 - log(psum)', local_dict={'chi2': chi2, 'psum': psum})
        else:
            lnp = chi2 - log(psum)
    else:
        lnp = chi2

    if __USE_NUMEXPR__:
        return numexpr.evaluate('where( (lnp < -1000), -1000, lnp)', local_dict={'lnp': lnp})
    else:
        lnp[ lnp < -1000] = -1000
        return lnp


def multi_job(lamb, flux, fluxerr, mask, fluxmod, extLaw, **kwargs):
    """ Shortcut to compute the log likelihood of multiple SEDs with the models
        for a given extinction parameter set.
    INPUTS:
        lamb    np.ndarray[float, ndim=1]   array of wavelengths in AA (Nfilters)
        flux    np.ndarray[float, ndim=2]   array of fluxes (Nfilters , Nobs)
        fluxerr np.ndarray[float, ndim=2]   array of flux errors (Nfilters , Nobs)
        mask    np.ndarray[bool, ndim=1]    mask array to apply during the calculations
                                mask.shape = flux.shape
        fluxmod np.ndarray[float, ndim=2]   array of modeled fluxes (Nfilters , Nmodels)
        extLaw  extinction.ExtinctionLaw    instance of extinction law
    KEYWORDS:
        **kwargs is forwarded to the getFluxAttenuation call
    OUTPUTS:
        ln(L)   np.ndarray[float, ndim=2]   array of ln(L) values (Nobs, Nmodels)
    """

    #get attenuation
    tau = getFluxAttenuation(extLaw, lamb, **kwargs)
    tau = exp(tau)  # less operations in the loop

    #compute lnp
    lnp = numpy.empty( ( flux.shape[0], fluxmod.shape[0] ), dtype=float)

    #This loop can be //
    for k in xrange(flux.shape[0]):
        deredflux = flux[k, :] * tau
        lnp[k, :] = computeLogLikelihood(deredflux, fluxerr[k, :], fluxmod, mask=mask)

    return lnp


def job(lamb, flux, fluxerr, mask, fluxmod, extLaw, **kwargs):
    """ Shortcut to compute the log likelihood of the SED with the models
        for a given extinction parameter set.
    INPUTS:
        lamb    np.ndarray[float, ndim=1]   array of wavelengths in AA
        flux    np.ndarray[float, ndim=1]   array of fluxes
        fluxerr np.ndarray[float, ndim=1]   array of flux errors
        mask    np.ndarray[bool, ndim=1]    mask array to apply during the calculations
                                mask.shape = flux.shape
        fluxmod np.ndarray[float, ndim=2]   array of modeled fluxes (Nfilters , Nmodels)
        extLaw  extinction.ExtinctionLaw    instance of extinction law
    KEYWORDS:
        **kwargs is forwarded to the getFluxAttenuation call
    OUTPUTS:
        ln(L)   np.ndarray[float, ndim=1]   array of ln(L) values (Nmodels)
    """

    #get attenuation
    tau = getFluxAttenuation(extLaw, lamb, **kwargs)

    #deredden the observed flux (faster than adding reddening to all models
    deredflux = flux * exp(tau)

    ind = (deredflux > 0.)
    deredflux[ind] = deredflux[ind]

    ind = (fluxmod > 0.)
    fluxmod[ind] = fluxmod[ind]
    #compute lnp
    lnp = computeLogLikelihood(deredflux, fluxerr, fluxmod, normed=False, mask=mask)
    #expchi2 = exp(lnp)
    #expchi2[numpy.isinf(expchi2)] = expchi2[ numpy.isfinite(expchi2) ].max()
    #psum = expchi2.sum()
    #lnp = lnp - log(psum)

    return lnp


def getSEDs(filter_names, lamb, specs):
    """
    Extract integrated fluxes through filters
    INPUTS:
        filter_names    list    list of filter names according to the filter lib
    """
    flist = photometry.load_filters(filter_names)

    r = numpy.empty( (len(specs), len(flist) ), dtype=float)
    lf = numpy.empty( len(flist), dtype=float )

    for kf in range(len(flist)):
        for ks in range(len(specs)):
            r[ks, kf] = flist[kf].getFlux(lamb, specs[ks])
            lf[kf] = flist[kf].cl

    return lf, r


def test_specs():
    from tools import figure

    osl = stellib.BaSeL()
    oAv = extinction.Cardelli()
    #fake DATA
    #fakein  = 2000 # random between 0 & 4523, no idea what this is :p
    idx      = osl.grid.where('(Teff >= 3.6) & (Teff <= 3.8) & (logG >= 4.5) & (logG <= 4.7) & (Z == 0.02)')
    fakein   = idx[0][0]
    fakesed  = numpy.copy(osl.spectra[fakein, :])
    Av0      = 0.1
    lamb     = osl.wavelength
    tau      = getFluxAttenuation(oAv, lamb, Av=Av0, Rv=3.1)
    fakesed *= exp(-tau)
    #magerr  = 0.05
    #fakeerr = fakesed * (1. - 10**(-0.4*magerr) )
    fakeerr = 0.5 * fakesed

    #get Models
    # will be replaced by broad-band SEDs but structure will be identical
    seds = numpy.copy(osl.spectra)
    #idx = osl.grid.where('(Z == 0.02)')
    #seds = osl.spectra[idx]
    lamb = osl.wavelength

    Av = numpy.arange(0, 1, 0.1)
    r = numpy.empty( (seds.shape[0], len(Av)), dtype=float )
    with timeit('Likelihood'):
        for k in range(len(Av)):
            r[:, k] = job(lamb, fakesed, fakeerr, seds, oAv, Av=Av[k], Rv=3.1)

    def plot(_r, idx=None):
        import pylab as plt
        if _r.ndim == 2:
            r = _r.sum(1)
        else:
            r = _r
        if idx is None:
            idx = numpy.arange(len(r))
        n0, bT, bg = numpy.histogram2d(osl.Teff[idx], osl.logg[idx], bins=[25, 11])
        n,  bT, bg = numpy.histogram2d(osl.Teff[idx], osl.logg[idx], bins=[bT, bg], weights=exp(r) )

        n0 = n0.astype(float) / n0.sum()
        n  = n.astype(float) / n.sum()
        n1 = numpy.copy(n[:])

        ind = (n0 > 0.)
        n[ind] /= n0[ind]

        n /= n.sum()

        n1 = numpy.ma.masked_where( n0 == 0, n1 )
        n  = numpy.ma.masked_where( n0 == 0, n  )

        plt.figure(1, figsize=(10, 10))
        plt.clf()

        ax0 = plt.subplot(221)
        ax0.imshow(n1.T, extent=[min(bT), max(bT), min(bg), max(bg)],
                vmin=0., vmax=numpy.max([n1.max(), n.max()]),
                origin='lower', aspect='auto')
        ax0.plot([osl.Teff[fakein]], [osl.logg[fakein]], 'o', mec='#ff0000', mfc='None', mew=2., ms=10.)
        ax0.set_xlabel('logT')
        ax0.set_ylabel('logg')
        ax0.set_xlim(ax0.get_xlim()[::-1])
        ax0.set_ylim(ax0.get_ylim()[::-1])
        ax0.set_title('Raw: P0(logT,logg) = K')

        ax1 = plt.subplot(222, sharex=ax0, sharey=ax0)
        ax1.imshow(n.T, extent=[min(bT), max(bT), min(bg), max(bg)],
                vmin=0., vmax=numpy.max([n1.max(), n.max()]),
                origin='lower', aspect='auto')
        ax1.plot([osl.Teff[fakein]], [osl.logg[fakein]], 'o', mec='#ff0000', mfc='None', mew=2., ms=10.)
        ax1.set_xlabel('logT')
        ax1.set_ylabel('logg')
        ax1.set_xlim(ax1.get_xlim()[::-1])
        ax1.set_ylim(ax1.get_ylim()[::-1])
        ax1.set_title('Corrected: P/P0(logT,logg) = K')

        ax2 = plt.subplot(223)
        x = 0.5 * (bT[1:] + bT[:-1])
        #ax2.step( bT[:-1], n1.sum(1), where='pre', lw=2.)
        #ax2.step( bT[:-1], n.sum(1),  where='pre', lw=2.)
        ax2.plot( x, n1.sum(1), lw=2., label='raw')
        ax2.plot( x, n.sum(1),  lw=2., label='cor')
        ylim = ax2.get_ylim()
        ax2.vlines([osl.Teff[fakein]], ylim[0], ylim[1])
        ax2.set_ylim(ylim)
        ax2.set_xlabel('log(Teff)')
        ax2.set_ylabel(r'P( data $\mid$ log(Teff) ) P(Teff) / P0(log(Teff))')
        ax2.legend(loc=0, frameon=False, borderaxespad=2., prop={'size': 14})

        ax3 = plt.subplot(224)
        x = 0.5 * (bg[1:] + bg[:-1])
        #ax3.step( bg[:-1], n1.sum(0), where='pre', lw=2.)
        #ax3.step( bg[:-1], n.sum(0),  where='pre', lw=2.)
        ax3.plot( x, n1.sum(0), lw=2.)
        ax3.plot( x, n.sum(0),  lw=2.)
        ylim = ax3.get_ylim()
        ax3.vlines([osl.logg[fakein]], ylim[0], ylim[1])
        ax3.set_ylim(ylim)
        ax3.set_xlabel('log(g)')
        ax3.set_ylabel(r'P( data $\mid$ log(g) ) P(log(g)) / P0(log(g))')

        figure.theme(ax=ax0)
        figure.theme(ax=ax1)
        figure.theme(ax=ax2)
        figure.theme(ax=ax3)

    plot(r)
    return r


if __name__ == '__main__':
    pass


def getFake(g, Av0=1., Rv0=3.1, err=0.1):
    oAv      = extinction.Cardelli()
    #idx      = g.grid.where('(logT >= 3.6) & (logT <= 3.8) & (logg >= 4.5) & (logg <= 4.7) & (Z == 0.02)')
    idx      = g.grid.where('(logT >= 3.95) & (logT <= 4.05) & (logg >= 1.85) & (logg <= 1.95) ')
    lamb     = g.lamb
    fakein   = idx[0][0]
    fakesed  = numpy.copy(g.seds[fakein, :])
    tau      = getFluxAttenuation(oAv, lamb, Av=Av0, Rv=Rv0)
    fakesed *= exp(-tau)
    #magerr  = 0.05
    #fakeerr = fakesed * (1. - 10**(-0.4*magerr) )
    fakeerr  = err * fakesed
    return fakein, lamb, fakesed, fakeerr


def test_seds(err=0.1, Av0=1., Z0=0.02):
    import grid
    #filters  = 'hst_wfc3_f225w hst_wfc3_f336w hst_acs_hrc_f475w hst_acs_hrc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    #osl = stellib.BaSeL()
    oAv = extinction.Cardelli()
    g = grid.FileSpectralGrid('libs/SEDs_basel_padovaiso.fits')
    lamb = g.lamb  # *1e6

    #fake DATA
    fakein, l, fakesed, fakeerr = getFake(g, Av0, 3.1, err=err)

    mask = numpy.zeros(fakesed.shape, dtype=bool)
    mask[3] = True
    mask[2] = True

    Av = numpy.arange(0., 3., 0.1)
    r = numpy.empty( (g.seds.shape[0], len(Av)), dtype=float )
    with timeit('Likelihood'):
        for k in range(len(Av)):
            r[:, k] = job(lamb[:], numpy.copy(fakesed), numpy.copy(fakeerr), mask, numpy.copy(g.seds), oAv, Av=Av[k], Rv=3.1)

    return g, r, Av, fakein, lamb, fakesed, fakeerr, Av0, Z0


def plotPDFs(g, r, Av, Av0, Z0, fakein, Q='logg logT logL logM logA Av Z' ):

    from tools import figure

    _q = Q.split()

    _r = exp(r)
    _r /= _r.sum()

    def plotPDF(ax, qk, *args, **kargs):
        if qk.upper() != 'AV':
            ax.hist(g.grid[qk], weights=_r.sum(1), bins=30)
            ax.set_xlabel(qk)
            ax.set_ylabel('P(data $\mid$ %s)' % qk )
            ylim = ax.get_ylim()
            ax.vlines(g.grid[qk][fakein], ylim[0], ylim[1], color='#ff0000')
            ax.set_ylim(ylim)
            ax.set_xlim(min(g.grid[qk]), max(g.grid[qk]))
        else:
            ax.hist(Av, weights=_r.sum(0), bins=numpy.linspace(min(Av), max(Av), len(Av) + 1))
            ax.set_xlabel(qk)
            ax.set_ylabel('P(data $\mid$ %s)' % qk )
            ax.set_xlim(min(Av), max(Av))
            ylim = ax.get_ylim()
            ax.vlines(Av0, ylim[0], ylim[1], color='#ff0000')
            ax.set_ylim(ylim)

    ncol = 3
    nl   = len(_q) / ncol + len(_q) % ncol
    k = 0
    for k in range(len(_q)):
        print _q[k] + ": " + str(g.grid[_q[k]][fakein])
        ax = figure.subplot(nl, ncol, k + 1)
        plotPDF(ax, _q[k])
    figure.show()
