""" Unit test for SEDfitter """
__version__ = '0.2dev'

import numpy
from numpy import exp
import inspect
import itertools
import mytables

from anased import computeLogLikelihood
from decorators import timeit
import extinction
import grid
import phot
import observations


#TODO: check the Alambda definition and update the calls is necessary
#TODO: cleanup and put meshgrid, iter_Av_grid into the grid package or a common package

def meshgrid(arrs, ravel=False):
    """ Generate a n-dim grid from a list of n arrays containing the points
    to use. The gridding occurs by varying the last value first.

    INPUTS:
        args    list[ndarray]   a,b,c containing the grid points of each dimension
    OUTPUTS:
        ans grid/array          output comparable to numpy.meshgrid or ndarray
    KEYWORDS:
        ravel   bool            return a 2d array if set, where each line is a n-d vector point
    """
    arrs = tuple(arrs)
    lens = map(len, arrs)
    dim = len(arrs)
    sz = 1
    for s in lens:
        sz *= s
    ans = []
    for i, arr in enumerate(arrs):
        slc = [1] * dim
        slc[i] = lens[i]
        arr2 = numpy.asarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j != i:
                arr2 = arr2.repeat(sz, axis=j)
        ans.append(arr2)
    if ravel:
        return numpy.vstack( map( numpy.ravel, tuple(ans) ) ).T
    else:
        return tuple(ans)


def iter_Av_grid(g0, oAv, **kwargs):
    """ generate an iterator that will create grids
    by applying extinction values on the spectra
    This does not include the integrated photometry step!

    INPUTS:
        g0  grid            initial spectral grid (unreddened)
        oAv ExtinctionLaw   extinction law to apply

    OUTPUTS:
        git iterator    iterator on SpectralGrid chunks

    KEYWORDS:
        **kwargs    grid point values.
                    each keyword will be used to call oAv.function
    """

    assert( isinstance(oAv, extinction.ExtinctionLaw)), 'Extinction must be an ExtinctionLaw instance, got %s' % type(oAv)
    assert( isinstance(g0,  grid.ModelGrid) ), 'Extinction must be an ExtinctionLaw instance, got %s' % type(oAv)

    #checking arguments of extinction.function
    av_args = [k for k in inspect.getargspec(extinction.RvFbumpLaw.function).args if k not in ['self', 'lamb', 'Alambda'] ]

    for k, v in kwargs.iteritems():
        assert( k in av_args ), 'Argument %s not in the extinction parameters' % k

    #generate the grid points
    if len(kwargs) > 1:
        gpts = meshgrid( kwargs.values(), ravel=True )
    else:
        gpts = kwargs.values()[0]

    #prepare output

    def gensubgrid( theta_av ):
        """
        **Partial function**
        Generate the chunk of grid corresponding to a given theta_av
        based on the initial grid g0
        """
        _args = {}
        for e, k in enumerate(kwargs):
            if hasattr(theta_av, '__iter__'):
                _args[k] = theta_av[e]
            else:
                _args[k] = theta_av
        tau        = oAv.function(g0.lamb * 1e-4, Alambda=True, **_args)
        outputSEDs = g0.seds * numpy.exp(-tau)[None, :]
        #copy original grid
        t = mytables.Table(g0.grid.data, header=g0.grid.header, name=g0.grid.header['NAME'] )
        for k, v in _args.iteritems():
            t.addCol( [ v ] * t.nrows, name=k)
        g = grid.SpectralGrid()
        g.lamb = g0.lamb[:]
        g.seds = outputSEDs
        g.grid = t
        return g

    # If the grid is small enough, you can compute all at once and
    # concatenate the result. However, most of the time the grid is big
    # enough that the maximum of element per array is reached (~10^9 on
    # 64bits arch). So that the default is to return an iterator on the grid
    # chunks
    # This means each chunk is independent from the others as well, ergo
    # parallel jobs are possible.
    return itertools.imap( gensubgrid, gpts )


def getFakeStar(g, idx, filts, err=0.1, oAv=None, **kwargs):
    """ Generate a fake sed from a model grid
    INPUTS:
        idx int                         index number on the grid
        filters list[filter]            list of filter names
    OUTPUTS:
        fakein  int                     the index of the model on the grid
        lamb    ndarray[float, ndim=1]  wavelength taken from the grid
        fakesed ndarray[float, ndim=1]  resulting SED

    KEYWORDS:
        err float                       proportional error to consider on the fluxes
        oAv ExtinctionLaw               if provided, apply extinction function using **kwargs

        **kwargs if provided, extra keywords are used to apply an extinction
    """
    lamb     = g.lamb
    fakein   = idx
    fakesed  = numpy.copy(g.seds[fakein, :])
    if (oAv is not None) & (len(kwargs) > 0):
        tau      = oAv(g.lamb * 1e-4, Alambda=True, **kwargs)
        fakesed *= exp(-tau)
    ## extract photometry
    fakecl, fakesed = phot.extractPhotometry(lamb, fakesed, filts, absFlux=True)

    #magerr  = 0.05
    #fakeerr = fakesed * (1. - 10**(-0.4*magerr) )
    fakeerr = err * fakesed

    return fakein, fakecl, fakesed, fakeerr


def getFakeCluster(g, age, Z, filts, err=0.1, suffix='0', Nstars=None, oAv=None, **kwargs):
    """ Generate a fake sed from a model grid
    INPUTS:
        idx     int                     index number on the grid
        age     float                   age of the desired population
        Z       float                   metallicity of the desired population
        filters list[filter]            list of filter names
    OUTPUTS:
        fakein  int                     the index of the model on the grid
        lamb    ndarray[float, ndim=1]  wavelength taken from the grid
        fakesed ndarray[float, ndim=1]  resulting SED

    KEYWORDS:
        suffix  str                 suffix used to store the cluster table
        err     float               proportional error to consider on the fluxes
        Nstars  int                 number of stars in the population (None = all possible stars)
        oAv     ExtinctionLaw       if provided, apply extinction function using **kwargs

        **kwargs if provided, extra keywords are used to apply an extinction
    """
    # find all the stars matching age and Z criteria
    glogA  = numpy.unique(g.logA)
    gZ     = numpy.unique(g.Z)
    _logA  = glogA[numpy.argmin( abs(glogA - numpy.log10(age)) )]
    _Z     = gZ[numpy.argmin( abs(gZ - Z) )]
    idx    = numpy.where((g.logA == _logA) & (g.Z == _Z))[0]
    # Restrict to Nstars is provided
    if Nstars is not None:
        idx = numpy.random.randint(0, len(idx), min([len(idx), Nstars]) )

    fakeseds = g.seds[idx]
    # Compute extinction is requested
    if (oAv is not None) & (len(kwargs) > 0):
        tau       = oAv(g.lamb * 1e-4, Alambda=True, **kwargs)
        fakeseds *= exp(-tau)[None, :]
    ## extract photometry
    gg      = grid.SpectralGrid()
    gg.lamb = g.lamb
    gg.seds = fakeseds
    gg.grid = g.grid[idx]
    seds    = gg.getSEDs(filts, absFlux=True)

    seds.grid.addCol(idx, name='idx')

    for k, v in kwargs.iteritems():
        seds.grid.header[k] = v

    filter_names = []
    for ek, nk in enumerate(filts):
        fname = nk.name
        seds.grid.addCol(seds.seds[:, ek], name=fname)
        seds.grid.addCol(seds.seds[:, ek] * err, name=fname + 'err')
        filter_names.append(nk.name)
    seds.grid.header['logA'] = _logA
    seds.grid.header['Z'] = _Z

    seds.grid.write('Tests/cl_%s.fits' % suffix)
    del seds

    obs = observations.Observations('Tests/cl_%s.fits' % suffix)
    obs.setFilters(filter_names)
    return obs


def test_observations_seds(err=0.1, age=1e7, Z=0.02, suffix='0'):
    filter_names  = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    #Load the initial model grid
    g0 = grid.FileSpectralGrid('libs/stellib_kurucz2004_padovaiso.spectralgrid.fits')

    #load the filters
    filts = phot.load_filters(filter_names, interp=True, lamb=g0.lamb)
    cls   = numpy.array( [ k.cl for k in filts ] )

    # define the extinction priors
    oAv = extinction.Cardelli()

    ##define Av with a step size
    #Av = numpy.arange(0.,3., 0.1)

    ##define Av with a number of points
    Av = numpy.linspace(0, 1, 3)
    #Av = numpy.array([0.0], dtype=float)
    Rv = numpy.array([3.1], dtype=float)
    fb = numpy.array([0.5], dtype=float)     # only one value still needs to be array/list/tuple

    ## Extinction parameters are randomly drawn from the extinction space
    Av0    = Av[numpy.random.randint(0, len(Av), 1)]
    Rv0    = Rv[numpy.random.randint(0, len(Rv), 1)]
    fb0    = fb[numpy.random.randint(0, len(fb), 1)]

    with timeit( 'Cluster generation: age=%g, Z=%0.4f' % (age, Z) ):
        obs = getFakeCluster(g0, age, Z, filts, err=0.1,
                            suffix=suffix, Nstars=None,
                            oAv=oAv, Av=float(Av0), Rv=float(Rv0), f_bump=float(fb0))

    #get the grid iterator
    iter_grid = iter_Av_grid(g0, oAv, Av=Av, Rv=Rv, f_bump=fb)

    for tn, (fakesed, fakeerr, mask) in obs.enumobs():
        #define output filename
        outname = 'Tests/cl_%s_t%d.fits' % ( suffix, tn )
        with timeit('Likelihood Object %d' % tn):
            for gk in iter_grid:
                with timeit('\t * Computing SEDS'):
                    seds = gk.getSEDs(filts, absFlux=True)
                with timeit('\t * Computing Lnp'):
                    lnp = computeLogLikelihood(fakesed, fakeerr, seds.seds, normed=False, mask=mask)
                with timeit('\t * Writing outputs'):
                    t = mytables.Table(name='SEDOUT')
                    #t = gk.grid
                    t.addCol(numpy.arange(seds.grid.nrows, dtype=int), name='idx')
                    t.addCol(lnp, name='lnp')
                    for ek, nk in enumerate(filter_names):
                        t.addCol(seds.seds[:, ek], name=nk)
                    t.header['Av'] = Av0[tn]
                    t.header['Rv'] = Rv0[tn]
                    t.header['f_bump'] = fb0[tn]
                    t.write( outname )

        with timeit('\t * Writing Original data'):
            d = dict(filters=filter_names, flux=fakesed, err=fakeerr, cl=cls)
            t = mytables.Table(d, name='INPUT')
            t.header['Av']     = Av0[tn]
            t.header['Rv']     = Rv0[tn]
            t.header['f_bump'] = fb0[tn]
            t.header['idx']    = tn
            t.write( outname )

        with timeit('\t * Writing Original grid'):
            g0.grid.header['NAME'] = 'MODELS'
            g0.grid.header['EXTNAME'] = 'MODELS'
            g0.grid.write( outname )

    return locals()


def test_seds(err=0.1):
    filter_names  = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    #Load the initial model grid
    g0 = grid.FileSpectralGrid('libs/stellib_kurucz2004_padovaiso.spectralgrid.fits')

    #load the filters
    filts = phot.load_filters(filter_names, interp=True, lamb=g0.lamb)

    # define the extinction priors
    oAv = extinction.Cardelli()

    ##define Av with a step size
    #Av = numpy.arange(0.,3., 0.1)

    ##define Av with a number of points
    Av = numpy.linspace(0, 1, 3)
    #Av = numpy.array([0.0], dtype=float)
    Rv = numpy.array([3.1], dtype=float)
    fb = numpy.array([0.5], dtype=float)     # only one value still needs to be array/list/tuple

    #get the grid iterator
    iter_grid = iter_Av_grid(g0, oAv, Av=Av, Rv=Rv, f_bump=fb)

    # Parameters from which we generate fake data
    ## Number of fake SEDs to play with
    N   = 1
    ## Extinction parameters are randomly drawn from the extinction space
    Av0    = Av[numpy.random.randint(0, len(Av), N)]
    Rv0    = Rv[numpy.random.randint(0, len(Rv), N)]
    fb0    = fb[numpy.random.randint(0, len(fb), N)]
    ## Initial SEDs are randomly drawn from the model space
    fakein = numpy.random.randint(0, g0.grid.nrows, N)

    for tn in range(N):
        #fake DATA
        idx, fakecl, fakesed, fakeerr = getFakeStar(g0, fakein[tn], filts, err=err, oAv=oAv, Av=Av0[tn], Rv=Rv0[tn], f_bump=fb0[tn])
        mask = numpy.zeros(fakesed.shape, dtype=bool)
        ## simulate non detection
        #mask[3] = True
        #mask[2] = True

        with timeit('Likelihood Object %d' % tn):
            #define output filename
            outname = 'Tests/t%d.fits' % ( tn )
            for gk in iter_grid:
                with timeit('\t * Computing SEDS'):
                    seds = gk.getSEDs(filts, absFlux=True)
                with timeit('\t * Computing Lnp'):
                    lnp = computeLogLikelihood(fakesed, fakeerr, seds.seds, normed=False, mask=mask)
                with timeit('\t * Writing outputs'):
                    t = mytables.Table(name='SEDOUT')
                    #t = gk.grid
                    t.addCol(numpy.arange(seds.grid.nrows, dtype=int), name='idx')
                    t.addCol(lnp, name='lnp')
                    for ek, nk in enumerate(filter_names):
                        t.addCol(seds.seds[:, ek], name=nk)
                    t.header['Av'] = Av0[tn]
                    t.header['Rv'] = Rv0[tn]
                    t.header['f_bump'] = fb0[tn]
                    t.write( outname )

        with timeit('\t * Writing Original data'):
            d = dict(filters=filter_names, flux=fakesed, err=fakeerr, cl=fakecl)
            t = mytables.Table(d, name='INPUT')
            t.header['Av']     = Av0[tn]
            t.header['Rv']     = Rv0[tn]
            t.header['f_bump'] = fb0[tn]
            t.header['idx']    = idx
            t.write( outname )

        with timeit('\t * Writing Original grid'):
            g0.grid.header['NAME'] = 'MODELS'
            g0.grid.header['EXTNAME'] = 'MODELS'
            g0.grid.write( outname )

    return locals()


def task(obs, idx, suffix, gridfile, oAv, Av, Rv, fb):
    #get the grid iterator
    g0 = grid.FileSpectralGrid(gridfile)
    iter_grid = iter_Av_grid(g0, oAv, Av=Av, Rv=Rv, f_bump=fb)

    #load the filters
    filts = phot.load_filters(obs.filters, interp=True, lamb=g0.lamb)
    cls   = numpy.array( [ k.cl for k in filts ] )

    for tn in idx:
        fakesed, fakeerr, mask = obs.getObs(idx)
        #define output filename
        outname = 'Tests/cl_%s_t%d.fits' % ( suffix, tn )
        with timeit('Likelihood Object %d' % tn):
            for gk in iter_grid:
                with timeit('\t * Computing SEDS'):
                    seds = gk.getSEDs(filts, absFlux=True)
                with timeit('\t * Computing Lnp'):
                    lnp = computeLogLikelihood(fakesed, fakeerr, seds.seds, normed=False, mask=mask)
                with timeit('\t * Writing outputs'):
                    t = mytables.Table(name='SEDOUT')
                    #t = gk.grid
                    t.addCol(numpy.arange(seds.grid.nrows, dtype=int), name='idx')
                    t.addCol(lnp, name='lnp')
                    for ek, nk in enumerate(obs.filters):
                        t.addCol(seds.seds[:, ek], name=nk)
                    t.write( outname )

        with timeit('\t * Writing Original data'):
            d = dict(filters=obs.filters, flux=fakesed, err=fakeerr, cl=cls)
            t = mytables.Table(d, name='INPUT')
            t.header['idx']    = tn
            t.write( outname )

        with timeit('\t * Writing Original grid'):
            g0.grid.header['NAME'] = 'MODELS'
            g0.grid.header['EXTNAME'] = 'MODELS'
            g0.grid.write( outname )


def test_parallel_observations_seds(err=0.1, age=1e7, Z=0.02, suffix='0', nthreads=1, pool=None):
    import multiprocessing as mp
    if pool is None:
        __pool__ = mp.Pool(nthreads)
    else:
        if isinstance(pool, mp.Pool):
            nthreads = getattr(pool, '_processes', None)

    filter_names  = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    #Load the initial model grid
    gridfile = 'libs/stellib_kurucz2004_padovaiso.spectralgrid.fits'
    g0 = grid.FileSpectralGrid(gridfile)

    #load the filters
    filts = phot.load_filters(filter_names, interp=True, lamb=g0.lamb)

    # define the extinction priors
    oAv = extinction.Cardelli()

    ##define Av with a step size
    #Av = numpy.arange(0.,3., 0.1)

    ##define Av with a number of points
    Av = numpy.linspace(0, 1, 3)
    #Av = numpy.array([0.0], dtype=float)
    Rv = numpy.array([3.1], dtype=float)
    fb = numpy.array([0.5], dtype=float)     # only one value still needs to be array/list/tuple

    ## Extinction parameters are randomly drawn from the extinction space
    Av0    = Av[numpy.random.randint(0, len(Av), 1)]
    Rv0    = Rv[numpy.random.randint(0, len(Rv), 1)]
    fb0    = fb[numpy.random.randint(0, len(fb), 1)]

    with timeit( 'Cluster generation: age=%g, Z=%0.4f' % (age, Z) ):
        obs = getFakeCluster(g0, age, Z, filts, err=0.1,
                            suffix=suffix, Nstars=None,
                            oAv=oAv, Av=float(Av0), Rv=float(Rv0), f_bump=float(fb0))

    def partition(list, n):
        """Partition list into n nearly equal length sublists"""
        division = len(list) / float(n)
        return [list[int(round(division * i)): int(round(division * (i + 1)))] for i in range(n)]


    #Scatter the work
    __pool__.map( task,  partition(range(len(obs)), nthreads) )

    return locals()
