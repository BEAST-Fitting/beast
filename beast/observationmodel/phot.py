"""
Photometric package
===================

Defines a Filter class and associated functions to extract photometry.

This also include functions to keep libraries up to date

.. note::

    integrations are done using :func:`trapz`
    Why not Simpsons? Simpsons principle is to take sequence of 3 points to
    make a quadratic interpolation. Which in the end, when filters have sharp
    edges, the error due to this "interpolation" are extremely large in
    comparison to the uncertainties induced by trapeze integration.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import numpy

import tables
from scipy.integrate import trapz

from ..tools.decorators import timeit
from ..config import __ROOT__

__default__      = __ROOT__ + '/filters.hd5'
__default_vega__ = __ROOT__ + '/vega.hd5'

# this is used to convert from bolometric luminosities to abs fluxes
# object to 10parsecs -- abs mag.
distc = 4. * numpy.pi * (3.0856775e19) ** 2

__all__ = ['Filter', 'IntegrationFilter', 'load_all_filters', 'load_filters',
           'load_Integrationfilters', 'extractPhotometry', 'extractSEDs',
           'STmag_to_flux', 'STmag_from_flux', 'fluxToMag', 'fluxErrTomag',
           'magToFlux','magErrToFlux', 'append_filter', 'appendVegaFilter']

class Filter(object):
    """Class filter
    Define a filter by its name, wavelength and transmission
    """
    #----------------------------------------------------------------------
    def info(self):

        print("""
Filter object information:
    name: {s.name:s}
    central wavelength: {s.cl:f}
    norm: {s.norm:f}
    pivot wavelength: {s.lpivot:f}
    definition contains {s.transmit.size:d} points""".format(s=self))
        """ display information about the current filter"""

    def __repr__(self):
        return "Filter: %s, %s" % (self.name, object.__repr__(self))

    def getFlux(self, slamb, sflux):
        """getFlux
        Integrate the flux within the filter and return the integrated energy
        If you consider applying the filter to many spectra, you might want to
        consider extractSEDs.

        Parameters
        ----------
        slamb: ndarray(dtype=float, ndim=1)
            spectrum wavelength definition domain

        sflux: ndarray(dtype=float, ndim=1)
            associated flux

        Returns
        -------
        flux: float
            Energy of the spectrum within the filter
        """
        if True in numpy.isinf(sflux):
            indinf = numpy.where(numpy.isinf(sflux))
            indfin = numpy.where(numpy.isfinite(sflux))
            sflux[indinf] = numpy.interp(slamb[indinf], slamb[indfin], sflux[indfin])
        ifT = numpy.interp(slamb, self.wavelength, self.transmit, left=0., right=0.)
        if True in (ifT > 0.):
            ind = numpy.where(ifT > 0.)
            a = numpy.trapz( slamb[ind] * ifT[ind] * sflux[ind], slamb[ind] )
            b = numpy.trapz( slamb[ind] * ifT[ind], slamb[ind] )
            if (numpy.isinf(a) | numpy.isinf(b)):
                print(self.name, "Warn for inf value")
            return a / b
        else:
            return 0.

    def __call__(self, slamb, sflux):
        return self.applyTo(slamb, sflux)

    def applyTo(self, slamb, sflux):
        """
        Apply filter to a spectrum

        Parameters
        ----------
        slamb: ndarray
            spectrum wavelength definition domain

        sflux: ndarray
            associated flux

        Returns
        -------
        flux: float
            new spectrum values accounting for the filter
        """
        ifT = numpy.interp(slamb, self.wavelength, self.transmit)
        return ifT * sflux

    def __init__(self, wavelength, transmit, name=''):
        """Constructor"""
        self.name       = name
        self.wavelength = wavelength
        self.transmit   = transmit
        self.norm       = trapz(transmit, wavelength)
        self.lT         = trapz(wavelength * transmit, wavelength)
        self.lpivot     = numpy.sqrt( self.lT / trapz(transmit / wavelength, wavelength) )
        self.cl         = self.lT / self.norm


class IntegrationFilter(object):
    """Class filter

    Define an integration filter from the range of integration
    """
    def info(self):

        print("""
    Integration Filter object information:
    name: {s.name:s}
    central wavelength: {s.cl:f}
    norm: {s.norm:f}
    pivot wavelength: {s.lpivot:f}
    definition contains {s.transmit.size:d} points""".format(s=self))
        """ display information about the current filter"""


    def __repr__(self):
        return "Filter: %s, %s" % (self.name, object.__repr__(self))

    def getFlux(self, slamb, sflux):
        """getFlux
        Integrate the flux within the filter and return the integrated energy
        If you consider applying the filter to many spectra, you might want to
        consider extractSEDs.

        Parameters
        ----------
        slamb: ndarray(dtype=float, ndim=1)
            spectrum wavelength definition domain

        sflux: ndarray(dtype=float, ndim=1)
            associated flux

        Returns
        -------
        flux: float
            Energy of the spectrum within the filter
        """
        if True in numpy.isinf(sflux):
            indinf = numpy.where(numpy.isinf(sflux))
            indfin = numpy.where(numpy.isfinite(sflux))
            sflux[indinf] = numpy.interp(slamb[indinf], slamb[indfin], sflux[indfin])

        # find common wavelength interval
        ind = ((slamb <= self.wavelength.max()) & (slamb >= self.wavelength.min()))

        if True in ind:
            _slamb = slamb[ind]
            a = numpy.trapz(_slamb * sflux[ind], _slamb)
            b = numpy.trapz(numpy.ones(_slamb.shape, dtype=float) * _slamb, _slamb)

            if (numpy.isinf(a) | numpy.isinf(b)):
                print(self.name, "Warn for inf value")
            return a / b
        else:
            return 0.

    def applyTo(self, slamb, sflux):
        """
        Apply filter to a spectrum

        Parameters
        ----------
        slamb: ndarray
            spectrum wavelength definition domain

        sflux: ndarray
            associated flux

        Returns
        -------
        flux: float
            new spectrum values accounting for the filter
        """
        # find common wavelength interval
        ind = ((slamb <= self.wavelength.max()) & (slamb >= self.wavelength.min()))

        r = numpy.sflux[:]
        r[~ind] = 0
        return r

    def __init__(self, wavelength, transmit, name=''):
        """Constructor"""
        self.name       = name
        self.wavelength = wavelength
        self.transmit   = transmit
        self.norm       = trapz(transmit, wavelength)
        self.lT         = trapz(transmit * wavelength, wavelength)
        self.lpivot     = numpy.sqrt( self.lT / trapz(1. / wavelength, wavelength) )
        self.cl         = self.lT / self.norm




def __load__(fname, ftab, interp=True, lamb=None, integrationFilter=False):
    """ Load a given filter from the library

    Parameters
    ----------
    fname: str
        normalized names according to filtersLib

    ftab: hd5root
        root from the filter library hd5 file

    interp: bool, optional
        reinterpolate the filters over given lambda points

    lamb: ndarray[float, ndim=1]
        desired wavelength definition of the filter

    integrationFilter: bool, optional
        set True for specail integraion filter such as Qion or E_uv
        if set, lamb should be given

    Returns
    -------
    filter: Filter instance
        filter object
    """
    if integrationFilter:
        ifT = numpy.interp(lamb, fname.wavelength, fname.transmit, left=0., right=0.)
        return IntegrationFilter(lamb, ifT, name=fname)

    else:
        fnode    = ftab.get_node('/filters/' + fname)
        flamb    = fnode[:]['WAVELENGTH']
        transmit = fnode[:]['THROUGHPUT']
        if interp & (lamb is not None):
            ifT = numpy.interp(lamb, flamb, transmit, left=0., right=0.)
            return Filter( lamb, ifT, name=fnode.name )
        else:
            return Filter( flamb, transmit, name=fnode.name )



def load_all_filters(interp=True, lamb=None, filterLib=None):
    """ load all filters from the library

    Parameters
    ----------
    interp: bool
        reinterpolate the filters over given lambda points

    lamb: ndarray[float, ndim=1]
        desired wavelength definition of the filter

    filterLib:  str
        path to the filter library hd5 file

    Returns
    -------
    filters: list[filter]
        list of filter objects
    """
    if filterLib is None:
        filterLib = __default__
    with tables.open_file(filterLib, 'r') as ftab:
        filters = [ __load__(fname, ftab, interp=interp, lamb=lamb) for fname in ftab.root.content.cols.TABLENAME ]
    return(filters)


def load_filters(names, interp=True, lamb=None, filterLib=None):
    """ load a limited set of filters

        Parameters
        ----------
        names: list[str]
            normalized names according to filtersLib

        interp: bool
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter

        filterLib: path
            path to the filter library hd5 file

        Returns
        -------
        filters: list[filter]
            list of filter objects
    """
    if filterLib is None:
        filterLib = __default__
    with tables.open_file(filterLib, 'r') as ftab:
        filters = [ __load__(fname, ftab, interp=interp, lamb=lamb) for fname in names ]
    return(filters)


def load_Integrationfilters(flist, interp=True, lamb=None):
    """ load a limited set of filters

        Parameters
        ----------
        flist: sequence(filter)
            list of filter object instances

        interp: bool
            reinterpolate the filters over given lambda points

        lamb: ndarray[float, ndim=1]
            desired wavelength definition of the filter


        Returns
        -------
        filters: list[filter]
            list of filter objects
    """
    filters = [ __load__(fname, ftab=None, interp=interp, lamb=lamb, integrationFilter=True) for fname in flist ]
    return(filters)



def extractPhotometry(lamb, spec, flist, absFlux=True):
    """Extract seds from a one single spectrum

    Parameters
    ----------
    lamb: ndarray[float,ndim=1]
        wavelength of spec

    spec: ndarray[float, ndim=1]
        spectrum

    flist: list[filter]
        list of filter objects

    absflux: bool
        return SEDs in absolute fluxes if set

    Returns
    -------
    cls: ndarray[float, ndim=1]
        filters central wavelength

    seds: ndarray[float, ndim=1]
        integrated sed
    """
    cls  = numpy.empty( len(flist), dtype=float)
    seds = numpy.empty( len(flist), dtype=float)
    for e, k in enumerate(flist):
        xl  = k.transmit > 0.
        tmp = lamb[xl] * k.transmit[xl]
        s0  = spec[:, xl]
        # apply absolute flux conversion if requested
        if absFlux:
            s0 /= distc
        a = trapz( tmp[None, :] * s0, lamb[xl], axis=1 )
        seds[e] = a / k.lT   # divide by integral (lambda T dlambda)
        cls[e]  = k.cl

    return cls, seds


def extractSEDs(g0, flist, absFlux=True):
    """ Extract seds from a grid

    Parameters
    ----------
    g0: ModelGrid instance
        initial spectral grid

    flist: sequence(filter)
        list of filter object instances

    absflux: bool
        return SEDs in absolute fluxes if set

    Returns
    -------
    cls: ndarray[float, ndim=1]
        filters central wavelength

    seds: ndarray[float, ndim=1]
        integrated sed

    grid: Table
        SED grid properties table from g0 (g0.grid)
    """
    lamb = g0.lamb
    seds = numpy.empty(( len(g0.grid), len(flist) ), dtype=float)
    cls  = numpy.empty( len(flist), dtype=float)
    for e, k in enumerate(flist):
        xl  = k.transmit > 0.
        tmp = lamb[xl] * k.transmit[xl]
        s0  = g0.seds[:, xl]
        # apply absolute flux conversion if requested
        if absFlux:
            s0 /= distc
        a = trapz( tmp[None, :] * s0, lamb[xl], axis=1 )
        seds[:, e] = a / k.lT
        cls[e] = k.cl

    #memgrid = grid.MemoryGrid(cls, seds, g0.grid)
    #return memgrid
    return cls, seds, g0.grid


def STmag_to_flux( v ):
    """
    Convert an ST magnitude to erg/s/cm2/AA (Flambda)

        mag = -2.5 \log_{10}(F) - 21.10

        M0 = 21.10
        F0 = 3.6307805477010028 10^{-9} erg/s/cm2/AA

    Parameters
    ----------
    v: np.ndarray[float, ndim=N] or float
        array of magnitudes

    Returns
    -------
    flux: np.ndarray[float, ndim=N], or float
        array of fluxes
    """
    v0 = 21.1
    return 10. ** ( -0.4 * (v - v0) )


def STmag_from_flux( v ):
    """
    Convert to ST magnitude from erg/s/cm2/AA (Flambda)

        mag = -2.5 \log_{10}(F) - 21.10

        M0 = 21.10
        F0 = 3.6307805477010028 10^{-9} erg/s/cm2/AA

    Parameters
    ----------
    v: np.ndarray[float, ndim=N], or float
        array of fluxes

    Returns
    -------
    mag: np.ndarray[float, ndim=N], or float
        array of magnitudes
    """
    v0 = 21.1
    return -2.5 * numpy.log10( v ) - v0


def fluxToMag(flux):
    """ Return the magnitudes from flux values

    Parameters
    ----------
    flux: np.ndarray[float, ndim=N]
        array of fluxes

    Returns
    -------
    mag: np.ndarray[float, ndim=N]
        array of magnitudes
    """
    return -2.5 * numpy.log10(flux)


def fluxErrTomag(flux, fluxerr):
    """ Return the magnitudes and associated errors from fluxes and flux error
    values

    Parameters
    ----------
    flux:    np.ndarray[float, ndim=1]
        array of fluxes

    fluxerr: np.ndarray[float, ndim=1]
        array of flux errors

    Returns
    -------
    mag: np.ndarray[float, ndim=1]
        array of magnitudes

    err: np.ndarray[float, ndim=1]
        array of magnitude errors
    """
    mag = fluxToMag(flux)
    return mag, -2.5 * numpy.log10( 1. - fluxerr / flux )


def magToFlux(mag):
    """ Return the flux from magnitude values

    Parameters
    ----------
    mag: np.ndarray[float, ndim=N]
        array of magnitudes

    Returns
    -------
    flux:  np.ndarray[float, ndim=N]
        array of fluxes
    """
    return 10 ** (-0.4 * mag)


def magErrToFlux(mag, err):
    """ Return the flux and associated errors from magnitude and mag error values

    Parameters
    ----------
    mag: np.ndarray[float, ndim=1]
        array of magnitudes

    err: np.ndarray[float, ndim=1]
        array of magnitude errors

    Returns
    -------
    flux:    np.ndarray[float, ndim=1]
        array of fluxes

    fluxerr: np.ndarray[float, ndim=1]
        array of flux errors
    """
    flux = magToFlux(mag)
    return flux, flux * ( 1. - magToFlux(err) )


class __newFilterTable__(tables.IsDescription):
    """ define table to store filter dataset """
    WAVELENGTH = tables.FloatCol(pos=0)
    THROUGHPUT = tables.FloatCol(pos=1)


def append_filter(lamb, flux, tablename, observatory, instrument, name,
                  comment=None, filterLib=__default__, updateVegaLib=True):
    """
    Edit the filter catalog and append a new one given by its transfer function

    Parameters
    ----------
    lamb: ndarray(dtype=float)
        wavelength of the filter definition

    flux: ndarray(dtype=float)
        transimission of the filter

    tablename: str
        table name in the library

    observatory: str
        observatory of the filter (Ground, HST, Spitzer, ...)

    instrument: str
        instrument associated with the filter

    name: str
        name of the filter

    comment: str, optional
        optinal comment to keep with the filter

    filterLib: str, optional
        filter library file to use

    updateVegaLib: bool
        if set calls the update function to the vega library
    """
    ftab = tables.open_file(filterLib, 'a')
    contentTab = ftab.get_node('/content')
    if contentTab.read_where('TABLENAME == "{0}"'.format(tablename)).size > 0:
        print('% {0}: Filter {1} already exists. Returning'.format(sys.argv[0], tablename))
        return

    # Gen Filter object including relevant details
    filtInst = list(filter(lamb, flux, name=name))
    # Add a new line in the content table
    newRow = contentTab.row
    newRow['TABLENAME'] = tablename
    newRow['OBSERVATORY'] = observatory
    newRow['INSTRUMENT'] = instrument
    newRow['NAME'] = filtInst.name
    newRow['NORM'] = filtInst.norm
    newRow['CWAVE'] = filtInst.cl
    newRow['PWAVE'] = filtInst.lpivot
    if comment is not None:
        newRow['COMMENT'] = comment
    newRow.append()
    contentTab.flush()
    # Create Table
    newTab = ftab.create_table('/filters', tablename, __newFilterTable__,
                              title=filtInst.name,
                              expectedrows=filtInst.wavelength.size)
    newRow = newTab.row
    for i in range(filtInst.wavelength.size):
        newRow["WAVELENGTH"] = filtInst.wavelength[i]
        newRow["THROUGHPUT"] = filtInst.transmit[i]
        newRow.append()
    newTab.flush()
    ftab.flush()
    ftab.close()
    print('% {0}: Filter {1} added to {2}'.format(sys.argv[0], name, filterLib))
    if updateVegaLib:
        appendVegaFilter(filtInst)


#-------------------------------------------------------------------------------
# VEGA SPECTRUM and VEGA ZEROPOINTS
#-------------------------------------------------------------------------------
def __analyseVegaSpectrum__(w, f, filters):
    """
    Returns property information from the application of a given set of filters

    Parameters
    ----------
    w: ndarray(dtype=float)
        wavelength definition of the spectrum

    f: ndarray(dtype=float)
        flux definition of the spectrum

    Returns
    -------
    props: dict
        properties to store with the filter that includes flux, magnitude
        values.
    """
    nFilters = len(filters)
    phot     = numpy.zeros((nFilters))
    cwave    = numpy.zeros((nFilters))
    fname    = []
    mag      = numpy.zeros((nFilters))
    for j in range(0, nFilters):
        fname.append(filters[j].name)
        cwave[j] = filters[j].cl
        phot[j] = filters[j].getFlux( w, f )
        mag[j] = -2.5 * numpy.log10(phot[j])
    return ({ 'fname': fname, 'cwave': cwave, 'lum': phot, 'mag': mag})


def appendVegaFilter(filtInst, VegaLib=__default_vega__):
    """
    Add filter properties to the Vega library

    Parameters
    ----------
    filtInst: Filter instance
        filter instance to get properties from and store information with Vega.

    VegaLib: str
        Vega Library
    """
#    import tables
    vtab = tables.open_file(VegaLib, 'a')
    vl = vtab.root.spectrum[:]['WAVELENGTH']
    vf = vtab.root.spectrum[:]['FLUX']
    sedTab = vtab.get_node('/sed')
    if sedTab.read_where('FNAME == "{0}"'.format(filtInst.name)).size > 0:
        print('% {0}: Filter {1} already exists. Returning'.format(sys.argv[0], filtInst.name))
        return

    data = __analyseVegaSpectrum__(vl, vf, [filtInst])
    newRow = sedTab.row
    newRow['FNAME'] = filtInst.name
    newRow['CWAVE'] = filtInst.cl
    newRow['LUM']   = data['lum'][0]
    newRow['MAG']   = data['mag'][0]
    newRow.append()
    sedTab.flush()
    vtab.close()
    print('% {0}: Filter {1} added to {2}'.format(sys.argv[0], filtInst.name, VegaLib))
