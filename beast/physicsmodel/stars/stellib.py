"""
Stellib class

Intent to implement a generic module to manage stellar library from various
sources.

The interpolation is implemented from the pegase.2 fortran converted algorithm.
(this may not be pythonic though)
"""
from __future__ import (absolute_import, division, print_function)

import numpy as np
from scipy.interpolate import interp1d
from numpy.lib import recfunctions
from astropy import constants, units

from ..grid import SpectralGrid
from ..helpers import Path

from ...external.eztables import Table
from ...config import __ROOT__, __NTHREADS__
from .include import __interp__
from ...tools.pbar import Pbar
from ...tools.helpers import nbytes

#lsun = 3.839e+26   # in W (Watts)
lsun = constants.L_sun.value
#sig_stefan = 5.67037321 * 1e-8  # W * m**-2 * K**-4
sig_stefan = constants.sigma_sb.value
#rsun = 6.955e8  # in meters
rsun = constants.R_sun.value

config = {
    'basel_2.2_pegase': __ROOT__ + '/BaSeL_v2.2.pegase.grid.fits',
    'elodie_3.1': __ROOT__ + '/Elodie_v3.1.grid.fits',
    'kurucz': __ROOT__ + '/kurucz2004.grid.fits',
    'tlusty': __ROOT__ + '/tlusty.lowres.grid.fits',
    'btsettl': __ROOT__ + '/bt-settl.lowres.grid.fits',
    'btsettl_medres': __ROOT__ + '/bt-settl.medres.grid.fits',
    'munari': __ROOT__ + '/atlas9-munari.hires.grid.fits'
}

__all__ = ['Stellib', 'CompositeStellib', 'Kurucz', 'Tlusty',
           'BTSettl', 'Munari', 'Elodie', 'BaSeL']


def isNestedInstance(obj, cl):
    """ Test for sub-classes types
        I could not find a universal test

        keywords
        --------
        obj: object instance
            object to test

        cl: Class
            top level class to test

        returns
        -------
        r: bool
            True if obj is indeed an instance or subclass instance of cl
    """
    tree = []
    for k in cl.__subclasses__():
        tree += k.__subclasses__()
    tree += cl.__subclasses__() + [ cl ]
    return issubclass(obj.__class__, tuple(tree))


def __interpSingle__(args):
    return np.asarray(interp(args[0], args[1], args[2], args[3], args[4],
                             args[5], args[6], args[7], args[8], args[9])).T


def __interpMany__(oSL, logT, logg, Z, logL, dT_max=0.1, eps=1e-06,
                   weights=None, pool=None, nthreads=__NTHREADS__):
    """ run interp on a list of inputs and returns reduced results

    Interpolation of the T,g grid at Z0 metallicity

    Interpolate on the grid and returns star indices and
    associated weights, and Z.
    3 to 12 stars are returned.
    It calls _interp_, but reduce the output to the relevant stars.

    keywords
    --------
    T0  ndarray(float)
        log(Teff) to obtain

    g0  ndarray(float)
        log(g) to obtain

    Z0 ndarray(float)
        metallicity values

    L0 ndarray(float)
        luminosity values

    dT_max: float
        If, T2 (resp. T1) is too far from T compared to T1 (resp. T2),
        i2 (resp. i1) is not used.
        (see below for namings)

    eps: foat
        temperature sensitivity under which points are considered to
        have the same temperature

    weights: ndarray(float)
        luminosity weigths to apply after interpolation

    pool: Pool-like object
        specify a multiprocessing pool for parallel processing

    nthreads: int
        number of processes to use by default

    returns
    -------

    Returns 3 to 12 star indexes and associated weights
    """
    if (pool is None) & (nthreads > 0):
        import multiprocessing as mp
        pool = mp.Pool(nthreads)

    if weights is None:
        seq = [ (logT[k], logg[k], Z[k], logL[k], oSL.Teff, oSL.logg,
                 oSL.Z, dT_max, eps, 1.) for k in range(len(logT)) ]
    else:
        seq = [ (logT[k], logg[k], Z[k], logL[k], oSL.Teff, oSL.logg,
                 oSL.Z, dT_max, eps, weights[k]) for k in range(len(logT)) ]

    if (pool is not None):
        r = pool.map( __interpSingle__, seq )
    else:
        r = list(map( __interpSingle__, seq ))

    return np.vstack(r)


def interp(T0, g0, Z0, L0, T, g, Z, dT_max=0.1, eps=1e-6, weight=1.):
    """ Interpolation of the T,g grid

    Interpolate on the grid and returns star indices and
    associated weights, and Z.
    3 to 12 stars are returned.
    It calls _interp_, but reduce the output to the relevant stars.

    keywords
    --------
    T0  double
        log(Teff) to obtain

    g0  double
        log(g) to obtain

    T   double
        log(Teff) of the grid

    g   double
        log(g) of the grid

    dT_max: float
        If, T2 (resp. T1) is too far from T compared to T1 (resp. T2),
        i2 (resp. i1) is not used.
        (see below for namings)

    eps: foat
        temperature sensitivity under which points are considered to
        have the same temperature

    returns
    -------

    Returns 3 to 12 star indexes and associated weights

    see __interp__

    TODO: compute new weights accounting for Z
    """
    _Z    = Z
    _Zv   = np.unique(_Z)
    _T    = np.asarray(T)
    _g    = np.asarray(g)

    bZ_m  = True in (_Zv == Z0)  # Z_match bool
    r     = np.where((_Zv < Z0))[0]
    Z_inf = _Zv[r.max()] if len(r) > 0 else -1.
    r     = np.where((_Zv > Z0))[0]
    Z_sup = _Zv[r.min()] if len(r) > 0 else -1.

    index   = np.zeros(4 * 3) - 1
    weights = np.zeros(4 * 3)
    Z       = np.zeros(4 * 3)

    if weight is None:
        weight = 1.

    if (bZ_m):
        ind         = np.where(_Z == Z0)
        i, w        = __interp__(T0, g0, _T[ind], _g[ind], dT_max, eps)
        index[8:]   = ind[0][i]
        weights[8:] = w
        Z[8:]       = [Z0] * 4
    else:
        if (Z_inf > 0.):
            ind         = np.where(_Z == Z_inf)
            i, w        = __interp__(T0, g0, _T[ind], _g[ind], dT_max, eps)
            index[:4]   = ind[0][i]
            weights[:4] = w
            Z[:4]       = [Z_inf] * 4

        if (Z_sup > 0.):
            ind          = np.where(_Z == Z_sup)
            i, w         = __interp__(T0, g0, _T[ind], _g[ind], dT_max, eps)
            index[4:8]   = ind[0][i]
            weights[4:8] = w
            Z[4:8]       = [Z_sup] * 4

        if ((Z_inf > 0.) & (Z_sup > 0.)):
            if ( Z_sup - Z_inf ) > 0.:
                fz = (Z0 - Z_inf) / ( Z_sup - Z_inf )
                weights[:4]  *= fz
                weights[4:8] *= ( 1. - fz )
            else:
                weights[:8]  *= 0.5

    ind = np.where(weights > 0)
    return index[ind].astype(int), 10 ** L0 * weight * weights[ind]  # / (weights[ind].sum()) #, Z[ind]


class Stellib(object):
    """ Basic stellar library class """
    def __init__(self, *args, **kargs):
        """ Contructor """
        pass

    def _load_(self):
        raise NotImplementedError

    @property
    def nbytes(self):
        """ return the number of bytes of the object """
        return nbytes(self)

    def interp(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6):
        """ Interpolation of the T,g grid

        Interpolate on the grid and returns star indices and
        associated weights, and Z.
        3 to 12 stars are returned.
        It calls _interp_, but reduce the output to the relevant stars.

        Parameters
        ----------
        T0  double
            log(Teff) to obtain

        g0  double
            log(g) to obtain

        T   double
            log(Teff) of the grid

        g   double
            log(g) of the grid

        dT_max: float
            If, T2 (resp. T1) is too far from T compared to T1 (resp. T2),
            i2 (resp. i1) is not used.
            (see below for namings)

        eps: foat
            temperature sensitivity under which points are considered to
            have the same temperature

        returns
        -------

        Returns 3 to 12 star indexes and associated weights

        see __interp__

        TODO: compute new weights accounting for Z
        """
        _Z    = self.Z
        _T    = np.asarray(self.grid['logT'], dtype=np.double)
        _g    = np.asarray(self.grid['logg'], dtype=np.double)
        return interp(T0, g0, Z0, L0, _T, _g, _Z, dT_max=0.1, eps=1e-6)

    def interpMany(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6, weights=None, pool=None, nthreads=__NTHREADS__):
        """ run interp on a list of inputs and returns reduced results

        Interpolation of the T,g grid at Z0 metallicity

        Interpolate on the grid and returns star indices and
        associated weights, and Z.
        3 to 12 stars are returned.
        It calls _interp_, but reduce the output to the relevant stars.

        Parameters
        ----------
        T0  ndarray(float)
            log(Teff) to obtain

        g0  ndarray(float)
            log(g) to obtain

        Z0 ndarray(float)
            metallicity values

        L0 ndarray(float)
            luminosity values

        dT_max: float
            If, T2 (resp. T1) is too far from T compared to T1 (resp. T2),
            i2 (resp. i1) is not used.
            (see below for namings)

        eps: float
            temperature sensitivity under which points are considered to
            have the same temperature

        weights: ndarray(float)
            luminosity weigths to apply after interpolation

        pool: Pool-like object
            specify a multiprocessing pool for parallel processing

        nthreads: int
            number of processes to use by default

        returns
        -------
        r: ndarray
            Returns 3 to 12 star indexes and associated weights

        see __interp__

        TODO: compute new weights accounting for Z
        """
        r = __interpMany__(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-06, weights=weights, pool=pool, nthreads=nthreads)
        idx = np.unique(r[:, 0])
        d = {}
        for idxk in idx:
            d[idxk] = 0.
        for k in range(len(r)):
            d[ r[k, 0] ] += r[k, 1]
        del r, idx
        return np.asarray(list(d.items()))

    def points_inside(self, xypoints, dlogT=0.1, dlogg=0.3):
        """
        Returns if a point is inside the polygon defined by the boundary of the library

        Parameters
        ----------
        xypoints: sequence
            a sequence of N logg, logT pairs.

        dlogT: float
            margin in logT

        dlogg: float
            margin in logg

        returns
        -------
        r: ndarray(dtype=bool)
            a boolean ndarray, True for points inside the polygon.
            A point on the boundary may be treated as inside or outside.
        """
        p = self.get_boundaries(dlogT=dlogT, dlogg=dlogg)
        return p.contains_points(xypoints)

    def get_radius(self, logl, logt):
        """ Returns the radius of a star given its luminosity and temperature

        Assuming a black body, it comes:
                R ^ 2 = L / ( 4 \pi \sigma T ^ 4 ),

        with:
            L, luminosity in W,
            pi, 3.141592...
            sig, Stephan constant in  W * m**-2 * K**-4
            T, temperature in K

        Parameters
        ----------
        logl: ndarray[float, ndim=1]
            log luminosities from the isochrones, in Lsun

        logt: ndarray[float, ndim=1]
            log temperatures from the isochrones, in K

        returns
        -------
        radii: ndarray[float, ndim=1]
            array of radii in m (SI units)
        """
        return np.sqrt( np.power(10.,logl) * lsun
                        / (4.0 * np.pi * sig_stefan
                           * np.power(np.power(10.,logt), 4.) ) )

    def get_boundaries(self, dlogT=0.1, dlogg=0.3, **kwargs):
        """ Returns the closed boundary polygon around the stellar library with
        given margins

        Parameters
        ----------
        s: Stellib
            Stellar library object

        dlogT: float
            margin in logT

        dlogg: float
            margin in logg

        returns
        -------
        b: ndarray[float, ndim=2]
            closed boundary edge points: [logT, logg]

        .. note::

            as computing the boundary could take time, it is saved in the object
            and only recomputed when parameters are updated
        """
        # if bbox is defined then assumes it is more precise and use it instead.
        if hasattr(self, 'bbox'):
            return Path(self.bbox(dlogT, dlogg))

        if getattr(self, '_bound', None) is not None:
            # check if recomputing is needed
            if ((self._bound[1] - dlogT) < 1e-3) and (abs(self._bound[2] - dlogg) < 1e-3):
                return self._bound[0]

        leftb   = [(np.max(self.logT[self.logg == k]) + dlogT, k ) for k in np.unique(self.logg)]
        leftb  += [(leftb[-1][1], leftb[-1][0] + dlogg)]
        leftb   = [(leftb[0][1], leftb[0][0] - dlogg)] + leftb

        rightb  = [(np.min(self.logT[self.logg == k]) - dlogT, k) for k in np.unique(self.logg)[::-1]]
        rightb += [(rightb[-1][1], rightb[-1][0] - dlogg)]
        rightb  = [(rightb[0][1], rightb[0][0] + dlogg)] + rightb

        b = leftb + rightb
        b += [b[0]]

        self._bound = (Path(np.array(b)), dlogT, dlogg)
        return self._bound[0]

    def genQ(self, qname, r, **kwargs):
        """ Generate a composite value from a previously calculated
            interpolation
            Works on 1 desired star or a population of stars

        Parameters
        ----------
        qname: str
            quantity name from self.grid

        r: ndarray
            the result from a previous interpolation

        Returns
        -------
        val: ndarray
            an array containing the value
        """
        return ( self.grid[qname][r[:, 0].astype(int)] * r[:, 1] ).sum()

    def genSpectrum(self, T0, g0=None, Z0=None, weights=None, **kwargs):
        """ Generate a composite sprectrum
        Does the interpolation or uses a previously calculated
        interpolation
        Works on 1 desired star or a population of stars

        if T0 and g0 are iterable, it calls interpMany

        Parameters
        ----------
        T0: float or sequence
            log(Teff) of each star or a 2d-array containing the result from a
            previous interpolation

        g0: float or sequence
            log(g) of each stars

        Z0: float or sequence
            metallicity

        weights: float or sequence
            individual weights of each star

        **kwargs: forwarded to interpMany

        Returns
        -------
        s: ndarray
            an array containing the composite spectrum
        """
        if Z0 is not None:
            if hasattr(T0, '__iter__'):
                _r = self.interpMany(T0, g0, Z0, weights=weights, **kwargs)
            else:
                _r = np.asarray(self.interp(T0, g0, Z0, weight=weights, **kwargs))
        else:
            assert( T0.ndim == 2), 'error expecting 2d-array'
            _r = T0
        return ( ( (self.spectra[_r[:, 0].astype(int)].T) * _r[:, 1]) ).sum(1)

    def gen_spectral_grid_from_given_points(self, pts, bounds=dict(dlogT=0.1,
                                                                   dlogg=0.3)):
        """
        Reinterpolate a given stellar spectral library on to an Isochrone grid

        Parameters
        ----------
        pts: dict like structure of points
            dictionary like or named data structure of points to interpolate at.
            pts must contain:
            logg  surface gravity in log-scale
            logT  log of effective temperatures (in Kelvins)
            logL  log of luminosity in Lsun units
            Z     metallicity

        bounds: dict
            sensitivity to extrapolation (see grid.get_stellib_boundaries)
            default: {dlogT:0.1, dlogg:0.3}

        Returns
        -------
        g: SpectralGrid
            Spectral grid (in memory) containing the requested list of
            stars and associated spectra
        """
        # Step 0: prepare outputs
        # =======================
        # Grid properties will be stored into a dictionary format until saved on disk
        # SEDs are kept into a ndarray
        ndata = len(pts)
        _grid  = {}
        _grid['radius'] = np.empty(ndata, dtype=float )

        # stores the grid+prior weights (initialize to 1)
        _grid['weight'] = np.full(ndata, 1.0, dtype=float)

        # stores the prior weights separately (initialize to 1)
        #   This will allow for adjustable priors and
        #   visualization of the priors themselves as weights include
        #   the grid weights to correct for the non-uniform grid spacing
        _grid['prior_weight'] = np.full(ndata, 1.0, dtype=float)

        # stores the grid weights separately (initialize to 1)
        #   these weights alone provide flat priors on all fit parameters
        _grid['grid_weight'] = np.full(ndata, 1.0, dtype=float)

        # index to the grid
        # useful to setup here as it will then be cleanly propagated
        #   to the SED grid
        _grid['specgrid_indx'] = np.full(ndata, 0.0, dtype=float)

        specs = np.empty( (ndata, len(self.wavelength)), dtype=float )

        # copy meta data of pts into the resulting structure
        if hasattr(pts, 'keys'):
            for key in list(pts.keys()):
                _grid[key] = np.asarray(pts[key])
        elif hasattr(pts, 'dtype'):
            if pts.dtype.names is not None:
                for key in pts.dtype.names:
                    _grid[key] = pts[key]
            else:
                raise AttributeError('pts is expected to have named items (keys or dtype.names)')

        # Step 1: Avoid Extrapolation
        # ===========================
        # check boundary conditions, keep the data but do not compute the sed if not needed
        bound_cond = self.points_inside(list(zip(pts['logT'], pts['logg'])))

        # Step 2: radii
        # =============
        # Stellar library models are given in cm^-2  ( 4 pi R)
        # Compute radii of each point using log(T) and log(L)
        radii = self.get_radius(pts['logL'], pts['logT'])
        _grid['radius'] = radii[:] / rsun

        # weights to apply during the interpolation
        # note that radii must be in cm
        weights = 4. * np.pi * (radii * 1e2) ** 2

        # Step 3: Interpolation
        # =====================
        # Do the actual interpolation, avoiding exptrapolations
        for mk, (rT, rg, rZ) in Pbar(ndata, desc='spectral grid').iterover(enumerate(zip(pts['logT'], pts['logg'], pts['Z']))):
            if bound_cond[mk]:
                s = np.array( self.interp(rT, rg, rZ, 0.) ).T
                specs[mk, :] = self.genSpectrum(s) * weights[mk]

        # Step 4: filter points without spectrum
        # ======================================
        idx = np.array(bound_cond)

        lamb = self.wavelength[:]
        specs = specs.compress(idx, axis=0)
        for k in list(_grid.keys()):
                _grid[k] = _grid[k].compress(idx, axis=0)

        # Step 5: Ship
        # ============
        header = {'stellib': self.source,
                  'comment': 'radius in Rsun',
                  'name': 'Reinterpolated stellib grid'}

        # populate the specgrid index
        _grid['specgrid_indx'] = np.arange(len(_grid['specgrid_indx']),
                                           dtype=np.int64)

        # ship
        g = SpectralGrid(lamb, seds=specs, grid=Table(_grid), header=header, backend='memory')

        return g

    def plot_boundary(self, ax=None, dlogT=0., dlogg=0., **kwargs):
        """
        Parameters
        ----------

        dlogT: float
            margin in logT (see get_boundaries)

        dlogg: float
            margin in logg (see get_boundaries)

        agg_filter: unknown
        alpha: float or None
        animated: [True | False]
        antialiased or aa: [True | False]  or None for default
        axes: an :class:`~matplotlib.axes.Axes` instance
        clip_box: a :class:`matplotlib.transforms.Bbox` instance
        clip_on: [True | False]
        clip_path: [ (:class:`~matplotlib.path.Path`, :class:`~matplotlib.transforms.Transform`) |  :class:`~matplotlib.patches.Patch` | None ]
        color: matplotlib color spec
        contains: a callable function
        edgecolor or ec: mpl color spec, or None for default, or 'none' for no color
        facecolor or fc: mpl color spec, or None for default, or 'none' for no color
        figure: a :class:`matplotlib.figure.Figure` instance
        fill: [True | False]
        gid: an id string
        hatch: [ '/' | '\\' | '|' | '-' | '+' | 'x' | 'o' | 'O' | '.' | '*' ]
        label: string or anything printable with '%s' conversion.
        linestyle or ls: ['solid' | 'dashed' | 'dashdot' | 'dotted']
        linewidth or lw: float or None for default
        lod: [True | False]
        path_effects: unknown
        picker: [None|float|boolean|callable]
        rasterized: [True | False | None]
        snap: unknown
        transform: :class:`~matplotlib.transforms.Transform` instance
        url: a url string
        visible: [True | False]
        zorder: any number

        .. seealso::

            :class:`Patch`
                For additional kwargs
        """
        import matplotlib.patches as patches
        from pylab import gca
        if ax is None:
            ax = gca()
        p = self.get_boundaries(dlogT=dlogT, dlogg=dlogg)
        ax.add_patch(patches.PathPatch(p, **kwargs))
        return p

    def __add__(self, other):
        if not isNestedInstance(other, Stellib):
            raise ValueError('expecting a Stellib object, got {0}'.format(type(other)))

        return CompositeStellib([self, other])

    def __repr__(self):
        return "{0:s}, ({1:s})\n{2:s}".format(self.name,
                                              nbytes(self, pprint=True),
                                              object.__repr__(self))


class CompositeStellib(Stellib):
    """ Generates an object from the union of multiple individual libraries """
    def __init__(self, osllist, *args, **kwargs):
        self._olist = osllist

    def __add__(self, other):
        if not isNestedInstance(other, Stellib):
            raise ValueError('expecting a Stellib object, got {0}'.format(type(other)))

        lst = [k for k in self._olist] + [other]
        return CompositeStellib(lst)

    def __radd__(self, other):
        if not isNestedInstance(other, Stellib):
            raise ValueError('expecting a Stellib object, got {0}'.format(type(other)))

        lst = [other] + [k for k in self._olist]
        return CompositeStellib(lst)

    @property
    def wavelength(self):
        """ return a common wavelength sampling to all libraries. This can be
        used to reinterpolate any spectrum onto a common definition """

        lambs = np.unique(np.asarray([ osl.wavelength[:]
                                       for osl in self._olist ]))
        return lambs

    @property
    def source(self):
        return ' + '.join([k.name for k in self._olist])

    def which_osl(self, xypoints, dlogT=0., dlogg=0.):
        """
        Returns the library indice that contains each point in xypoints

        The decision is made from a two step search:
            * first, each point is checked against the strict boundary of each
              library (i.e., dlogT = 0, dlogg = 0).
            * second, if points are not found in strict mode, the boundary is
              relaxed and a new search is made.

        Each point is associated to the first library matching the above
            conditions.

        Parameters
        ----------
        xypoints: sequence
            a sequence of N logg, logT pairs.

        dlogT: float
            margin in logT

        dlogg: float
            margin in logg

        returns
        -------
        res: ndarray(dtype=int)
            a ndarray, 0 meaning no library covers the point, and 1, ... n,
               for the n-th library
        """
        xy = np.asarray(xypoints)

        # check that all points are in the full boundary area
        # MF: testing why points_inside does not agree on all computers...
        # as we do not keep individual results, no need to store then all
        # first, collapse directly

        #res_temp = np.zeros((len(xy),len(self._olist)))
        #for ek,ok in enumerate(self._olist):
        #    res_temp[:, ek] = ok.points_inside(xy, dlogT=dlogT,
        #                                       dlogg=dlogg).astype(int)
        res_temp = np.zeros(len(xy), dtype=int)
        for ek, ok in enumerate(self._olist):
            res_temp += ok.points_inside(xy, dlogT=dlogT,
                                         dlogg=dlogg).astype(int)

        ind = res_temp > 0
        res = np.zeros(len(xy), dtype=int)
        res[ind] = 1
        res = res - 1

        #res = self.points_inside(xy, dlogT=dlogT, dlogg=dlogg).astype(int) - 1
        # if res == -1: invalid point, res == 0: proceed

        if max(res) < 0:
            # DEBUG: should generate an exeception in further functions
            # TODO: get rid and replace
            return
            # return res

        # Strict mode
        # ===========
        # Not extrapolation allowed >> dlogT = 0, dlogg = 0
        # 0 is used to flag points without a matching library yet
        # libraries are then indexed from 1 to n
        # -1 means point outside the compound library
        for ek, ok in enumerate(self._olist):
            if 0 in res:
                ind = np.atleast_1d(np.squeeze(np.where(res == 0)))
                r = ok.points_inside(xy[ind], dlogT=0., dlogg=0.)
                res[ind[r]] = ek + 1

        # Relaxed mode
        # ============
        # In this case we accept some flexibility in the boundary limits,
        # which allows limited extrapolation ranges.
        # this only affects points not already matched
        if 0 in res:
            for ek, ok in enumerate(self._olist):
                if 0 in res:
                    ind = np.atleast_1d(np.squeeze(np.where(res == 0)))
                    r = ok.points_inside(xy[ind], dlogT=dlogT, dlogg=dlogg)
                    res[ind[r]] = ek + 1
        return res

    def __repr__(self):
        return "CompositeStellib, {0}\n{1}".format(object.__repr__(self),
                                                [k.name for k in self._olist])

    def get_boundaries(self, dlogT=0.1, dlogg=0.3, **kwargs):
        """ Returns the closed boundary polygon around the stellar library with
        given margins

        Parameters
        ----------
        s: Stellib
            Stellar library object

        dlogT: float
            margin in logT

        dlogg: float
            margin in logg

        returns
        -------
        b: ndarray[float, ndim=2]
            (closed) boundary points: [logg, Teff] (or [Teff, logg] is swap
            is True)

        .. note::
            as computing the boundary could take time, it is saved in the object
            and only recomputed when parameters are updated
        """
        if getattr(self, '_bound', None) is not None:
            if (((self._bound[1] - dlogT) < 1e-3) and
                (abs(self._bound[2] - dlogg) < 1e-3)):
                return self._bound[0]

        b = [osl.get_boundaries(dlogT=dlogT, dlogg=dlogg, **kwargs)
             for osl in self._olist]
        self._bound = (Path.make_compound_path(*b), dlogT, dlogg)
        return self._bound[0]

    def interp(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6, bounds={}):
        """ Interpolation of the T,g grid

        Interpolate on the grid and returns star indices and
        associated weights, and Z.
        3 to 12 stars are returned.
        It calls _interp_, but reduce the output to the relevant stars.

        Parameters
        ----------
        T0: double
            log(Teff) to obtain

        g0: double
            log(g) to obtain

        T:  double
            log(Teff) of the grid

        g:  double
            log(g) of the grid

        dT_max: float
            If, T2 (resp. T1) is too far from T compared to T1 (resp. T2),
            i2 (resp. i1) is not used.
            (see below for namings)

        eps: foat
            temperature sensitivity under which points are considered to
            have the same temperature

        bounds: dict
            sensitivity to extrapolation (see `:func: Stellib.get_boundaries`)
            default: {dlogT:0.1, dlogg:0.3}

        returns
        -------
        (osl, r): tuple
            osl: is the library index starting from 1. 0 means no coverage.
            r: is the result from interp call on the corresponding library.
            a 3 to 12 star indexes and associated weights
        """
        dlogT = bounds.get('dlogT', 0.1)
        dlogg = bounds.get('dlogg', 0.3)

        osl_index = self.which_osl(np.atleast_2d([T0, g0]),
                                   dlogT=dlogT, dlogg=dlogg)[0]

        if osl_index > 0:
            return (osl_index,
                    self._olist[osl_index - 1].interp(self, T0, g0, Z0, L0,
                                                      dT_max=0.1, eps=1e-6))
        else:
            return [(0, None)]

    def interpMany(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6, weights=None,
                   bounds={}, pool=None, nthreads=__NTHREADS__):
        """ run interp on a list of inputs and returns reduced results

        Interpolation of the T,g grid at Z0 metallicity

        Interpolate on the grid and returns star indices and
        associated weights, and Z.
        3 to 12 stars are returned.
        It calls _interp_, but reduce the output to the relevant stars.

        Parameters
        ----------
        T0: ndarray(float)
            log(Teff) to obtain

        g0: ndarray(float)
            log(g) to obtain

        Z0: ndarray(float)
            metallicity values

        L0: ndarray(float)
            luminosity values

        dT_max: float
            If, T2 (resp. T1) is too far from T compared to T1 (resp. T2),
            i2 (resp. i1) is not used.
            (see below for namings)

        eps: foat
            temperature sensitivity under which points are considered to
            have the same temperature

        weights: ndarray(float)
            luminosity weigths to apply after interpolation

        bounds: dict
            sensitivity to extrapolation (see `:func: Stellib.get_boundaries`)
            default: {dlogT:0.1, dlogg:0.3}

        pool: Pool-like object
            specify a multiprocessing pool for parallel processing

        nthreads: int
            number of processes to use by default

        returns
        -------
        (osl, r): tuple
            osl is the library index starting from 1. 0 means no coverage.
            r is the result from interp call on the corresponding library.
            a 3 to 12 star indexes and associated weights
        """
        dlogT = bounds.get('dlogT', 0.1)
        dlogg = bounds.get('dlogg', 0.3)

        osl_index = self.which_osl(list(zip(T0, g0)), dlogT=dlogT, dlogg=dlogg)

        g = []
        for oslk, osl in enumerate(self._olist):
            # make a generator to avoid keeping all in memory
            ind = np.where(osl_index - 1 == oslk)
            if np.squeeze(ind).size is not 0:
                g.append( [oslk + 1, osl.interpMany(T0[ind], g0[ind],
                                                    Z0[ind], L0[ind],
                                                    dT_max=dT_max, eps=eps,
                                                    weights=weights,
                                                    pool=pool,
                                                    nthreads=nthreads)] )

        return g

    def genQ(self, qname, r, **kwargs):
        """ Generate a composite value from a previously calculated
            interpolation
            Works on 1 desired star or a population of stars

        Parameters
        ----------
        qname: str
            quantity name from self.grid

        r: (osl, r) tuple
            osl: is the library index starting from 1. 0 means no coverage.
            r: is the result from interp call on the corresponding library.

        returns
        -------
        q: float
            value (from weighted sum)
        """
        vals = 0.
        for _osl, _r in r:
            if _osl > 0:
                vals += self._olist[_osl - 1].genQ(qname, _r, **kwargs)
        return vals

    def genSpectrum(self, T0, g0=None, Z0=None, weights=None, **kwargs):
        """ Generate a composite sprectrum
            Does the interpolation or uses a previously calculated interpolation
            Works on 1 desired star or a population of stars

        Parameters
        ----------
        T0: ndarray(float)
            log(Teff) to obtain

        g0: ndarray(float)
            log(g) to obtain

        Z0: ndarray(float)
            metallicity values

        weights: ndarray(float)
            individual weights of each star

        **kwargs forwarded to interp(Many)


        returns
        -------
        s: ndarray
            an array containing the composite spectrum reinterpolated onto
            self.wavelength

        .. note::

            if T0 and g0 are iterable, it calls interpMany
        """
        if Z0 is not None:
            if hasattr(T0, '__iter__'):
                _r = self.interpMany(T0, g0, Z0, weights=weights, **kwargs)
            else:
                _r = np.asarray(self.interp(T0, g0, Z0, weight=weights,
                                            **kwargs))
        else:
            _r = T0

        l0 = self.wavelength
        s = np.zeros(len(l0), dtype=float)
        for osl, _rk in _r:
            if osl > 0:
                sp = (((self._olist[osl - 1].spectra[_r[:, 0].astype(int)].T)
                       * _r[:, 1]) ).sum(1)
                lamb = self._olist[osl - 1].wavelength
                s += np.interp(l0, lamb, sp)
        return s

    def gen_spectral_grid_from_given_points(self, pts,
                                            bounds=dict(dlogT=0.1, dlogg=0.3)):
        """
        Reinterpolate a given stellar spectral library on to an Isochrone grid

        Parameters
        ----------
        pts: dict like structure of points
            dictionary like or named data structure of points to interpolate at.
            pts must contain:
            logg  surface gravity in log-scale
            logT  log of effective temperatures (in Kelvins)
            logL  log of luminosity in Lsun units
            Z     metallicity

        bounds: dict
            sensitivity to extrapolation (see `:func: Stellib.get_boundaries`)
            default: {dlogT:0.1, dlogg:0.3}

        Returns
        -------
        g: SpectralGrid
            Spectral grid (in memory) containing the requested list of stars
            and associated spectra
        """
        dlogT = bounds.get('dlogT', 0.1)
        dlogg = bounds.get('dlogg', 0.3)

        osl_index = self.which_osl(list(zip(pts['logT'], pts['logg'])),
                                   dlogT=dlogT, dlogg=dlogg)

        seds = []
        grid = []
        l0 = self.wavelength
        for oslk, osl in enumerate(self._olist):
            # make a generator to avoid keeping all in memory
            _pts = {}
            #oslk + 1 since 0 corresponds to "not covered by any osl"
            ind = (osl_index == (oslk + 1))
            #print sum(ind)
            if np.sum(ind) > 0:
                if hasattr(pts, 'keys'):
                    keys = list(pts.keys())
                elif hasattr(pts, 'dtype'):
                    keys = pts.dtypes.names
                else:
                    raise AttributeError('Input pts is expected to have named fields')
                for k in keys:
                    _pts[k] = pts[k][ind]
                #keep track of the spectra library that is selected
                _pts['osl'] = osl_index[ind]
                _pts = Table(_pts)
            #_pts = [ (logg, logT, logL, Z) for (logg, logT, logL, Z, ok) in zip(pts['logg'], pts['logT'], pts['logL'], pts['Z'], osl_index) if (ok - 1 == oslk) ]
                gk = osl.gen_spectral_grid_from_given_points(_pts, bounds=dict(dlogT=0.1, dlogg=0.3))
                gk.seds = interp1d(gk.lamb, gk.seds, axis=1)(l0)
                seds.append(gk.seds)
                grid.append(gk.grid)

        header = {'stellib': self.source,
                  'comment': 'radius in Rsun',
                  'name': 'Reinterpolated stellib grid'}

        _grid = recfunctions.stack_arrays( grid, defaults=None, usemask=False, asrecarray=True)

        # populate the specgrid index
        _grid['specgrid_indx'] = np.arange(len(_grid['specgrid_indx']),
                                           dtype=np.int64)

        g = SpectralGrid(l0, seds=np.vstack(seds), grid=Table(_grid), header=header, backend='memory')

        return g


class Elodie(Stellib):
    """ Elodie 3.1 stellar library derived class

    This library matches BaSeL 2.2 grid definition

    References
    ----------
    Prugniel et al 2007, astro-ph/703658
    """
    def __init__(self, *args, **kwargs):
        self.name = 'ELODIE v3.1'
        self.source = config['elodie_3.1']
        self._load_()

    def _load_(self):
        g = SpectralGrid(self.source, backend='memory')
        self.wavelength = g.lamb
        self.grid = g.grid
        self.grid.header['NAME'] = self.name
        self.spectra = g.seds

    def bbox(self, dlogT=0.05, dlogg=0.25):
        """ Boundary of Elodie library

        Parameters
        ----------
        dlogT: float
            log-temperature tolerance before extrapolation limit

        dlogg: float
            log-g tolerance before extrapolation limit

        Returns
        -------
        bbox: ndarray
            (logT, logg) edges of the bounding polygon
        """
        bbox = [(3.301 - dlogT, 5.500 + dlogg),
                (3.301 - dlogT, 3.500 - dlogg),
                (3.544 - dlogT, 3.500 - dlogg),
                (3.544 - dlogT, 1.000),
                (3.477, 0.600 + dlogg),
                (3.447 - dlogT, 0.600 + dlogg),
                (3.398 - dlogT, 0.280 + dlogg),
                (3.398 - dlogT, -1.020 - dlogg),
                (3.398, -1.020 - dlogg),
                (3.447, -1.020 - dlogg),
                (3.505 + dlogT, -0.700 - dlogg),
                (3.544 + dlogT, -0.510 - dlogg),
                (3.574 + dlogT, -0.290 - dlogg),
                (3.602 + dlogT, 0.000 - dlogg),
                (3.778, 0.000 - dlogg),
                (3.778 + dlogT, 0.000),
                (3.875 + dlogT, 0.500),
                (3.929 + dlogT, 1.000),
                (3.954 + dlogT, 1.500),
                (4.021 + dlogT, 2.000 - dlogg),
                (4.146, 2.000 - dlogg),
                (4.146 + dlogT, 2.000),
                (4.279 + dlogT, 2.500),
                (4.415 + dlogT, 3.000),
                (4.491 + dlogT, 3.500),
                (4.544 + dlogT, 4.000),
                (4.602 + dlogT, 4.500),
                (4.699 + dlogT, 5.000 - dlogg),
                (4.699 + dlogT, 5.000 + dlogg),
                (3.525 + dlogT, 5.000 + dlogg),
                (3.525 + dlogT, 5.500 + dlogg),
                (3.301 - dlogT, 5.500 + dlogg) ]

        return np.array(bbox)

    @property
    def logg(self):
        return self.grid['logg']

    @property
    def logT(self):
        return self.grid['logT']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']

    @property
    def NHI(self):
        return self.grid['NHI']

    @property
    def NHeI(self):
        return self.grid['NHeI']

    @property
    def NHeII(self):
        return self.grid['NHeII']


class BaSeL(Stellib):
    """ BaSeL 2.2 (This library is used in Pegase.2)
        This library is used in Pegase.2

    The BaSeL stellar spectral energy distribution (SED) libraries are libraries of
    theoretical stellar SEDs recalibrated using empirical photometric data.
    Therefore, we call them semi-empirical libraries.

    The BaSeL 2.2 library was calibrated using photometric data from solar
    metallicity stars.

    References
    ----------
    * Lejeune, Cuisiner, and Buser, 1998 A&AS, 130, 65
    * can be downloaded http://www.astro.unibas.ch/BaSeL_files/BaSeL2_2.tar.gz
    """
    def __init__(self, *args, **kwargs):
        self.name = 'BaSeL 2.2 (Pegase.2 version)'
        self.source = config['basel_2.2_pegase']
        self._load_()

    def _load_(self):
        g = SpectralGrid(self.source, backend='memory')
        self.wavelength = g.lamb
        self.grid = g.grid
        self.grid.header['NAME'] = 'Basel 2.2 (pegase)'
        self.spectra = g.seds

    def bbox(self, dlogT=0.05, dlogg=0.25):
        """ Boundary of Basel 2.2 library

        Parameters
        ----------
        dlogT: float
            log-temperature tolerance before extrapolation limit

        dlogg: float
            log-g tolerance before extrapolation limit

        Returns
        -------
        bbox: ndarray
            (logT, logg) edges of the bounding polygon
        """
        bbox = [(3.301 - dlogT, 5.500 + dlogg),
                (3.301 - dlogT, 3.500 - dlogg),
                (3.544 - dlogT, 3.500 - dlogg),
                (3.544 - dlogT, 1.000),
                (3.477, 0.600 + dlogg),
                (3.447 - dlogT, 0.600 + dlogg),
                (3.398 - dlogT, 0.280 + dlogg),
                (3.398 - dlogT, -1.020 - dlogg),
                (3.398, -1.020 - dlogg),
                (3.447, -1.020 - dlogg),
                (3.505 + dlogT, -0.700 - dlogg),
                (3.544 + dlogT, -0.510 - dlogg),
                (3.574 + dlogT, -0.290 - dlogg),
                (3.602 + dlogT, 0.000 - dlogg),
                (3.778, 0.000 - dlogg),
                (3.778 + dlogT, 0.000),
                (3.875 + dlogT, 0.500),
                (3.929 + dlogT, 1.000),
                (3.954 + dlogT, 1.500),
                (4.021 + dlogT, 2.000 - dlogg),
                (4.146, 2.000 - dlogg),
                (4.146 + dlogT, 2.000),
                (4.279 + dlogT, 2.500),
                (4.415 + dlogT, 3.000),
                (4.491 + dlogT, 3.500),
                (4.544 + dlogT, 4.000),
                (4.602 + dlogT, 4.500),
                (4.699 + dlogT, 5.000 - dlogg),
                (4.699 + dlogT, 5.000 + dlogg),
                (3.525 + dlogT, 5.000 + dlogg),
                (3.525 + dlogT, 5.500 + dlogg),
                (3.301 - dlogT, 5.500 + dlogg) ]

        return np.array(bbox)

    @property
    def logg(self):
        return self.grid['logg']

    @property
    def logT(self):
        return self.grid['logT']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']

    @property
    def logZ(self):
        return self.grid['logZ']

    @property
    def NHI(self):
        return self.grid['NHI']

    @property
    def NHeI(self):
        return self.grid['NHeI']

    @property
    def NHeII(self):
        return self.grid['NHeII']


class Kurucz(Stellib):
    """
    The stellar atmosphere models by Castelli and Kurucz 2004 or ATLAS9

    * LTE
    * PP
    * line blanketing
    """
    def __init__(self, filename=None, *args, **kwargs):
        self.name = 'Kurucz 2004'
        if filename is None:
            self.source = config['kurucz']
        else:
            self.source = filename
        self._load_()

    def _load_(self,):
        g = SpectralGrid(self.source, backend='memory')
        self.wavelength = g.lamb
        self.grid = g.grid
        self.grid.header['NAME'] = self.name
        self.spectra = g.seds

    def bbox(self, dlogT=0.05, dlogg=0.25):
        """ Boundary of Kurucz 2004 library

        Parameters
        ----------
        dlogT: float
            log-temperature tolerance before extrapolation limit

        dlogg: float
            log-g tolerance before extrapolation limit

        Returns
        -------
        bbox: ndarray
            (logT, logg) edges of the bounding polygon
        """
        bbox = [(3.54406 - dlogT, 5.000 + dlogg),
                (3.55403 - dlogT, 0.000 - dlogg),
                (3.778, 0.000 - dlogg),
                (3.778 + dlogT, 0.000),
                (3.875 + dlogT, 0.500),
                (3.929 + dlogT, 1.000),
                (3.954 + dlogT, 1.500),
                (4.146, 2.000 - dlogg),
                (4.146 + dlogT, 2.000),
                (4.279 + dlogT, 2.500),
                (4.415 + dlogT, 3.000),
                (4.491 + dlogT, 3.500),
                (4.591 + dlogT, 4.000),
                (4.689 + dlogT, 4.500),
                (4.699 + dlogT, 5.000 + dlogg),
                (3.544 - dlogT, 5.000 + dlogg) ]

        return np.array(bbox)

    @property
    def logT(self):
        return self.grid['logT']

    @property
    def logg(self):
        return self.grid['logg']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']

    @property
    def logZ(self):
        return self.grid['logz']


class Tlusty(Stellib):
    """
    Tlusty O and B stellar atmospheres

    * NLTE
    * Parallel Planes
    * line blanketing

    References
    ----------
    Hubeny 1988 for initial reference
    Lanz, T., & Hubeny, I. (2003) for more recent (NL TE) developments

    * **OSTAR2002 Grid**: O-type stars, 27500 K <= Teff <= 55000 K
        * Reference: Lanz & Hubeny (2003)

    * **BSTAR2006 Grid**: Early B-type stars, 15000 K <= Teff <= 30000 K
            * Reference: Lanz & Hubeny (2007)

    files are available at: http://nova.astro.umd.edu/Tlusty2002/database/

    O and B stars rebinned to nearly 20,000 frequency points (for CLOUDY usage)
    http://nova.astro.umd.edu/Tlusty2002/database/obstar_merged_3d.ascii.gz
    """
    def __init__(self, filename=None, *args, **kwargs):
        self.name = 'Tlusty'
        if filename is None:
            self.source = config['tlusty']
        else:
            self.source = filename
        self._load_()

    def _load_(self):
        g = SpectralGrid(self.source, backend='memory')
        self.wavelength = g.lamb
        self.grid = g.grid
        self.grid.header['NAME'] = 'tlusty'
        self.spectra = g.seds

    def bbox(self, dlogT=0.05, dlogg=0.25):
        """ Boundary of Tlusty library

        Parameters
        ----------
        dlogT: float
            log-temperature tolerance before extrapolation limit

        dlogg: float
            log-g tolerance before extrapolation limit

        Returns
        -------
        bbox: ndarray
            (logT, logg) edges of the bounding polygon
        """
        bbox = [(4.176 - dlogT, 4.749 + dlogg),
                (4.176 - dlogT, 1.750 - dlogg),
                (4.176 + dlogT, 1.750 - dlogg),
                (4.255 + dlogT, 2.000 - dlogg),
                (4.447 + dlogT, 2.750 - dlogg),
                (4.478 + dlogT, 3.000 - dlogg),
                (4.544 + dlogT, 3.250 - dlogg),
                (4.740 + dlogT, 4.000 - dlogg),
                (4.740 + dlogT, 4.749 + dlogg),
                (4.176 - dlogT, 4.749 + dlogg) ]

        return np.array(bbox)

    @property
    def logT(self):
        return self.grid['logT']

    @property
    def logg(self):
        return self.grid['logg']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']

    @property
    def logZ(self):
        return self.grid['logZ']


class BTSettl(Stellib):
    """
    BT-Settl Library

    References
    ----------

    Paper: Few refereed publications
      Older Ref = http://adsabs.harvard.edu/abs/2000ApJ...539..366A

    Conference Proceedings:
      http://adsabs.harvard.edu/abs/2016sf2a.conf..223A
      http://adsabs.harvard.edu/abs/2012RSPTA.370.2765A

    Files available at: https://phoenix.ens-lyon.fr/Grids/BT-Settl/

    Current Library: AGSS2009 Abundances (due to grid availability)
    Spectra rebinned to match Kurucz, and custom 2 Ang medium resolution
    """
    def __init__(self, medres=True, *args, **kwargs):
        self.name = 'BTSettl'
        if medres:
            self.source = config['btsettl_medres']
        else:
            self.source = config['btsettl']
        self._load_()

    def _load_(self):
        g = SpectralGrid(self.source, backend='memory')
        self.wavelength = g.lamb
        self.grid = g.grid
        self.grid.header['NAME'] = self.name
        self.spectra = g.seds

    def bbox(self, dlogT=0.05, dlogg=0.25):
        """ Boundary of BT-Settl library

        Parameters
        ----------
        dlogT: float
            log-temperature tolerance before extrapolation limit

        dlogg: float
            log-g tolerance before extrapolation limit

        Returns
        -------
        bbox: ndarray
            (logT, logg) edges of the bounding polygon
        """
        bbox = [(3.41497 - dlogT, 6.0 + dlogg),
                (3.41497 - dlogT,-0.5 - dlogg),
                (3.84510 + dlogT,-0.5 - dlogg),
                (4.07918 + dlogT, 0.0 - dlogg),
                (4.17609 + dlogT, 0.5 - dlogg),
                (4.30103 + dlogT, 1.0 - dlogg),
                (4.39794 + dlogT, 1.5 - dlogg),
                (4.47712 + dlogT, 2.0 - dlogg),
                (4.60206 + dlogT, 2.5 - dlogg),
                (4.60206 + dlogT, 3.0 - dlogg),
                (4.69897 + dlogT, 3.5 - dlogg),
                (4.84510 + dlogT, 4.0 - dlogg),
                (4.84510 + dlogT, 4.5 + dlogg),
                (4.00000 + dlogT, 4.5 + dlogg),
                (4.00000 + dlogT, 5.0 + dlogg),
                (3.69897 + dlogT, 5.0 + dlogg),
                (3.69897 + dlogT, 5.5 + dlogg),
                (3.60206 + dlogT, 5.5 + dlogg),
                (3.60206 + dlogT, 6.0 + dlogg),
                (3.41497 - dlogT, 6.0 + dlogg)]

        return np.array(bbox)

    @property
    def logT(self):
        return self.grid['logT']

    @property
    def logg(self):
        return self.grid['logg']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']

    @property
    def logZ(self):
        return self.grid['logZ']


class Munari(Stellib):
    """
    ATLAS9 stellar atmospheres providing higher res than Kurucz
    medium resolution (1 Ang/pix) in optical (2500-10500 Ang)

    References
    ----------

    Paper: Munari et al. 2005 A&A 442 1127
    http://adsabs.harvard.edu/abs/2005A%26A...442.1127M

    Files available at: http://archives.pd.astro.it/2500-10500/
    """
    def __init__(self, *args, **kwargs):
        self.name = 'Munari'
        self.source = config['munari']
        self._load_()

    def _load_(self):
        g = SpectralGrid(self.source, backend='memory')
        self.wavelength = g.lamb
        self.grid = g.grid
        self.grid.header['NAME'] = self.name
        self.spectra = g.seds

    def bbox(self, dlogT=0.05, dlogg=0.25):
        """ Boundary of Munari library

        Parameters
        ----------
        dlogT: float
            log-temperature tolerance before extrapolation limit

        dlogg: float
            log-g tolerance before extrapolation limit

        Returns
        -------
        bbox: ndarray
            (logT, logg) edges of the bounding polygon
        """
        bbox = [(3.54407 - dlogT, 5.0 + dlogg),
                (3.54407 - dlogT, 0.0 - dlogg),
                (3.77815 + dlogT, 0.0 - dlogg),
                (3.87506 + dlogT, 0.5 - dlogg),
                (3.91645 + dlogT, 1.0 - dlogg),
                (3.95424 + dlogT, 1.5 - dlogg),
                (3.98900 + dlogT, 2.0 - dlogg),
                (3.98900 + dlogT, 5.0 + dlogg),
                (3.54407 - dlogT, 5.0 + dlogg)]

        return np.array(bbox)

    @property
    def logT(self):
        return self.grid['logT']

    @property
    def logg(self):
        return self.grid['logg']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']

    @property
    def logZ(self):
        return self.grid['logZ']
