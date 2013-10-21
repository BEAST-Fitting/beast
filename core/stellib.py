"""
Stellib class

Intent to implement a generic module to manage stellar library from various
sources.

The interpolation is implemented from the pegase.2 fortran converted algorithm.
(this may not be pythonic though)
"""
import numpy as np
import pyfits
from .grid import SpectralGrid
from ..external.eztables import Table
from ..external import ezunits
from ..config import __ROOT__, __NTHREADS__
from ..include import __interp__
from ..tools import progressbar
from matplotlib.path import Path
from scipy.interpolate import interp1d
from numpy.lib import recfunctions


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
    return np.asarray(interp(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9])).T


def __interpMany__(oSL, logT, logg, Z, logL, dT_max=0.1, eps=1e-06, weights=None, pool=None, nthreads=__NTHREADS__):
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
        seq = [ (logT[k], logg[k], Z, logL[k], oSL.Teff, oSL.logg, oSL.Z, dT_max, eps, 1.) for k in range(len(logT)) ]
    else:
        seq = [ (logT[k], logg[k], Z, logL[k], oSL.Teff, oSL.logg, oSL.Z, dT_max, eps, weights[k]) for k in range(len(logT)) ]

    if (pool is not None):
        r = pool.map( __interpSingle__, seq )
    else:
        r = map( __interpSingle__, seq )

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
        pass

    def interp(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6):
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
        _Z    = self.Z
        _T    = np.asarray(self.grid['Teff'], dtype=np.double)
        _g    = np.asarray(self.grid['logG'], dtype=np.double)
        return interp(T0, g0, Z0, L0, _T, _g, _Z, dT_max=0.1, eps=1e-6)

    def interpMany(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6, weights=None, pool=None, nthreads=__NTHREADS__):
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

        see __interp__

        TODO: compute new weights accounting for Z
        """
        #_t = np.asarray(T0)
        #_g = np.asarray(g0)
        r = __interpMany__(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-06, weights=weights, pool=pool, nthreads=nthreads)
        idx = np.unique(r[:, 0])
        # d = { idxk:0. for idxk in idx }
        d = {}
        for idxk in idx:
            d[idxk] = 0.
        for k in xrange(len(r)):
            d[ r[k, 0] ] += r[k, 1]
        del r, idx
        return np.asarray(d.items())

    def points_inside(self, xypoints, dlogT=0.1, dlogg=0.3):
        """
        Returns if a point is inside the polygon defined by the boundary of the library

        keywords
        --------

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
        p = self.get_boundaries(dlogT=dlogT, dlogg=dlogg, closed=True)
        return p.contains_points(xypoints)

    def get_radius(self, logl, logt):
        """ Returns the radius of a star given its luminosity and temperature

        Assuming a black body, it comes:
                R ^ 2 = L / ( 4 pi sig T ^ 4 ),

        with:
            L, luminosity in W,
            pi, 3.141592...
            sig, Stephan constant in  W * m**-2 * K**-4
            T, temperature in K

        keywords
        --------

        logl: ndarray[float, ndim=1]
            log luminosities from the isochrones, in Lsun

        logt: ndarray[float, ndim=1]
            log temperatures from the isochrones, in K

        returns
        -------
        radii: ndarray[float, ndim=1]
            array of radii in m (SI units)
        """
        lsun = 1. * ezunits.unit['lsun'].to('W').magnitude  # 3.839e26 W
        sig  = 5.67037321 * 1e-8 * ezunits.unit[' W * m**-2 * K**-4'].magnitude
        return np.sqrt( (10 ** logl) * lsun / (4.0 * np.pi * sig * ((10 ** logt) ** 4)) )

    def get_boundaries(self, dlogT=0.1, dlogg=0.3, closed=True):
        """ Returns the closed boundary polygon around the stellar library with
        given margins

        keywords
        --------
        s: Stellib
            Stellar library object

        dlogT: float
            margin in logT

        dlogg: float
            margin in logg

        closed: bool
            if set, close the polygon

        returns
        -------
        b: ndarray[float, ndim=2]
            (closed) boundary points: [logg, Teff]

        Note
        ----
        as computing the boundary could take time, it is saved in the object
        and only recomputed when parameters are updated
        """
        if getattr(self, '_bound', None) is not None:
            if ((self._bound[1] - dlogT) < 1e-3) and (abs(self._bound[2] - dlogg) < 1e-3):
                return self._bound[0]

        leftb   = [(k, np.max(self.logT[self.logg == k]) + dlogT ) for k in np.unique(self.logg)]
        leftb  += [ (leftb[-1][0] + dlogg, leftb[-1][1]) ]
        leftb   = [ (leftb[0][0] - dlogg, leftb[0][1]) ] + leftb
        rightb  = [(k, np.min(self.logT[self.logg == k]) - dlogT ) for k in np.unique(self.logg)[::-1]]
        rightb += [ (rightb[-1][0] - dlogg, rightb[-1][1]) ]
        rightb  = [ (rightb[0][0] + dlogg, rightb[0][1]) ] + rightb
        b = leftb + rightb
        if closed:
            b += [b[0]]
        self._bound = (Path(np.array(b)), dlogT, dlogg)
        return self._bound[0]

    def genQ(self, qname, r, **kwargs):
        """ Generate a composite value from a previously calculated
            interpolation
            Works on 1 desired star or a population of stars

        Inputs:
            qname   quantity name from self.grid
            r   the result from a previous interpolation

        Outputs:
            an array containing the value
        """
        return ( self.grid[qname][r[:, 0].astype(int)] * r[:, 1] ).sum()

    def genSpectrum(self, T0, g0=None, Z0=None, weights=None, **kwargs):
        """ Generate a composite sprectrum
            Does the interpolation or uses a previously calculated
            interpolation
            Works on 1 desired star or a population of stars

        Inputs:
            T0  log(Teff) of each star or a 2d-array containing
                the result from a previous interpolation
            g0  log(g) of each stars
            Z0  metallicity

            if T0 and g0 are iterable, it calls interpMany

        Keywords:
            weights individual weights of each star
            **kwargs forwarded to interp(Many)

        Outputs:
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

    def gen_spectral_grid_from_given_points(self, pts, bounds=dict(dlogT=0.1, dlogg=0.3)):
        """ Reinterpolate a given stellar spectral library on to an Isochrone grid

        keywords
        --------
        pts: dict like structure of points
            dictionary like or named data structure of points to interpolate at.
            must contain:
                * logg  surface gravity in log-scale
                * logT  log of effective temperatures (in Kelvins)
                * logL  log of luminosity in Lsun units
                * Z     metallicity

        bounds:  dict
            sensitivity to extrapolation (see grid.get_stellib_boundaries)
            default: {dlogT:0.1, dlogg:0.3}

        Returns
        -------
        g: SpectralGrid
            Spectral grid (in memory) containing the requested list of stars and associated spectra
        """
        # Step 0: prepare outputs
        # =======================
        # Grid properties will be stored into a dictionary format until saved on disk
        # SEDs are kept into a ndarray
        ndata = len(pts)
        _grid  = {}
        _grid['radius'] = np.empty(ndata, dtype=float )
        _grid['keep'] = np.empty(ndata, dtype=bool )

        specs = np.empty( (ndata, len(self.wavelength)), dtype=float )

        # copy meta data of pts into the resulting structure
        for key in pts.keys():
            _grid[key] = np.asarray(pts[key])

        # Step 1: Avoid Extrapolation
        # ===========================
        # check boundary conditions, keep the data but do not compute the sed if not needed
        bound_cond = self.points_inside(zip(pts['logg'], pts['logT']))

        # Step 2: radii
        # =============
        # Stellar library models are given in cm^-2  ( 4 pi R)
        # Compute radii of each point using log(T) and log(L)
        radii = self.get_radius(pts['logL'], pts['logT'])
        rsun = ezunits.unit['Rsun'].to('m').magnitude  # 6.955e8 m
        _grid['radius'] = radii[:] / rsun

        # weights to apply during the interpolation
        # note that radii must be in cm
        weights = 4. * np.pi * (radii * 1e2) ** 2

        # Step 3: Interpolation
        # =====================
        # Do the actual interpolation, avoiding exptrapolations
        with progressbar.PBar(ndata, txt='spectral grid') as Pbar:
            for mk, rT, rg, rZ in enumerate(zip(pts['logT'], pts['logg'], pts['Z'])):
                Pbar.update(mk)
                if bound_cond[mk]:
                    s = np.array( self.interp(rT, rg, rZ, 0.) ).T
                    specs[mk, :] = self.genSpectrum(s) * weights[mk]

        # Step 4: filter points without spectrum
        # ======================================
        idx = np.array(bound_cond)

        lamb = self.wavelength[:]
        specs = specs.compress(idx, axis=0)
        for k in _grid.keys():
                _grid[k] = _grid[k].compress(idx, axis=0)

        # Step 5: Ship
        # ============
        header = {'stellib': self.source,
                  'comment': 'radius in Rsun',
                  'name': 'Reinterpolated stellib grid'}

        g = SpectralGrid(lamb, seds=specs, grid=Table(_grid), header=header, backend='memory')
        return g

    def plot_boundary(self, ax, **kwargs):
        """
        Valid kwargs are:
            agg_filter: unknown
            alpha: float or None
            animated: [True | False]
            antialiased or aa: [True | False]  or None for default
            axes: an :class:`~matplotlib.axes.Axes` instance
            clip_box: a :class:`matplotlib.transforms.Bbox` instance
            clip_on: [True | False]
            clip_path: [ (:class:`~matplotlib.path.Path`,         :class:`~matplotlib.transforms.Transform`) |         :class:`~matplotlib.patches.Patch` | None ]
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
        p = self.get_boundaries(dlogT=0.1, dlogg=0.3, closed=True)
        ax.add_patch(patches.PathPatch(**kwargs))
        return p

    def __add__(self, other):
        if not isNestedInstance(other, Stellib):
            raise ValueError('expecting a Stellib object, got {}'.format(type(other)))

        return CompositeStellib([self, other])


class CompositeStellib(Stellib):
    def __init__(self, osllist, *args, **kwargs):
        self._olist = osllist

    def __add__(self, other):
        if not isNestedInstance(other, Stellib):
            raise ValueError('expecting a Stellib object, got {}'.format(type(other)))

        lst = [k for k in self._olist] + [other]
        return CompositeStellib(lst)

    def __radd__(self, other):
        if not isNestedInstance(other, Stellib):
            raise ValueError('expecting a Stellib object, got {}'.format(type(other)))

        lst = [other] + [k for k in self._olist]
        return CompositeStellib(lst)

    @property
    def wavelength(self):
        """ return a common wavelength sampling to all libraries. This can be
        used to reinterpolate any spectrum onto a common definition """

        lambs = np.unique([ osl.wavelength for osl in self._olist ])
        return lambs

    @property
    def source(self):
        return ' + '.join([k.name for k in self._olist])

    def which_osl(self, xypoints, dlogT=0.1, dlogg=0.3):
        """
        Returns the first library indice that contains each point in xypoints

        keywords
        --------

        xypoints: sequence
            a sequence of N logg, logT pairs.

        dlogT: float
            margin in logT

        dlogg: float
            margin in logg

        returns
        -------
        r: ndarray(dtype=int)
            a ndarray, 0 meaning no library covers the point, and 1, ... n, for the n-th library
        """
        xy = np.asarray(xypoints)
        res = np.zeros(len(xy), dtype=int)
        res[self._olist[0].point_inside(xy, dlotT=dlogT, dlogg=dlogg)] = 1

        for ek, ok in enumerate(self._olist[1:]):
            if 0 in res:
                ind = np.squeeze(np.where(res == 0))
                r = ok.point_inside(xy[ind], dlotT=dlogT, dlogg=dlogg)
                res[ind[r]] = ek
        return r

    def __repr__(self):
        return "CompositeStellib, {}\n{}".format(object.__repr__(self), [k.name for k in self._olist])

    def get_boundaries(self, dlogT=0.1, dlogg=0.3, closed=True):
        """ Returns the closed boundary polygon around the stellar library with
        given margins

        keywords
        --------
        s: Stellib
            Stellar library object

        dlogT: float
            margin in logT

        dlogg: float
            margin in logg

        closed: bool
            if set, close the polygon

        returns
        -------
        b: ndarray[float, ndim=2]
            (closed) boundary points: [logg, Teff]

        Note
        ----
        as computing the boundary could take time, it is saved in the object
        and only recomputed when parameters are updated
        """
        if getattr(self, '_bound', None) is not None:
            if ((self._bound[1] - dlogT) < 1e-3) and (abs(self._bound[2] - dlogg) < 1e-3):
                return self._bound[0]

        b = [osl.get_boundaries(dlogT=dlogT, dlogg=dlogg, closed=closed) for osl in self._olist]
        self._bound = (Path.make_compound_path(*b), dlogT, dlogg)
        return self._bound[0]

    def interp(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6):
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
        (osl, r): tuple
            osl: is the library index starting from 1. 0 means no coverage.
            r: is the result from interp call on the corresponding library.
                a 3 to 12 star indexes and associated weights
        """
        osl_index = self.which_osl(np.atleast_2d([g0, T0]))[0]

        if osl_index > 0:
            return (osl_index, self._olist[osl_index - 1].interp(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6))
        else:
            return [(0, None)]

    def interpMany(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6, weights=None, pool=None, nthreads=__NTHREADS__):
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
        returns
        -------
        (osl, r): tuple
            osl: is the library index starting from 1. 0 means no coverage.
            r: is the result from interp call on the corresponding library.
                a 3 to 12 star indexes and associated weights
        """
        osl_index = self.which_osl(zip(g0, T0))

        g = []
        for oslk, osl in enumerate(self._olist):
            # make a generator to avoid keeping all in memory
            ind = np.where(osl_index - 1 == oslk)
            g.append( oslk + 1, osl.interpMany(T0[ind], g0[ind], Z0[ind], L0[ind], dT_max=dT_max, eps=eps, weights=weights, pool=pool, nthreads=nthreads) )

        return g

    def genQ(self, qname, r, **kwargs):
        """ Generate a composite value from a previously calculated
            interpolation
            Works on 1 desired star or a population of stars

        keywords
        --------
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

        keywords
        --------
        T0  ndarray(float)
            log(Teff) to obtain

        g0  ndarray(float)
            log(g) to obtain

        Z0 ndarray(float)
            metallicity values

        weights: ndarray(float)
            individual weights of each star

        **kwargs forwarded to interp(Many)

        Note: if T0 and g0 are iterable, it calls interpMany

        returns
        -------
        s: ndarray
            an array containing the composite spectrum reinterpolated onto self.wavelength
        """
        if Z0 is not None:
            if hasattr(T0, '__iter__'):
                _r = self.interpMany(T0, g0, Z0, weights=weights, **kwargs)
            else:
                _r = np.asarray(self.interp(T0, g0, Z0, weight=weights, **kwargs))
        else:
            _r = T0

        l0 = self.wavelength
        s = np.zeros(len(l0), dtype=float)
        for osl, _rk in _r:
            if osl > 0:
                sp = (((self._olist[osl - 1].spectra[_r[:, 0].astype(int)].T) * _r[:, 1]) ).sum(1)
                lamb = self._olist[osl - 1].wavelength
                s += np.interp(l0, lamb, sp)
        return s

    def gen_spectral_grid_from_given_points(self, pts, bounds=dict(dlogT=0.1, dlogg=0.3)):
        """ Reinterpolate a given stellar spectral library on to an Isochrone grid

        keywords
        --------
        pts: dict like structure of points
            dictionary like or named data structure of points to interpolate at.
            must contain:
                * logg  surface gravity in log-scale
                * logT  log of effective temperatures (in Kelvins)
                * logL  log of luminosity in Lsun units
                * Z     metallicity

        bounds:  dict
            sensitivity to extrapolation (see grid.get_stellib_boundaries)
            default: {dlogT:0.1, dlogg:0.3}

        Returns
        -------
        g: SpectralGrid
            Spectral grid (in memory) containing the requested list of stars and associated spectra
        """
        osl_index = self.which_osl(zip(pts['logg'], pts['logT']))

        seds = []
        grid = []
        l0 = self.wavelength
        for oslk, osl in enumerate(self._olist):
            # make a generator to avoid keeping all in memory
            _pts = ( (logg, logT) for (logg, logT, ok) in zip(pts['logg'], pts['logT'], osl_index) if (ok - 1 == oslk) )
            gk = osl.gen_spectral_grid_from_given_points(_pts, bounds=dict(dlogT=0.1, dlogg=0.3))
            gk.seds = interp1d(gk.lamb, gk.seds, axis=1)(l0)
            seds.append(gk.seds)
            grid.append(gk.grid)

        header = {'stellib': self.source,
                  'comment': 'radius in Rsun',
                  'name': 'Reinterpolated stellib grid'}

        _grid = recfunctions.stack_arrays( grid, defaults=None, usemask=False, asrecarray=True)

        g = SpectralGrid(l0, seds=np.vstack(seds), grid=Table(_grid), header=header, backend='memory')

        return g


class Elodie(Stellib):
    """ Elodie 3.1 stellar library derived class """
    def __init__(self, *args, **kwargs):
        self.name = 'ELODIE v3.1 (Prugniel et al 2007, astro-ph/703658)'
        self.source = __ROOT__ + '/libs/stellib_ELODIE_3.1.fits'
        self._load_()

    def _load_(self):
        with pyfits.open(self.source) as f:
            #load data
            self._getWaveLength_(f)
            self._getTGZ_(f)
            self._getSpectra_(f)

    def _getWaveLength_(self, f):
        self.wavelength = np.asarray((f[2].data[:]).tolist()).ravel()

    def _getTGZ_(self, f):
        cols = f[1].columns.names
        d = {}
        #d = { k: f[1].data.field(k) for k in cols }
        for k in cols:
            d[k] = f[1].data.field(k)
        self.grid = Table(d)
        self.grid.header['NAME'] = 'TGZ'
        del d, cols

    def _getSpectra_(self, f):
        self.spectra = f[0].data[:]

    @property
    def logg(self):
        return self.grid['logG']

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
    """ BaSeL 2.2 + Rauch stellar library derived class
        This library is used in Pegase.2
    """
    def __init__(self, *args, **kwargs):
        self.name = 'BaSeL 2.2 + Rauch (Pegase.2 version)'
        self.source = __ROOT__ + '/libs/stellib_BaSeL_2.2_Rauch.fits'
        self._load_()

    def _load_(self):
        with pyfits.open(self.source) as f:
            #load data
            self._getWaveLength_(f)
            self._getTGZ_(f)
            self._getSpectra_(f)

    def _getWaveLength_(self, f):
        self.wavelength = np.asarray((f[2].data[:]).tolist()).ravel()

    def _getTGZ_(self, f):
        cols = f[1].columns.names
        d = {}
        #d = { k: f[1].data.field(k) for k in cols }
        for k in cols:
            d[k] = f[1].data.field(k)
        self.grid = Table(d)
        self.grid.header['NAME'] = 'TGZ'
        del d, cols

    def _getSpectra_(self, f):
        self.spectra = f[0].data[:]

    @property
    def logg(self):
        return self.grid['logG']

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


class Kurucz(Stellib):
    """
    The stellar atmosphere models by Castelli and Kurucz 2004
    """
    def __init__(self, *args, **kwargs):
        self.name = 'Kurucz 2004'
        self.source = __ROOT__ + '/libs/kurucz2004.grid.fits'
        self._load_()

    def _load_(self):
        g = SpectralGrid(self.source, backend='memory')
        self.wavelength = g.lamb
        self.grid = g.grid
        self.grid.header['NAME'] = 'TGZ'
        self.spectra = g.seds

    @property
    def logT(self):
        return np.log10(self.grid['T0'])

    @property
    def logg(self):
        return self.grid['logG']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']


class Tlusty(Stellib):
    """
    Tlusty O and B stellar atmospheres
    """
    def __init__(self, *args, **kwargs):
        self.name = 'Tlusty'
        self.source = __ROOT__ + '/libs/tlusty.grid.fits'
        self._load_()

    def _load_(self):
        g = SpectralGrid(self.source, backend='memory')
        self.wavelength = g.lamb
        self.grid = g.grid
        self.grid.header['NAME'] = 'tlusty'
        self.spectra = g.seds

    @property
    def logT(self):
        return np.log10(self.grid['Teff'])

    @property
    def logg(self):
        return self.grid['logG']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']
