""" All to make simple figures

TODO: make a prior estimate function somewhere
TODO: make grid object more memory friendly,
        especially we deeply optimize reading memory with HDF storage using below classes
"""
import numpy as np


from beast.proba import getNorm_lnP as _getNorm_lnP
from beast.core.grid import SpectralGrid, ModelGrid
from beast.core.odict import odict
from beast.external.eztables import Table
from beast.proba import expectation, percentile
from beast.core.hdfstore import HDFStore


class HDF_Sparse_Storage(HDFStore):
    """ Handling in and out of the HDF5 file of sparse lnp format
        This class can be used as a context manager

        Any attribute of the HDF will be directly available transparently if
        the source if opened
    """

    def get_star(self, obj, inputs=False):
        """get_star  -- quick access to a star in the sparse lnp storage

        keywords
        --------

        obj: int
            object id in the file

        inputs: bool (optional)
            if set, add input values to the returned arguments

        returns
        -------

        indx: ndarray[ndim=1, dtype=int]
            model coordinates from the sparse lnp storage

        lnps: ndarray[ndim=1, dtype=float]
            lnp associated to model coordinates from the sparse lnp storage

        inpt: ndarray[ndim=2, dtype=float]
            only if input argument is set
            input table: flux, errp, errm, mask
        """
        if inputs is True:
            with self as s:
                indx = s.getNode('/star_{:d}/idx'.format(obj)).read().astype(int)
                lnps = s.getNode('/star_{:d}/lnp'.format(obj)).read().astype(float)
                inpt = s.getNode('/star_{:d}/input'.format(obj)).read()
                return indx, lnps, inpt
        else:
            with self as s:
                indx = s.getNode('/star_{:d}/idx'.format(obj)).read().astype(int)
                lnps = s.getNode('/star_{:d}/lnp'.format(obj)).read().astype(float)
            return indx, lnps

    def set_star(self, obj, inputs, lnp, lnp_threshold=-40, flush=True):
        """set_star -- write a star into the storage

        keywords
        --------

        obj: int
            object id number

        inputs: ndarray[ndim=2, dtype=float]
            input information stored as is, e.g., np.array([sed, errp, errm, mask]).T

        lnp: ndarray[ndim=1, dtype=float]
            lnp values computed from the likelihood

        lnp_threshold: float (default -40)
            toss out grid points where lnp - lnp_max > threshold (lnps < 0)

        flush: bool (default True)
            commit to the storage file if set

        """
        #Need ragged arrays rather than uniform table
        if self._mode[0] == 'r':
            raise AttributeError('Storage opened in read-only mode')

        with self as storage:
            # manual attack avoids multiple flush calls, ie., faster this way!
            star_group = storage.createGroup('/', 'star_%d'  % obj, title="star %d" % obj)
            indx = np.squeeze(np.where((lnp - max(lnp[np.isfinite(lnp)])) > lnp_threshold))
            storage.createArray(star_group, 'input', np.asarray(inputs) )
            storage.createArray(star_group, 'idx', np.array(indx, dtype=np.int32))
            storage.createArray(star_group, 'lnp', np.array(lnp[indx], dtype=np.float32))
            #commit changes
            if flush is True:
                self.flush()


class PosteriorResult(object):
    """posterior of a single star that can be evaluated over parameters"""
    def __init__(self, indx, lnps, inputs, weights, proxy):
        """__init__

        keywords
        --------

        indx: ndarray[ndim=1, dtype=int]
            model coordinates from the sparse lnp storage

        lnps: ndarray[ndim=1, dtype=float]
            lnp associated to model coordinates from the sparse lnp storage

        inputs: ndarray[ndim=2, dtype=float]
            input information stored as is, e.g., np.array([sed, errp, errm, mask]).T

        weights: ndarray[ndim=1, dtype=float]
            weights associated to model coordinates, ie, priors
            (default uniform weights)

        proxy: PosteriorProxy
            allow extraction of model parameters to evaluate the posterior on
        """
        # check arguments
        if len(indx) != len(lnps):
            raise AttributeError('indx and lnps dimensions mismatched')

        # make default weights
        if weights is None:
            weights = np.ones(len(indx), dtype=float)

        self.indx = indx
        self.lnps = lnps
        self.weights = weights
        self.log_norm = self.__lognorm__()
        self.proxy = proxy
        self.inputs = inputs

    @property
    def norm(self):
        return np.exp(self.log_norm)

    def __lognorm__(self):
        """__lognorm__ -- normalization constant of the lnps

        note: uses beast's internal getNorm_lnp optimized to deal with the sum
        of exponentials.

        returns
        -------
        lognorm: float
            log of the normalization constant of the lnps
            if this number ends up non-finite, returns max(lnps)
        """
        log_norm = np.log(_getNorm_lnP(self.lnps))
        if not np.isfinite(log_norm):
            log_norm = self.lnps.max()
        return log_norm

    def evaluate(self, qname):
        """evaluate -- returns the ND coordinates and posterior values (not log-posterior)

        keywords
        --------

        qname: str or list of str
            if str:  name of the quantity or expression to evaluate from the grid table
            if list: list of qquantities or expresions

        returns
        -------
        q: ndarray like
            array of qname at which the sparse likelihood is defined

        weights: ndarray[dtype=float, ndim=1]
            posterior probability values for each q point
        """
        g0 = self.proxy.models
        if hasattr(qname, '__iter__'):
            q = Table(g0.grid[self.indx])
            q = np.asarray([ q.evalexpr(qk) for qk in qname])
        else:
            q = g0.grid.evalexpr(qname)[self.indx]

        weights = np.exp(self.lnps - self.log_norm)
        if self.weights is not None:
            weights *= self.weights

        return q, weights

    def Q_best(self, qname=None):
        """ Best Property values:

        keywords
        --------

        qname: str or list of str or None
            if str:  name of the quantity or expression to evaluate from the grid table
            if list: list of quantities or expresions
            if None: all properties are returned

        returns
        -------
        e_dict: dict or ndarray[float, ndim=1]
            if qname is iterable, returns a dict with a (qname, ndarray) pairs
            else returns only the ndarray of best values (one per obj in cllist)
        """
        #get quantities
        g0 = self.proxy.models
        lnps = self.lnps.astype(float)
        indx = self.indx.astype(int)
        weights = np.exp(lnps - self.log_norm)
        if self.weights is not None:
            weights *= self.weights
        data = g0.grid[indx[weights.argmax()]]
        if qname is None:
            qname = g0.grid.keys()
        if hasattr(qname, '__iter__'):
            r = odict()
            for qk in qname:
                lbl = '{:s}_Best'.format(qk)
                r[lbl] = data[qk]
        else:
            r = data[qname]
        return r

    def Q_percentile(self, qname=None, p=[16., 50., 84.]):
        """ Percentile values of any given grid property (incl. expression) but seds,

        keywords
        --------

        qname: str or list of str or None
            if str:  name of the quantity or expression to evaluate from the grid table
            if list: list of quantities or expresions
            if None: all properties are returned

        p: array-like
            list of percentile values

        OUTPUT
        ------
        e_dict: dict or ndarray[float, ndim=1]
            if qname is iterable, returns a dict with a (qname, ndarray) pair
            else returns only the ndarray of percentile values
        """
        #get quantities
        g0 = self.proxy.models
        lnps = self.lnps.astype(float)
        indx = self.indx.astype(int)
        weights = np.exp(lnps - self.log_norm)
        if self.weights is not None:
            weights *= self.weights
        _p = np.asarray(p, dtype=float)
        data = g0.grid[indx]
        if qname is None:
            qname = g0.grid.keys()

        if hasattr(qname, '__iter__'):
            r = odict()
            for qk in qname:
                tmp = percentile(data[qk], _p, weights=weights)
                for ek, pk in enumerate(p):
                    lbl = '{:s}_p{:d}'.format(qk, int(pk))
                    r[lbl] = tmp[ek]
        else:
            r = percentile(data[qname], _p, weights=weights)
        return r

    def Q_expect(self, qname=None):
        """ Expectation values of any given grid property (incl. expression) but seds,
        i.e.:
                integral(p(q) * q dq) / integral(p(q) dq),
        which in a discrete world becomes
                sum(p(q_i) * q_i) / sum(p(q_i)

        see sed_expect for sed expectation values

        keywords
        --------

        qname: str or list of str or None
            if str:  name of the quantity or expression to evaluate from the grid table
            if list: list of quantities or expresions
            if None: all properties are returned

        returns
        -------
        e_dict: dict or ndarray[float, ndim=1]
            if qname is iterable, returns a dict with a (qname, ndarray) pair
            else returns only the ndarray of expected values (one per obj in objlist)
        """
        g0 = self.proxy.models
        lnps = self.lnps.astype(float)
        indx = self.indx.astype(int)
        weights = np.exp(lnps - self.log_norm)
        if self.weights is not None:
            weights *= self.weights
        data = g0.grid[indx]
        if qname is None:
            qname = g0.grid.keys()

        #get quantities
        if hasattr(qname, '__iter__'):
            r = odict()
            for qk in qname:
                lbl = '{:s}_E'.format(qk)
                r[lbl] = expectation(data[qname], weights=weights)
        else:
            r = expectation(data[qname], weights=weights)

        return r

    def get_2dmap(self, q1, q2, bins=100, log=True, smooth=0, order=1, *args, **kwargs):
        """get_2dmap - grid the (q1, q2) joined posterior

        keywords
        --------
        q1: str
            first axis quantity. Can be an expression

        q2: str
            second axis quantity. Can be an expression

        bins: int, tuple of lists
            grid definition based on np.histogram defition

        log: bool
            apply log-scale transformation

        smooth: float
            apply a smoothing factor in pixel units (0 = none)

        order: int
            convolution order to use in the smoothing

        returns
        -------
        zi: ndarray[dtype=float, ndim=2]
            image map

        e: list of 4 floats
            extent of the map
        """

        from scipy.ndimage.filters import gaussian_filter
        from scipy.ndimage import zoom
        from matplotlib.mlab import griddata

        if bins:
            if np.shape(bins) == ():
                xbin = ybin = bins
            else:
                assert(np.size(bins) == 2)
                xbin, ybin = bins

        (x, y), z = self.evaluate([q1, q2])

        if log is True:
            z = np.log10(z)

        xi = np.linspace(x.min(), x.max(), xbin)
        yi = np.linspace(y.min(), y.max(), ybin)
        zi = griddata(x, y, z, xi, yi)

        extent = [xi.min(), xi.max(), yi.min(), yi.max()]
        if smooth > 0:
            zi = zoom(zi, smooth, order=order)
            zi = gaussian_filter(zi, smooth / 5.)

        return zi, extent

    def scatter_map(self, q1, q2, log=True, ax=None, **kwargs):
        """get_2dmap - grid the (q1, q2) joined posterior

        keywords
        --------
        q1: str
            first axis quantity. Can be an expression

        q2: str
            second axis quantity. Can be an expression

        log: bool
            apply log-scale transformation

        ax: matplotlib.axes
            axis to plot into (default: plt.gca())

        **kwars: forwarded to plt.scatter

        returns
        -------
        sc: plt.scatter result
        """
        import pylab as plt

        if ax is None:
            ax = plt.gca()

        (x, y), z = self.evaluate([q1, q2])
        if log is True:
            z = np.log10(z)

        ax.set_xlabel('$' + q1 + '$')
        ax.set_ylabel('$' + q2 + '$')
        return ax.scatter(x, y, c=z, **kwargs)

    def plot_contour_map(self, q1, q2, bins=100, levels=None, log=False, smooth=0, order=1, ax=None, **kwargs):
        """plot_contour_map

        keywords
        --------
        q1: str
            first axis quantity. Can be an expression

        q2: str
            second axis quantity. Can be an expression

        bins: int, tuple of lists
            grid definition based on np.histogram defition

        levels: iterable
            levels to use in the contours in quartile units
            (default [0.95, 0.75, 0.5, 0.25, 0.05])

        log: bool
            apply log-scale transformation

        smooth: float
            apply a smoothing factor in pixel units (0 = none)

        order: int
            convolution order to use in the smoothing

        ax: matplotlib.axes
            axis to plot into (default: plt.gca())

        **kwargs forwarded to plt.contour

        returns
        -------
        cs: matplotlib.contour.QuadContourSet
            result from plt.contour
        """
        import pylab as plt

        if ax is None:
            ax = plt.gca()

        if levels is None:
            levels = [0.95, 0.75, 0.5, 0.25, 0.05]

        h, e = self.get_2dmap(q1, q2, bins=bins, log=log, smooth=smooth, order=order)
        ax.set_xlabel('$' + q1 + '$')
        ax.set_ylabel('$' + q2 + '$')
        return ax.contour(h / h.max(), levels=levels, extent=e, **kwargs)

    def imshow_map(self, q1, q2, bins=100, log=True, smooth=0, order=1, origin='lower', aspect='auto', ax=None, **kwargs):
        """imshow_map

        keywords
        --------
        q1: str
            first axis quantity. Can be an expression

        q2: str
            second axis quantity. Can be an expression

        bins: int, tuple of lists
            grid definition based on np.histogram defition

        log: bool
            apply log-scale transformation

        smooth: float
            apply a smoothing factor in pixel units (0 = none)

        order: int
            convolution order to use in the smoothing

        ax: matplotlib.axes
            axis to plot into (default: plt.gca())

        **kwargs forwarded to plt.contour

        returns
        -------
        im: matplotlib.image.AxesImage
            result from plt.imshow
        """
        import pylab as plt

        if ax is None:
            ax = plt.gca()

        h, e = self.get_2dmap(q1, q2, bins=bins, log=log, smooth=smooth, order=order)
        ax.set_xlabel('$' + q1 + '$')
        ax.set_ylabel('$' + q2 + '$')
        return ax.imshow(h, extent=e, origin=origin, aspect=aspect, **kwargs)


class PosteriorProxy(object):
    """ Layer on top of a storage
        Priors will be eventually allowed to be set from either individual
        weights or from functions.
    """
    def __init__(self, lnps, models):
        """__init__

        keywords
        --------

        lnps: str or HDF_Sparse_Storage
            where sparse log-likelihood are stored

        models: str or beast.core.grid.Grid
            grid of models
        """
        if isinstance(lnps, HDF_Sparse_Storage):
            self.lnps = lnps
        elif type(lnps) == str:
            self.lnps = HDF_Sparse_Storage(lnps)
        else:
            raise TypeError('lnps is expected to be str or HDF_Sparse_Storage, got {}'.format(type(lnps)))

        if isinstance(models, ModelGrid) or isinstance(models, SpectralGrid):
            self.models = models
        elif type(models) == str:
            self.models = HDF_Sparse_Storage(models)
        else:
            raise TypeError('models is expected to be str or FileSEDGrid, got {}'.format(type(models)))

        self.priors = None
        self.pvalues = None

    def set_priors(self, pvalues):
        """set_priors -- set the priors from individual model weights

        keywords
        --------

        pvalues: ndarray[ndim=1, dtype=float]
            individual weights to apply to each models.
            assumes len(pvalues) == len(models)
        """
        self.pvalues = pvalues

    def get_posterior(self, obj):
        """get_posterior -- generate a PosteriorResult Object from which all
        PDF can be generated for a given star

        keywords
        --------

        obj: int
            object id in the file

        returns
        -------
        p: PosteriorResult
            PosteriorResult object of the requested star
        """
        indx, lnps, inpt = self.lnps.get_star(obj, inputs=True)
        if self.pvalues is not None:
            prior = self.pvalues[indx]
        else:
            prior = None
        return PosteriorResult(indx, lnps, inpt, prior, self)
