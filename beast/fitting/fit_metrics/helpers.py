"""
HDF5 Storage optimized functions
================================

    arange                      Return evenly spaced values within a given interval.
    best_bins                   get non-uniform binning of unique values of Q.
    compute_uniform_prior       Compute the individual model weights to make a given flat prior.
    get_Q_from_node             returns a quantity from a HDF5 Node given its math expression.
    get_centers_from_bins       returns the bin centers from a list of edges.
    get_nclusters               returns the number of likelihoods stored in a file.
    nice_bins                   Define a prior bin width on a quantity q and
                                    iterate until none of the bins are purely empty.
    Q_expect                    Expectation values of any given grid property Q (incl. expression) but seds
    Q_percentile                Percentile values of any given grid property Q (incl. expression) but seds
    sed_expect                  Sed expectation values
    sed_percentile              SED percentile values

:author: MF
:last update: Fri Jun 14 17:35:21 PDT 2013
"""
import numpy as np
import tables

#from .progressbar import PBar
from .likelihood import getNorm_lnP
from .common import percentile, expectation

__all__ = [ 'arange', 'best_bins', 'compute_uniform_prior',
    'get_Q_from_node', 'get_centers_from_bins', 'get_nclusters',
    'nice_bins', 'Q_expect', 'sed_expect', 'Q_percentile', 'sed_percentile']


def get_Q_from_node(node, expr, condvars={}):
    """ returns a quantity from a HDF5 Node given its math expression.
        Assuming that all quantities are either from the node of in condvars

        all np function can be used (log, exp, pi...)

        method
        ------
        Pure eval expression was too slow when using disk access and HD5 nodes.
        Instead, we use the tables.expression feature that parse "expr" and
        extracts the variable names for interest.
        Based on this list we build the context of evaluation by add missing
        values from the node.

        INPUTS
        ------
        node: tables.table.Table
            the Table node to get data from (need named columns)

        expr: str
            the mathematical expression which could include numpy functions or
            simply be a column name

        KEYWORDS
        --------
        condvars: dict
            The condvars mapping may be used to define the variable names
            appearing in the condition.

        OUTPUTS
        -------
        q: ndarray like
            the column values requested with type and shape defined in the table.
    """
    colnames = node.colnames
    if expr in colnames:
        q = node.read(field=expr)
    else:
        #parse qname
        names = tables.expression.getExprNames(expr, {})[0]
        for k in names:
            if k in colnames:
                condvars[k] = node.read(field=k)
        q = eval(expr, condvars, np.__dict__)
    return q


def get_centers_from_bins(bins):
    """returns the bin centers
    INTPUTS
    -------
    bins: ndarray[float, ndim=1]
        the bin edges to get the centers from

    OUTPUTS
    -------
    center: ndarray[float, ndim=1]
        center of the bins
    """
    return 0.5 * (bins[:-1] + bins[1:])


def arange(v0, vn, dv, inflate=0.0):
    """
    Return evenly spaced values within a given interval.
    A more robust version than numpy.arange
    When using a non-integer step, such as 0.1, the results from numpy are
    often not consistent.  Instead we use linspace as a correction here.

    INPUTS
    ------
    v0: float
        Start of interval.

    vn: float
        Stop of interval.

    dv: float
        Initial spacing between values. The final value could be slighly larger
        when inflate is set to non-zero value

    KEYWORDS
    --------
    inflate: float
        margin in dv units to keep around the data

    OUTPUTS
    -------
    val: ndarray[float, ndim=1]
        returns evenly spaced values
    """
    if inflate > 0.:
        tmp = inflate * dv
        nb = int((vn - v0 + 2 * tmp) / dv + 1.)
        return np.linspace(v0, vn + tmp, nb) - 0.5 * tmp
    else:
        nb = int((vn - v0) / dv + 1.)
        return np.linspace(v0, vn, nb)


def nice_bins(q, dq0, ddq=0.1, centers=False):
    """ Define a prior bin width on a quantity q and iterate until none of the bins are purely empty
    This is not optimal binning, but only make sure that the final bin width is
    big enough that we never exhibit holes in the unique values of Q.

    INPUTS
    ------
    q: ndarray[float, ndim=1]
        the quantity you want to bin
    dq0: float
        initial bin width

    KEYWORDS
    --------
    ddq: float
        fractional variation to explore at each iteration step

    OUTPUTS
    -------
    bins: ndarray[float, ndim=1]
        final uniformely spaced bins
    dq: float
        final bin width
    """
    val = np.unique(q)
    b0 = 1.5 * val[0] - 0.5 * val[1]
    bn = 1.5 * val[-1] - 0.5 * val[-2]
    dq = dq0

    b = arange(b0, bn, dq)
    h, _ = np.histogram(val, b)
    if 0 in h:
        _ddq = ddq * dq
        while 0 in h:
            dq += _ddq
            b = arange(b0, bn, dq)
            h, _ = np.histogram(val, b)
    if centers:
        return b, dq, get_centers_from_bins(b)
    else:
        return b, dq


def best_bins(Q, centers=False):
    """ get non-uniform binning of unique values of Q
    INPUTS
    ------
    Q: ndarray[float, ndim=1]
        the list of values to bin

    KEYWORDS
    --------
    centers: bool
        if set returns the centers as well

    OUTPUTS
    -------
    if centers is set
        bins, centers: tuple of ndarray[float, ndim=1]
            selected bin edges and centers
    else
        bins: ndarray[float, ndim=1]
            bin edges

    """
    val = np.unique(Q)
    bins = np.hstack([[1.5 * val[0] - 0.5 * val[1]], val[:-1] + 0.5 * np.diff(val), [1.5 * val[-1] - 0.5 * val[-2]]])
    if not centers:
        return bins
    else:
        return bins, get_centers_from_bins(bins)


def compute_uniform_prior(node, qnames=[], method='nice', correlation=True, dq=0.1):
    """ Compute the corrected weights of each model on the mods_grid node in
    order to obtain a uniform prior.  The prior can be a combination of N
    parameters, correlated or not.

    The method implemented is based on density estimations in N x 1d for
    uncorrelated dimensions or N-d for correlated dimensions.

    When method is set to nice:
        The density estimation makes use of binning the models into evenly
        space bins on each dimensions and the binning is defined such that none
        of the bins are empty
        see also: nice_bins

    When method is set to best:
        The bins are made from individual unique values in each dimension.
        Hence the resulting partitions are non even.

    Once the bins are defined, the final N-d volumes are made by the cartesian
    product of the dimension bins.  The number of models in each bins (could be
    0 in N-d space) gives the non-normalized intrisic density distribution of
    the models. The uniform prior correction is then the inverse of the density
    for each model.

    Note: the optimal binning is not trivial.
        In fact the best way to correct for a non-uniform prior is do random
        samples of the prior distribution and weight the models accordingly.
        But this method is expensive for high number of models and complex
        intrinsic priors.

        the nice binning method is useful for plots because it avoids most of
        visual artifacts while not being computed for each plot.

    INPUTS
    ------
    node: tables.table.Table
        the mods_grid table node containing the grid point parameters.

    qnames: list
        the property expression set that defines the prior

    KEYWORDS
    --------
    method: 'nice' or 'best'
        the binning method. See nice_bins, best_bins

    correlation: bool
        if set the prior is estimated in the N dimensional space define by the qnames set.

    dq: float
        used in the binning method when set to 'nice'

    OUTPUT
    ------
    prior: ndarray(float, ndim=1)
        the normalized probability of each model to apply in order to have the
        requested uniform prior.
    """
    if (method not in ['best', 'nice']):
        raise ValueError('Method must be in [best, nice]')

    if type(qnames) == str:
        qnames = list(qnames)

    if len(qnames) == 0:
        return 1.

    if correlation is False or len(qnames) == 1:
        final = np.ones(node.shape, dtype=float)
        for qk in qnames:
            q = get_Q_from_node(node, qk)
            if method.lower() == 'nice':
                bins, _ = nice_bins(q, dq, centers=False)
            else:
                bins = best_bins(q, centers=False)
            # uniform prior correction
            prior, _ = np.histogram(q, bins=bins)
            prior = prior.astype(float)
            prior[ prior == 0 ] = 1.
            prior /= prior.sum()
            final *= prior[np.digitize(q, bins) - 1]
    else:   # correlated priors
        sample = [ get_Q_from_node(node, qk) for qk in qnames ]
        if method.lower() == 'nice':
            bins = [ nice_bins(qk, dq, centers=False)[0] for qk in sample ]
        else:
            bins = [ best_bins(qk, centers=False) for qk in sample ]

        dd = np.vstack([(np.digitize(sk, bk) - 1) for (sk, bk) in zip(sample, bins)]).T
        prior, _ = np.histogramdd(np.asarray(sample).T, bins=bins)
        final = np.asarray( [prior[tuple(vk)] for vk in dd] )

    final = final.astype(float)
    ind = final > 0
    final[ind] = 1. / final[ind]
    return final


def get_nclusters(lnpfile):
    """ Return the number of clusters in lnpfile """
    if type(lnpfile) == str:
        f = tables.openFile(lnpfile)
    elif isinstance(lnpfile, tables.file.File):
        f = lnpfile

    nobs = f.root.flux.shape[1]

    if not isinstance(lnpfile, tables.file.File):
        f.close()

    return nobs


def Q_expect(lnpfile, qname, cllist=None, prior=None):
    """ Expectation values of any given grid property (incl. expression) but seds,
    i.e.:
            integral(p(q) * q dq) / integral(p(q) dq),
    which in a discrete world becomes
            sum(p(q_i) * q_i) / sum(p(q_i)

    see sed_expect for sed expectation values

    INPUTS
    ------

    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    KEYWORDS
    --------

    cllist: list or array like
        index numbers of clusters to extract

    prior: ndarray[float, ndim=1]
        prior probability of each model on the grid to apply (default 1/Nmodel)
        see: compute_uniform_prior

    OUTPUT
    ------
    e_dict: dict or ndarray[float, ndim=1]
        if qname is iterable, returns a dict with a (qname, ndarray) pair
        else returns only the ndarray of expected values (one per obj in cllist)
    """
    if type(lnpfile) == str:
        f = tables.openFile(lnpfile)
    elif isinstance(lnpfile, tables.file.File):
        f = lnpfile

    if cllist is None:
        nobs = get_nclusters(f)
        cllist = range(nobs)
    else:
        nobs = len(cllist)

    #get quantities
    if hasattr(qname, '__iter__'):
        r = {}
        for qk in qname:
            r[qk] = Q_expect(f, qk, cllist, prior)
    else:
        #get grid node
        node = f.getNode('/mods_grid')
        q = get_Q_from_node(node, qname)
        r = np.empty(len(cllist), dtype=float)
        with PBar(nobs, txt='Expectations') as pb:
            for e, clk in enumerate(cllist):
                pb.update(e, txt='E({})'.format(qname))
                clnode = f.getNode('/lnps_cl%d' % clk)
                indx = clnode[:, 0].astype(int)
                lnps = clnode[:, 1].astype(float)
                log_norm = np.log(getNorm_lnP(lnps))
                if not np.isfinite(log_norm):
                    log_norm = lnps.max()
                weights = np.exp(lnps - log_norm)
                if prior is not None:
                    weights *= prior[indx]
                r[e] = expectation(q[indx], weights=weights)

    if not isinstance(lnpfile, tables.file.File):
        f.close()

    return r


def Q_percentile(lnpfile, qname, p=[16., 50., 84.], cllist=None, prior=None):
    """ Percentile values of any given grid property (incl. expression) but seds,

    see also:    sed_percentile, percentile

    INPUTS
    ------

    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    KEYWORDS
    --------
    p: array-like
        list of percentile values

    cllist: list or array like
        index numbers of clusters to extract

    prior: ndarray[float, ndim=1]
        prior probability of each model on the grid to apply (default 1/Nmodel)
        see: compute_uniform_prior

    OUTPUT
    ------
    e_dict: dict or ndarray[float, ndim=2]
        if qname is iterable, returns a dict with a (qname, ndarray) pair
        else returns only the ndarray of percentile values (one per obj in cllist)
    """
    if type(lnpfile) == str:
        f = tables.openFile(lnpfile)
    elif isinstance(lnpfile, tables.file.File):
        f = lnpfile

    if cllist is None:
        nobs = get_nclusters(f)
        cllist = range(nobs)
    else:
        nobs = len(cllist)

    #get quantities
    if hasattr(qname, '__iter__'):
        r = {}
        for qk in qname:
            r[qk] = Q_percentile(f, qk, p, cllist, prior)
    else:
        #get grid node
        node = f.getNode('/mods_grid')
        nval = len(p)
        _p = np.asarray(p, dtype=float)
        q = get_Q_from_node(node, qname)
        r = np.empty((len(cllist), nval), dtype=float)
        with PBar(nobs, txt='Percentiles') as pb:
            for e, clk in enumerate(cllist):
                pb.update(e, txt='Percentiles({})'.format(qname))
                clnode = f.getNode('/lnps_cl%d' % clk)
                indx = clnode[:, 0].astype(int)
                lnps = clnode[:, 1].astype(float)
                log_norm = np.log(getNorm_lnP(lnps))
                if not np.isfinite(log_norm):
                    log_norm = lnps.max()
                weights = np.exp(lnps - log_norm)
                if prior is not None:
                    weights *= prior[indx]
                r[e, :] = percentile(q[indx], _p, weights=weights)

    if not isinstance(lnpfile, tables.file.File):
        f.close()

    return r


def sed_expect(lnpfile, cllist=None, prior=None):
    """ Sed Expectation values, i.e.:
            integral(p(sed) * sed dsed) / integral(p(sed) dsed),
    which in a discrete world becomes
            sum(p(sed_i) * sed_i) / sum(p(sed_i)

    see Q_expect for other property expectation values

    INPUTS
    ------

    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    KEYWORDS
    --------

    cllist: list or array like
        index numbers of clusters to extract

    prior: ndarray[float, ndim=1]
        prior probability of each model on the grid to apply (default 1/Nmodel)
        see: compute_uniform_prior

    OUTPUT
    ------
    e_sed: ndarray[float, ndim=2]
        Nobs x Nfilters array of expected seds
    """
    if type(lnpfile) == str:
        f = tables.openFile(lnpfile)
    elif isinstance(lnpfile, tables.file.File):
        f = lnpfile

    #get grid node
    node = f.getNode('/mods_seds')
    nfilters = node.shape[1]
    #names = f.root.mods_grid.attrs['filters']

    if cllist is None:
        nobs = get_nclusters(f)
        cllist = range(nobs)
    else:
        nobs = len(cllist)

    #get quantities
    r = np.empty((nobs, nfilters), dtype=float)
    with PBar(nfilters, txt='SED Expectations') as pb:
        for fk in range(nfilters):
            q = node[:, fk]
            for e, clk in enumerate(cllist):
                pb.update(e, txt='E(filter {})'.format(fk))
                clnode = f.getNode('/lnps_cl%d' % clk)
                indx = clnode[:, 0].astype(int)
                lnps = clnode[:, 1].astype(float)
                log_norm = np.log(getNorm_lnP(lnps))
                if not np.isfinite(log_norm):
                    log_norm = lnps.max()
                weights = np.exp(lnps - log_norm)
                if prior is not None:
                    weights *= prior[indx]
                r[e, fk] = expectation(q[indx], weights=weights)

    if not isinstance(lnpfile, tables.file.File):
        f.close()

    return r


def sed_percentile(lnpfile, p=[16., 50., 84.], cllist=None, prior=None):
    """ Sed percentile values:

    see Q_percentile for other property expectation values

    INPUTS
    ------

    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    KEYWORDS
    --------
    p: array-like
        list of percentile values

    cllist: list or array like
        index numbers of clusters to extract

    prior: ndarray[float, ndim=1]
        prior probability of each model on the grid to apply (default 1/Nmodel)
        see: compute_uniform_prior

    OUTPUT
    ------
    e_sed: ndarray[float, ndim=3]
        N(obs) x N(filters) x N(p) array of percentile seds
    """
    if type(lnpfile) == str:
        f = tables.openFile(lnpfile)
    elif isinstance(lnpfile, tables.file.File):
        f = lnpfile

    #get grid node
    node = f.getNode('/mods_seds')
    nfilters = node.shape[1]
    #names = f.root.mods_grid.attrs['filters']

    if cllist is None:
        nobs = get_nclusters(f)
        cllist = range(nobs)
    else:
        nobs = len(cllist)

    nval = len(p)
    _p = np.asarray(p, dtype=float)

    #get quantities
    r = np.empty((nobs, nfilters, nval), dtype=float)
    with PBar(nobs, txt='SED percentile') as pb:
        for fk in range(nfilters):
            q = node[:, fk]
            for e, clk in enumerate(cllist):
                pb.update(e, txt='Percentiles(filter {})'.format(fk))
                clnode = f.getNode('/lnps_cl%d' % clk)
                indx = clnode[:, 0].astype(int)
                lnps = clnode[:, 1].astype(float)
                log_norm = np.log(getNorm_lnP(lnps))
                if not np.isfinite(log_norm):
                    log_norm = lnps.max()
                weights = np.exp(lnps - log_norm)
                if prior is not None:
                    weights *= prior[indx]
                r[e, fk, :] = percentile(q[indx], _p, weights=weights)

    if not isinstance(lnpfile, tables.file.File):
        f.close()

    return r


def Q_best(lnpfile, qname, cllist=None, prior=None):
    """ Best Property values:

    Note: external loop on Q whereas cllist to extract data only once!

    Note on parallel version:
        # from ised.external.ezmap import map
        # def job(qk):
        #    return (qk, Q_best(lnpfile, qk, cllist, prior))
        # r.update( map(job, qname, ncpu=len(qname)) )

    INPUTS
    ------

    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    KEYWORDS
    --------

    cllist: list or array like
        index numbers of clusters to extract

    prior: ndarray[float, ndim=1]
        prior probability of each model on the grid to apply (default 1/Nmodel)
        see: compute_uniform_prior

    OUTPUT
    ------
    e_dict: dict or ndarray[float, ndim=1]
        if qname is iterable, returns a dict with a (qname, ndarray) pairs
        else returns only the ndarray of best values (one per obj in cllist)
    """
    if type(lnpfile) == str:
        f = tables.openFile(lnpfile)
    elif isinstance(lnpfile, tables.file.File):
        f = lnpfile

    if cllist is None:
        nobs = get_nclusters(f)
        cllist = range(nobs)
    else:
        nobs = len(cllist)

    #get quantities
    if hasattr(qname, '__iter__'):
        r = {}
        for qk in qname:
            r[qk] = Q_best(f, qk, cllist, prior)
    else:
        #get grid node
        node = f.getNode('/mods_grid')
        q = get_Q_from_node(node, qname)
        r = np.empty(len(cllist), dtype=float)
        with PBar(nobs, txt='Best') as pb:
            for e, clk in enumerate(cllist):
                pb.update(e, txt='Best({})'.format(qname))
                clnode = f.getNode('/lnps_cl%d' % clk)
                indx = clnode[:, 0].astype(int)
                lnps = clnode[:, 1].astype(float)
                log_norm = np.log(getNorm_lnP(lnps))
                if not np.isfinite(log_norm):
                    log_norm = lnps.max()
                weights = np.exp(lnps - log_norm)
                if prior is not None:
                    weights *= prior[indx]
                r[e] = q[indx[weights.argmax()]]

    if not isinstance(lnpfile, tables.file.File):
        f.close()

    return r

def calc_covar_mat(flux_in,flux_rec):
    """Computes the covariance matrix, the inverse of the Cholesky decomposition and lnQ for a given star.

    INPUTS
    -----
    flux_in: ndarray[float, ndim=1]
       input SED corresponding to a modeled SED
    flux_rec: ndarray[float, ndim=1]
       recovered SED from the Artificial Star Tests

    OUTPUTS
    -------
    cov_mat, inv_cholesky_decomposition, lnQ: Tuple (ndarray[float, ndim=2], ndarray[float, ndim=2], float)
    inv_cholesky_covar = C^-1 where cov_mat=CC^T (lower triangular matrix)
    ln(Q) = ln(sqrt(|covariance_matrix|))
    """
    
    diffs = flux_in - flux_rec
    diffs -= diffs.mean(axis=0)

    cov_mat = np.cov(diffs.T) 
    lnQ = -0.5*np.log(np.linalg.det(cov_mat))
    cholesky_decomposition = np.linalg.cholesky(cov_mat)
    inv_cholesky_decomposition = np.linalg.inv(cholesky_decomposition)
    return (cov_mat, inv_cholesky_decomposition, lnQ)
