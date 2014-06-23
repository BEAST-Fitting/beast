"""
Fitting Pipeline
=============

I use `ezpipe`, a pipeline package I wrote in order to clean the syntax and
allow more flexibilities. In particular it will simplifies the management of
intermediate results or broken jobs.

do the fit
----------
the pipeline is sequence of tasks
    tasks = ( t_fit(g, **fit_kwargs), t_summary_table(g, **stat_kwargs) )
    fit = Pipeline('fit_fake', tasks_fit)

    fit(project, obs)

The pipeline is equivalent to:
table = (project, obs) | t_fit(g, **fit_kwargs) | t_summary_table(g, **stat_kwargs)


Todo: adding knobs to turn disk writing on-off

"""

import sys
import numpy as np
import tables
from beast.core import grid
from beast.core.odict import odict
from beast.proba import N_logLikelihood, SN_logLikelihood
from beast.proba import expectation, percentile, getNorm_lnP
from beast.tools.pbar import Pbar
from beast.external.eztables import Table
from beast.external.ezpipe.helpers import RequiredFile, task_decorator


__all__ = ['fit_model_seds_pytables', 't_fit', 'summary_table', 't_summary_table']


from functools import wraps


def generator(func):
    """ Allow to clearly read when a function is a generator in the code
    Do nothing more.
    """
    @wraps(func)
    def wrap(*args, **kwargs):
        return func(*args, **kwargs)
    return wrap


@generator
def tee(args, func_seq):
    """ apply a sequence of functions onto the same input and returns the
    corresponding sequence

    Parameters
    ----------
    args: object
        object corresponding to the common inputs to the functions

    func_seq: sequence
        sequence of callable objects that will take args as input and return values
    """
    for func in func_seq:
        yield func(args)


def compute_sparse_lnp(obk, model_seds, threshold=40, **kwargs):
    """ compute the likelihood of one observation and return necessary
    information filtered according to the threshold value. Any lnp such that we return
    only (models, values) that verify:

    .. math::
        |max(lnp) - lnp| <= threshold

    The likelihood itself will be either N_logLikelihood or SN_logLikelihood
    depending on the obk input length.

    Parameters
    ----------
    obk: sequence (len= 3 or 4)
        sequence representing an observed sed with necessary information
        obk is expected to be either
        3 elements: (sed, err, mask) or
        4 elements: (sed, errp, errm, mask)

    model_seds: ndarray, dtype=float, shape=(nfilters, nmodels)
        Model grid of SEDs to compute the likelihood on

    threshold: float, optional (default 40)
        threshold on the delta lnp values to keep.

    Returns
    -------
    input_arr: ndarray, dtype=float, shape=(nfilters, 3) or (nfilters, 4)
        exact input values seen by the likelihood function into one single array

    idx_arr: ndarray, dtype=int64
        index of SEDs in model_seds that satisfy the sparse selection

    lnp_arr: ndarray, dtype=float
        values of the log-likelihood of selected models

    model_seds: ndarray, dtype=float, shape=(nfilters, nmodels)
        Model grid of SEDs to compute the likelihood on
    """
    if len(obk) == 3:
        (sed, err, mask) = obk
        lnp = N_logLikelihood(  sed, err, model_seds, mask=mask.astype(np.int32), lnp_threshold=abs(threshold) )
    elif len(obk) == 4:
        (sed, errp, errm, mask) = obk
        lnp = SN_logLikelihood(  sed, errp, errm, model_seds, mask=mask.astype(np.int32), lnp_threshold=abs(threshold) )
    else:
        raise AttributeError('getObs is expected to return 3 or 4 values, got {0}'.format(len(obk)))

    indx = np.where((lnp - max(lnp[np.isfinite(lnp)])) > threshold)
    if len(obk) == 4:
        input_arr = np.array([sed, errp, errm, mask]).T
    else:
        input_arr = np.array([sed, err, mask]).T

    idx_arr = np.array(indx[0], dtype=np.int64)
    lnp_arr = np.array(lnp[indx[0]], dtype=np.float32)

    return input_arr, idx_arr, lnp_arr, model_seds


def write_to(outfile, input_arr, idx_arr, lnp_arr, tn, flush=True):
    """ Export results from compute_sparse_lnp onto disk (HDF file)

    for every star, it generates a group named 'star <tn>' in which it saves:
        input: input_arr,
        idx: idx_arr,
        lnp: lnp_arr

    Parameters
    ----------
    outfile: tables.File instance
        opened file with writing rights in which results will be stored

    input_arr: ndarray, dtype=float, shape=(nfilters, 3) or (nfilters, 4)
        exact input values seen by the likelihood function into one single array

    idx_arr: ndarray, dtype=int64
        index of SEDs in model_seds that satisfy the sparse selection

    lnp_arr: ndarray, dtype=float
        values of the log-likelihood of selected models

    tn: int
        star index used as reference to the star (hdf group)

    flush: bool, optional (default=True)
        if set, trigger a flush on the File instance after writing.

    Returns
    -------
    None
    """
    #Need ragged arrays rather than uniform table
    star_group = outfile.createGroup('/', 'star_%d'  % tn, title="star %d" % tn)
    outfile.createArray(star_group, 'input', input_arr)
    outfile.createArray(star_group, 'idx', idx_arr)
    outfile.createArray(star_group, 'lnp', lnp_arr)

    if flush:
        #commit changes
        outfile.flush()


def fit_model_seds_pytables(obs, sedgrid, threshold=-40, outname='lnp.hd5', gridbackend='cache', **kwargs):
    """
    Fit model seds with noise for sensitivity tests

    parameters
    ----------
    obs: Observations
        Observation object to analyze

    sedgrid: SpectralGrid
        stellar model SEDs (luminosities)

    threshold: float
        toss out grid points where lnp - lnp_max < threshold
        This value defined how sparse the final storage will be

    outname: str, optional (default='lnp.hd5')
        path to the file that will store all the computations

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    returns
    -------
        None
    """
    filters = obs.getFilters()

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid, backend=gridbackend)
    else:
        g0 = sedgrid

    with tables.openFile(outname, 'w') as outfile:
        #Save wavelengths in root, remember #n_stars = root._v_nchildren -1
        outfile.createArray(outfile.root, 'grid_waves', g0.lamb[:])
        outfile.createArray(outfile.root, 'obs_filters', filters[:])

        #loop over the obs and do the work
        if hasattr(g0.seds, 'read'):
            _seds = g0.seds.read()
        else:
            _seds = g0.seds

        for tn, obk in Pbar(len(obs), desc='Calculating Lnp').iterover(obs.enumobs()):

            #compute the log-likelihood
            input_arr, idx_arr, lnp_arr, model_seds = compute_sparse_lnp(obk, _seds, threshold)

            # include grid sampling prior
            #lnp = lnp - np.log(g0['Density'] / g0['Density'].sum())

            #Need ragged arrays rather than uniform table
            write_to(outfile, input_arr, idx_arr, lnp_arr, model_seds, flush=True)

    return outname


@generator
def fit_model_seds_memory(obs, sedgrid, threshold=-40, gridbackend='cache', **kwargs):
    """
    Fit model seds with noise for sensitivity tests

    .. note::
        this function is a generator that iterates over each observation

    parameters
    ----------
    obs: Observations
        Observation object to analyze

    sedgrid: SpectralGrid
        stellar model SEDs (luminosities)

    threshold: float
        toss out grid points where lnp - lnp_max < threshold
        This value defined how sparse the final storage will be

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    returns
    -------
        None

    yields
    ------
    input_arr: ndarray, dtype=float, shape=(nfilters, 3) or (nfilters, 4)
        exact input values seen by the likelihood function into one single array

    idx_arr: ndarray, dtype=int64
        index of SEDs in model_seds that satisfy the sparse selection

    lnp_arr: ndarray, dtype=float
        values of the log-likelihood of selected models

    model_seds: ndarray, dtype=float, shape=(nfilters, nmodels)
        Model grid of SEDs to compute the likelihood on
    """
    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid, backend=gridbackend)
    else:
        g0 = sedgrid

    #loop over the obs and do the work
    if hasattr(g0.seds, 'read'):
        _seds = g0.seds.read()
    else:
        _seds = g0.seds

    for tn, obk in Pbar(len(obs), desc='Calculating Lnp').iterover(obs.enumobs()):

        #compute the log-likelihood
        input_arr, idx_arr, lnp_arr, model_seds = compute_sparse_lnp(obk, _seds, threshold)

        # include grid sampling prior
        #lnp = lnp - np.log(g0['Density'] / g0['Density'].sum())

        yield input_arr, idx_arr, lnp_arr, model_seds


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

        keywords
        --------
        node: tables.table.Table
            the Table node to get data from (need named columns)

        expr: str
            the mathematical expression which could include numpy functions or
            simply be a column name

        condvars: dict
            The condvars mapping may be used to define the variable names
            appearing in the condition.

        returns
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


def prepare_for_stats(input_arr, idx_arr, lnp_arr, sedgrid, qname, prior=None, gridbackend='cache'):
    """ Common initialization to computing statistics:
    Prepare input arguments, compute normalizations and apply priors

    Parameters
    ----------
    input_arr: ndarray, dtype=float, shape=(nfilters, 3) or (nfilters, 4)
        exact input values seen by the likelihood function into one single array

    idx_arr: ndarray, dtype=int64
        index of SEDs in model_seds that satisfy the sparse selection

    lnp_arr: ndarray, dtype=float
        values of the log-likelihood of selected models

    model_seds: ndarray, dtype=float, shape=(nfilters, nmodels)
        Model grid of SEDs to compute the likelihood on

    sedgrid: str or grid.SEDgrid instance
        model grid

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    prior: ndarray[float, ndim=1]
        prior probability of each model on the grid to apply (default 1/Nmodel)
        see: compute_uniform_prior

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    Returns
    -------
    q: recarray
        parameters of the models given by idx_arr

    seds: ndarray
        seds of the selected models

    weights: ndarray
        posterior probability values of each model, weights to apply when doing statistics
    """

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid, backend=gridbackend)
    else:
        g0 = sedgrid

    # prepare extraction:
    ## extract models
    q = g0[idx_arr]
    seds = g0.seds[idx_arr, :]

    ## compute log-likelihood norm
    log_norm = np.log(getNorm_lnP(lnp_arr))
    if not np.isfinite(log_norm):
        log_norm = lnp_arr.max()

    ## convert to weights incl. priors
    weights = np.exp(lnp_arr - log_norm)
    if prior is not None:
        weights *= prior[idx_arr]

    return q, seds, weights


def Q_best(input_arr, idx_arr, lnp_arr, sedgrid, qname, prior=None, gridbackend='cache'):
    """ Best Property values:

    Parameters
    ----------
    input_arr: ndarray, dtype=float, shape=(nfilters, 3) or (nfilters, 4)
        exact input values seen by the likelihood function into one single array

    idx_arr: ndarray, dtype=int64
        index of SEDs in model_seds that satisfy the sparse selection

    lnp_arr: ndarray, dtype=float
        values of the log-likelihood of selected models

    model_seds: ndarray, dtype=float, shape=(nfilters, nmodels)
        Model grid of SEDs to compute the likelihood on

    sedgrid: str or grid.SEDgrid instance
        model grid

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    prior: ndarray[float, ndim=1]
        prior probability of each model on the grid to apply (default 1/Nmodel)
        see: compute_uniform_prior

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    returns
    -------
    e_dict: dict or ndarray[float, ndim=1]
        if qname is iterable, returns a dict with a (qname, ndarray) pairs
        else returns only the ndarray of best values (one per obj in cllist)

    sed: ndarray, shape=(nfilters,1)
        best sed, i.e., most probable when given priors or most likely otherwise
    """
    q, seds, weights = prepare_for_stats(input_arr, idx_arr, lnp_arr, sedgrid,
                                         qname, prior=prior,
                                         gridbackend=gridbackend)

    # specific statistic
    best_idx = weights.argmax()
    q = q[best_idx]

    #get quantities
    if hasattr(qname, '__iter__'):
        r = odict()
        for qk in qname:
            lbl = '{0:s}_Best'.format(qk)
            r[lbl] = q[qk]
    else:
        r = q[qname]

    return r, seds[best_idx]


def Q_expect(input_arr, idx_arr, lnp_arr, sedgrid, qname, prior=None, gridbackend='cache'):
    """ Expectation values of any given grid property (incl. expression) but seds,
    i.e.:
            integral(p(q) * q dq) / integral(p(q) dq),
    which in a discrete world becomes
            sum(p(q_i) * q_i) / sum(p(q_i)

    see sed_expect for sed expectation values

    Parameters
    ----------
    input_arr: ndarray, dtype=float, shape=(nfilters, 3) or (nfilters, 4)
        exact input values seen by the likelihood function into one single array

    idx_arr: ndarray, dtype=int64
        index of SEDs in model_seds that satisfy the sparse selection

    lnp_arr: ndarray, dtype=float
        values of the log-likelihood of selected models

    model_seds: ndarray, dtype=float, shape=(nfilters, nmodels)
        Model grid of SEDs to compute the likelihood on

    sedgrid: str or grid.SEDgrid instance
        model grid

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    prior: ndarray[float, ndim=1]
        prior probability of each model on the grid to apply (default 1/Nmodel)
        see: compute_uniform_prior

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    returns
    ------
    e_dict: dict or ndarray[float, ndim=1]
        if qname is iterable, returns a dict with a (qname, ndarray) pair
        else returns only the ndarray of expected values

    sed: ndarray, shape=(nfilters, 1)
        expected SED from the probability distribution
    """
    q, seds, weights = prepare_for_stats(input_arr, idx_arr, lnp_arr, sedgrid,
                                         qname, prior=prior,
                                         gridbackend=gridbackend)

    #get quantities
    if hasattr(qname, '__iter__'):
        r = odict()
        for qk in qname:
            lbl = '{0:s}_E'.format(qk)
            r[lbl] = expectation(q[qk], weights=weights)
    else:
        r = expectation(q[qname], weights=weights)

    #get expected SED
    sed = np.zeros(seds.shape[1], dtype=float)
    for fk in range(seds.shape[1]):
        sed[fk] = expectation(seds[:, fk], weights=weights)

    return r, sed


def Q_percentile(lnpfile, sedgrid, qname, p=[16., 50., 84.], objlist=None, prior=None, gridbackend='cache'):
    """ Percentile values of any given grid property (incl. expression) but seds,

    see also:    sed_percentile, percentile

    keywords
    --------

    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    sedgrid: str or grid.SEDgrid instance
        model grid

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    p: array-like
        list of percentile values

    objlist: list or array like
        index numbers of objects to extract

    prior: ndarray[float, ndim=1]
        prior probability of each model on the grid to apply (default 1/Nmodel)
        see: compute_uniform_prior

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    returns
    -------
    e_dict: dict or ndarray[float, ndim=2]
        if qname is iterable, returns a dict with a (qname, ndarray) pair
        else returns only the ndarray of percentile values (one per obj in cllist)
    """
    if type(lnpfile) == str:
        f = tables.openFile(lnpfile)
    elif isinstance(lnpfile, tables.file.File):
        f = lnpfile
    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid, backend=gridbackend)
    else:
        g0 = sedgrid

    if objlist is None:
        # nchildren - 2 since wavelength, and filters are also saved there
        nobs = f.root._v_nchildren - 2
        objlist = range(nobs)
    else:
        nobs = len(objlist)

    #get quantities
    if hasattr(qname, '__iter__'):
        r = odict()
        for qk in qname:
            tmp = Q_percentile(f, g0, qk, p, objlist, prior)
            for ek, pk in enumerate(p):
                lbl = '{0:s}_p{1:d}'.format(qk, int(pk))
                r[lbl] = tmp[:, ek]
    else:
        #get grid node
        #would be nicer in memory to work on disk here
        #node = f.getNode('/mods_grid')
        #q = get_Q_from_node(node, qname)
        q = g0[qname]
        nval = len(p)
        _p = np.asarray(p, dtype=float)
        r = np.empty((len(objlist), nval), dtype=float)
        with Pbar(nobs, desc='Percentiles') as pb:
            for e, obj in pb.iterover(enumerate(objlist)):
                pb.desc = 'Percentiles({0})'.format(qname)
                lnps = f.getNode('/star_{0:d}/lnp'.format(obj)).read().astype(float)
                indx = f.getNode('/star_{0:d}/idx'.format(obj)).read().astype(int)
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


def summary_table(lnpfname, obs, sedgrid, keys=None, method=None, outname=None, gridbackend='cache'):
    """
    keywords
    --------

    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    sedgrid: str or grid.SEDgrid instance
        model grid

    keys: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    method: str or list of str
        method must be in ['expectation', 'best', 'percentile']

    outname: str
        if set, save the table into this file

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    returns
    -------
    tab: eztable.Table
        table object containing all the statistics
    """
    if type(lnpfname) == str:
        lnpfile = tables.openFile(lnpfname, 'r')
    else:
        lnpfile = lnpfname

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid, backend=gridbackend)
    else:
        g0 = sedgrid

    if keys is None:
        keys = g0.keys()

    if method is None:
        method = 'best expectation percentile'.split()

    #make sure keys are real keys
    for key in keys:
        if not (key in g0.keys()):
            raise KeyError('Key "{0}" not recognized'.format(key))

    r = odict()
    if ('expectation' in method):
        r.update( Q_expect(lnpfile, g0, keys) )

    if ('best' in method):
        r.update( Q_best(lnpfile, g0, keys) )

    if ('percentile' in method):
        r.update(Q_percentile(lnpfile, g0, keys, p=[16., 50., 84.]))

    summary_tab = Table(r, name="Summary Table")

    if not isinstance(lnpfname, tables.file.File):
        lnpfile.close()

    if outname is not None:
        summary_tab.write(outname)

    return summary_tab


#---------------------------------------------------------
# Pipeline interface                        [sec:pipeline]
#---------------------------------------------------------

@task_decorator(logger=sys.stdout)
def t_fit(project, obs, g, threshold=-40, gridbackend='cache', outname=None):
    """t_fit -- run the fitting part

    keywords
    --------

    project: str
        token of the project this task belongs to

    obs: Observation object instance
        observation catalog

    g: grid.SpectralGrid instance
        SED model grid instance

    threshold: float
        toss out grid points where lnp - lnp_max < threshold
        This value defined how sparse the final storage will be

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    outname: str, optional
        optional filename and path to store the outputs

    returns
    -------
    project: str
        token of the project this task belongs to

    lnp_source: str
        file in which sparse lnp values are stored

    obs: Observation object instance
        observation catalog
    """
    if outname is None:
        outname = '{0}_lnp.hd5'.format(project)
    else:
        outname = '{0}_lnp.hd5'.format(outname)
    lnp_source = RequiredFile(outname, fit_model_seds_pytables, obs, g, threshold=threshold, outname=outname, gridbackend=gridbackend)
    return project, lnp_source(), obs


@task_decorator(logger=sys.stdout)
def t_summary_table(project, lnpfname, obs, sedgrid, keys=None, method=None, gridbackend='cache', outname=None):
    """t_summary_table -- task to generate the summary table

    keywords
    --------
    project: str
        token of the project this task belongs to

    lnpfname: str
        file in which sparse lnp values are stored

    obs: Observation object instance
        observation catalog

    sedgrid: grid.SpectralGrid instance
        SED model grid instance

    keys: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    method: str or list of str
        method must be in ['expectation', 'best', 'percentile']

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    outname: str, optional
        optional filename and path to store the outputs

    returns
    -------
    project: str
        token of the project this task belongs to

    stats: eztable.Table
        statistics table

    obs: Observation object instance
        observation catalog

    sedgrid: grid.SpectralGrid instance
        SED model grid instance
    """
    if outname is None:
        outname = '{0}_stats.fits'.format(project)
    else:
        outname = '{0}_stats.fits'.format(outname)
    stat_source = RequiredFile(outname, summary_table, lnpfname, obs, sedgrid, keys=keys, method=method, outname=outname, gridbackend=gridbackend)
    return project, stat_source(), obs, sedgrid
