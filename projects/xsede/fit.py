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
"""

import sys
import numpy as np
import tables
from beast.core import grid
from beast.core.odict import odict
from beast.proba.likelihood import *
from beast.proba import expectation, percentile, getNorm_lnP
from beast.tools.pbar import Pbar
from beast.external.eztables import Table
from beast.external.ezpipe.helpers import RequiredFile, task_decorator


__all__ = ['fit_model_seds_pytables', 't_fit', 'summary_table', 't_summary_table']


def fit_model_seds_pytables(obs, sedgrid, ast, threshold=-40, outname='lnp.hd5', gridbackend='cache'):
    """
    Fit model seds with noise for sensitivity tests

    keywords
    --------
    obs: Observations
        Observation object to analyze

    sedgrid: SpectralGrid
        stellar model SEDs (luminosities)

    threshold: float
        toss out grid points where lnp - lnp_max < threshold
        This value defined how sparse the final storage will be

    outname: string
        output file directory for results

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

    ast_error = ast.root.error[:]
    ast_bias = ast.root.bias[:]

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
            #if len(obk) == 3:
            #    (sed, err, mask) = obk
            #    lnp = N_logLikelihood(  sed, err, _seds, mask=mask.astype(np.int32), lnp_threshold=abs(threshold) )
            #elif len(obk) == 4:
            #    (sed, errp, errm, mask) = obk
            #    lnp = SN_logLikelihood(  sed, errp, errm, _seds, mask=mask.astype(np.int32), lnp_threshold=abs(threshold) )
            #else:
            #    raise AttributeError('getObs is expected to return 3 or 4 values, got {0}'.format(len(obk)))
            (sed) = obk
            (lnp,chi2) = N_logLikelihood_NM(sed,_seds,ast_error,ast_bias,mask=None, lnp_threshold=abs(threshold) )
            # include grid sampling prior
            #lnp = lnp - np.log(g0['Density'] / g0['Density'].sum())
            #print len(lnp)
            #Need ragged arrays rather than uniform table
            star_group = outfile.createGroup('/', 'star_%d'  % tn, title="star %d" % tn)
            indx = np.where((lnp - max(lnp[np.isfinite(lnp)])) > threshold)
            #if len(obk) == 4:
            #    outfile.createArray(star_group, 'input', np.array([sed, errp, errm, mask]).T)
            #else:
            #    outfile.createArray(star_group, 'input', np.array([sed, err, mask]).T)
            outfile.createArray(star_group, 'input', np.array([sed]).T)
            outfile.createArray(star_group, 'idx', np.array(indx[0], dtype=np.int64))
            outfile.createArray(star_group, 'lnp', np.array(lnp[indx[0]], dtype=np.float32))
            outfile.createArray(star_group, 'chi2', np.array(chi2[indx[0]], dtype=np.float32))
            #commit changes
            outfile.flush()

    return outname


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


def Q_expect(lnpfile, sedgrid, qname, objlist=None, prior=None, gridbackend='cache'):
    """ Expectation values of any given grid property (incl. expression) but seds,
    i.e.:
            integral(p(q) * q dq) / integral(p(q) dq),
    which in a discrete world becomes
            sum(p(q_i) * q_i) / sum(p(q_i)

    see sed_expect for sed expectation values

    keywords
    --------

    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    sedgrid: str or grid.SEDgrid instance
        model grid

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    objlist: list or array like
        index numbers of objects to extract

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
        else returns only the ndarray of expected values (one per obj in objlist)
    """
    if type(lnpfile) == str:
        f = tables.openFile(lnpfile)
    else:
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
            lbl = '{0:s}_E'.format(qk)
            r[lbl] = Q_expect(f, g0, qk, objlist, prior)
    else:
        #make sure keys are useful keys
        if not (qname in g0.keys()):
            raise KeyError('Key "{0}" not recognized'.format(qname))
        #get grid node
        #would be nicer in memory to work on disk here
        #node = f.getNode('/mods_grid')
        #q = get_Q_from_node(node, qname)
        q = g0[qname]
        r = np.empty(len(objlist), dtype=float)

        with Pbar(nobs, txt='Expectations') as pb:
            for e, obj in pb.iterover(enumerate(objlist)):
                pb.desc = 'E({0})'.format(qname)
                lnps = f.getNode('/star_{0:d}/lnp'.format(obj)).read().astype(float)
                indx = f.getNode('/star_{0:d}/idx'.format(obj)).read().astype(int)
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


def Q_best(lnpfile, sedgrid, qname, objlist=None, prior=None, gridbackend='cache'):
    """ Best Property values:

    Note: external loop on Q whereas cllist to extract data only once!

    Note on parallel version:
        # from ised.external.ezmap import map
        # def job(qk):
        #    return (qk, Q_best(lnpfile, qk, cllist, prior))
        # r.update( map(job, qname, ncpu=len(qname)) )

    keywords
    --------

    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    sedgrid: str or grid.SEDgrid instance
        model grid

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

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
    e_dict: dict or ndarray[float, ndim=1]
        if qname is iterable, returns a dict with a (qname, ndarray) pairs
        else returns only the ndarray of best values (one per obj in cllist)
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
            lbl = '{0:s}_Best'.format(qk)
            r[lbl] = Q_best(f, g0, qk, objlist, prior)
    else:
        #get grid node
        #would be nicer in memory to work on disk here
        #node = f.getNode('/mods_grid')
        #q = get_Q_from_node(node, qname)
        q = g0[qname]
        r = np.empty(len(objlist), dtype=float)
        with Pbar(nobs, txt='Best') as pb:
            for e, obj in pb.iterover(enumerate(objlist)):
                pb.desc = 'Best({0})'.format(qname)
                lnps = f.getNode('/star_{0:d}/lnp'.format(obj)).read().astype(float)
                indx = f.getNode('/star_{0:d}/idx'.format(obj)).read().astype(int)
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
    keys.remove('osl')
    keys.remove('keep')
    keys.remove('weight')
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
def t_fit(project, obs, g, ast, threshold=-40, gridbackend='cache', outname=None):
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
    lnp_source = RequiredFile(outname, fit_model_seds_pytables, obs, g, ast, threshold=threshold, outname=outname, gridbackend=gridbackend)
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
