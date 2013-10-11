import numpy as np
import tables
from beast.core import grid
from beast.tools import progressbar
from beast.proba import N_logLikelihood
from beast.proba import expectation, percentile, getNorm_lnP
from beast.external.eztables import Table
from beast.external.eztables.core.odict import odict
from ezpipe.helpers import RequiredFile, task_decorator
import sys


__all__ = ['fit_model_seds_pytables', 't_fit', 'summary_table', 't_summary_table']


def fit_model_seds_pytables(obs, sedgrid, threshold=-40, outname='lnp.hd5'):

    """
    Fit model seds with noise for sensitivity tests
    INPUTS:
        obs         Observtions     Observation object to analyze
        sedgrid     SEDgrid         stellar model SEDs (luminosities)

    KEYWORDS:
        threshold   float          toss out grid points where lnp - lnp_max < threshold
        outname     string         output file directory for results

    TODO: Clean up the disk output structures
    HDF is not a bad choice since it compress and keep all in one file!
    """
    filters = obs.getFilters()

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid)
    else:
        g0 = sedgrid

    with tables.openFile(outname, 'w') as outfile:
        #Save wavelengths in root, remember #n_stars = root._v_nchildren -1
        outfile.createArray(outfile.root, 'grid_waves', g0.lamb)
        outfile.createArray(outfile.root, 'obs_filters', filters)

        #loop over the obs and do the work
        with progressbar.PBar(len(obs), txt="Calculating lnp") as pbar:

            for tn, (sed, err, mask) in obs.enumobs():

                lnp = N_logLikelihood(  sed, err, g0.seds, mask=mask, lnp_threshold=abs(threshold) )

                #Need ragged arrays rather than uniform table
                star_group = outfile.createGroup('/', 'star_%d'  % tn, title="star %d" % tn)
                indx = np.where((lnp - max(lnp[np.isfinite(lnp)])) > -40.)
                outfile.createArray(star_group, 'input', np.array([sed, err, mask]).T)
                outfile.createArray(star_group, 'idx', np.array(indx[0], dtype=np.int32))
                outfile.createArray(star_group, 'lnp', np.array(lnp[indx[0]], dtype=np.float32))
                #commit changes
                outfile.flush()

                pbar.update(tn, force=True)  # Forcing because it can be long to show the first ETA


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


def Q_expect(lnpfile, sedgrid, qname, objlist=None, prior=None):
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

    sedgrid: str or grid.SEDgrid instance
        model grid

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    KEYWORDS
    --------

    objlist: list or array like
        index numbers of objects to extract

    prior: ndarray[float, ndim=1]
        prior probability of each model on the grid to apply (default 1/Nmodel)
        see: compute_uniform_prior

    OUTPUT
    ------
    e_dict: dict or ndarray[float, ndim=1]
        if qname is iterable, returns a dict with a (qname, ndarray) pair
        else returns only the ndarray of expected values (one per obj in objlist)
    """
    if type(lnpfile) == str:
        f = tables.openFile(lnpfile)
    elif isinstance(lnpfile, tables.file.File):
        f = lnpfile

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid)
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
            lbl = '{:s}_E'.format(qk)
            r[lbl] = Q_expect(f, g0, qk, objlist, prior)
    else:
        #make sure keys are real keys
        if not (qname in g0.keys()):
            raise KeyError('Key "{}" not recognized'.format(qname))
        #get grid node
        #would be nicer in memory to work on disk here
        #node = f.getNode('/mods_grid')
        #q = get_Q_from_node(node, qname)
        q = g0.grid[qname]
        r = np.empty(len(objlist), dtype=float)
        with progressbar.PBar(nobs, txt='Expectations') as pb:
            for e, obj in enumerate(objlist):
                pb.update(e, txt='E({})'.format(qname))
                lnps = f.getNode('/star_{:d}/lnp'.format(obj)).read().astype(float)
                indx = f.getNode('/star_{:d}/idx'.format(obj)).read().astype(int)
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


def Q_best(lnpfile, sedgrid, qname, objlist=None, prior=None):
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

    sedgrid: str or grid.SEDgrid instance
        model grid

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    KEYWORDS
    --------

    objlist: list or array like
        index numbers of objects to extract

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

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid)
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
            lbl = '{:s}_Best'.format(qk)
            r[lbl] = Q_best(f, g0, qk, objlist, prior)
    else:
        #get grid node
        #would be nicer in memory to work on disk here
        #node = f.getNode('/mods_grid')
        #q = get_Q_from_node(node, qname)
        q = g0.grid[qname]
        r = np.empty(len(objlist), dtype=float)
        with progressbar.PBar(nobs, txt='Best') as pb:
            for e, obj in enumerate(objlist):
                pb.update(e, txt='Best({})'.format(qname))
                lnps = f.getNode('/star_{:d}/lnp'.format(obj)).read().astype(float)
                indx = f.getNode('/star_{:d}/idx'.format(obj)).read().astype(int)
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


def Q_percentile(lnpfile, sedgrid, qname, p=[16., 50., 84.], objlist=None, prior=None):
    """ Percentile values of any given grid property (incl. expression) but seds,

    see also:    sed_percentile, percentile

    INPUTS
    ------

    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    sedgrid: str or grid.SEDgrid instance
        model grid

    qname: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    KEYWORDS
    --------
    p: array-like
        list of percentile values

    objlist: list or array like
        index numbers of objects to extract

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
    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid)
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
                lbl = '{:s}_p{:d}'.format(qk, int(pk))
                r[lbl] = tmp[:, ek]
    else:
        #get grid node
        #would be nicer in memory to work on disk here
        #node = f.getNode('/mods_grid')
        #q = get_Q_from_node(node, qname)
        q = g0.grid[qname]
        nval = len(p)
        _p = np.asarray(p, dtype=float)
        r = np.empty((len(objlist), nval), dtype=float)
        with progressbar.PBar(nobs, txt='Percentiles') as pb:
            for e, obj in enumerate(objlist):
                pb.update(e, txt='Percentiles({})'.format(qname))
                lnps = f.getNode('/star_{:d}/lnp'.format(obj)).read().astype(float)
                indx = f.getNode('/star_{:d}/idx'.format(obj)).read().astype(int)
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


def summary_table(lnpfname, obs, sedgrid, keys=None, method=None, outname=None):
    """
    INPUTS
    ------
    lnpfile: str or tables.file.File
        lnp file to use given by its path or the open file

    sedgrid: str or grid.SEDgrid instance
        model grid

    keys: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    KEYWORDS
    --------
    keys: str or list of str
        if str:  name of the quantity or expression to evaluate from the grid table
        if list: list of qquantities or expresions

    method: str or list of str
        method must be in ['expectation', 'best', 'percentile']

    outname: str
        if set, save the table into this file
    """
    if type(lnpfname) == str:
        lnpfile = tables.openFile(lnpfname, 'r')
    else:
        lnpfile = lnpfname

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid)
    else:
        g0 = sedgrid

    if keys is None:
        keys = g0.grid.keys()

    if method is None:
        method = 'best expectation percentile'.split()

    #make sure keys are real keys
    for key in keys:
        if not (key in g0.keys()):
            raise KeyError('Key "{}" not recognized'.format(key))

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
def t_fit(project, obs, g, threshold=-40):
    outname = '{}_lnp.hd5'.format(project)
    lnp_source = RequiredFile(outname, fit_model_seds_pytables, obs, g, threshold=threshold, outname=outname)
    return project, lnp_source(), obs


@task_decorator(logger=sys.stdout)
def t_summary_table(project, lnpfname, obs, sedgrid, keys=None, method=None):
    outname = '{}_stats.fits'.format(project)
    stat_source = RequiredFile(outname, summary_table, lnpfname, obs, sedgrid, keys=keys, method=method, outname=outname)
    return project, stat_source(), obs, sedgrid
