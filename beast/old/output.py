from beast.external.eztables import Table
import numpy as np
from beast.old import grid


def weighted_percentile(data, wt, percentiles):
    """Compute weighted percentiles.
    If the weights are equal, this is the same as normal percentiles.
    Elements of the C{data} and C{wt} arrays correspond to
    each other and must have equal length (unless C{wt} is C{None}).

    @param data: The data.
    @type data: A L{np.ndarray} array or a C{list} of numbers.
    @param wt: How important is a given piece of data.
    @type wt: C{None} or a L{np.ndarray} array or a C{list} of numbers.
        All the weights must be non-negative and the sum must be
        greater than zero.
    @param percentiles: what percentiles to use.  (Not really percentiles,
        as the range is 0-1 rather than 0-100.)
    @type percentiles: a C{list} of numbers between 0 and 1.
    @rtype: [ C{float}, ... ]
    @return: the weighted percentiles of the data.
    """
    assert np.greater_equal(percentiles, 0.0).all(), "Percentiles less than zero"
    assert np.less_equal(percentiles, 1.0).all(), "Percentiles greater than one"
    data = np.asarray(data)
    assert len(data.shape) == 1
    if wt is None:
        wt = np.ones(data.shape, np.float)
    else:
        wt = np.asarray(wt, np.float)
        assert wt.shape == data.shape
        assert np.greater_equal(wt, 0.0).all(), "Not all weights are non-negative."
    assert len(wt.shape) == 1
    n = data.shape[0]
    assert n > 0
    i = np.argsort(data)
    sd = np.take(data, i, axis=0)
    sw = np.take(wt, i, axis=0)
    aw = np.add.accumulate(sw)
    if not aw[-1] > 0:
        raise ValueError("Nonpositive weight sum")
    w = (aw - 0.5 * sw) / aw[-1]
    spots = np.searchsorted(w, percentiles)
    o = []
    for (s, p) in zip(spots, percentiles):
        if s == 0:
            o.append(sd[0])
        elif s == n:
            o.append(sd[n - 1])
        else:
            f1 = (w[s] - p) / (w[s] - w[s - 1])
            f2 = (p - w[s - 1]) / (w[s] - w[s - 1])
            assert (f1 >= 0) & (f2 >= 0) & (f1 <= 1) & (f2 <= 1)
            assert abs(f1 + f2 - 1.0) < 1e-6
            o.append(sd[s - 1] * f1 + sd[s] * f2)
    return o


def summaryStats(g, indxs, lnp, filters, Q='logA logM Z Av Rv f_bump logT logg logL'):
    """ Generate a dictionary of summary statistics per quantity Q contained in
    the grid and for which the log-likelihood where computed
    INPUTS:
        g       grid                    grid of model
        indxs   list                    model indices (for sparse computations)
        lnp     ndarray[float, ndim=1]  log-likelihood values on the grid
        filters list                    list of filter names to use for sed predicting
    KEYWORDS:
        Q       str                     list of quantities separated by spaces
    OUTPUTS:
        d       dict                    dictionary of statistics
    """
    _q = Q.split()
    _r = np.exp(lnp)
    _r /= _r.sum()

    d = {}

    for k in _q:
        try:
            m0, m1, m2 = weighted_percentile( g.grid[k][indxs], _r, [0.16, 0.5, 0.84] )
        except:
            m0, m1, m2 = np.nan, np.nan, np.nan

        d['%s_p16' % k ] = m0
        d['%s_p50' % k ] = m1
        d['%s_p84' % k ] = m2

    for k in range(len(filters)):
        try:
            m0, m1, m2 = weighted_percentile( g.seds[:, k][indxs], _r, [0.16, 0.5, 0.84] )
        except:
            m0, m1, m2 = np.nan, np.nan, np.nan
        d['%s_p16' % (filters[k]) ] = m0
        d['%s_p50' % (filters[k]) ] = m1
        d['%s_p84' % (filters[k]) ] = m2

    return d


def summaryTable(files, grid_filename, filters, outfilename, Q='logA logM Z Av Rv f_bump logT logg logL', **kwargs):

    # get the nD stellar/dust SED grid
    ext_grid = grid.FileSEDGrid(grid_filename)

    d = {}
    for sfile in files:

        print 'working on ', sfile

        # get the nD likelihood function for each star
        lnp = Table(sfile, extension='LNP')

        # get the summary for each parameter
        indxs = lnp['idx'].astype(int)
        print len(indxs)

        _d = summaryStats(ext_grid, indxs, lnp['lnp'], filters, Q=Q)
        if len(d) == 0:
            for k, v in _d.items():
                d[k] = np.array([v])
        else:
            for k, v in _d.items():
                d[k] = np.hstack([d[k], np.array([v])])

    t = Table(d, name='SED analysis summury table')

    # add header information if passed
    for k in kwargs:
        t.header[k] = kwargs[k]

    t.write(outfilename)


# probably not needed or needing modifcation
def fullTable(fname, g, r, Av, lamb, fakesed, fakeerr, filters, **kwargs):
    import mypickleshare
    d = mypickleshare.PickleShareDB(fname)
    d['lnp'] = r
    d['Av'] = Av
    d['filters'] = filters
    d['GRID'] = g.grid.header['SOURCE']
    d['LAMB'] = lamb
    d['INPUTSED'] = fakesed
    d['INPUTERR'] = fakeerr

    for k in kwargs:
        d[k] = kwargs[k]

    del d

#Tchernyshyov contribution starts here:

# MF: Because it was spectific to the purpose, I moved this to the sensitivity
#     test script
#
#
#import progressbar
#import numexpr
#import tables
#import re
#import glob
#
#
#def get_nstars_mtests(test_dir):
#    #Figure out how many fake stars and realizations there are in a directory
#    fakes = glob.glob(test_dir + 'fake_star*')
#    max_n = 0
#    max_m = 0
#    for fake in fakes:
#        n, m = re.findall("fake_star_\d+_\d+", fake)[0].split('_')[2:]
#        if int(n) >= max_n:
#            max_n = int(n)
#        if int(m) >= max_m:
#            max_m = int(m)
#    return max_n + 1, max_m + 1
#
#
#def expectation_values_eztables(grid, keys=['Av', 'Rv', 'f_bump'], dir='Tests/fake_many_0/'):
#    """
#    args:
#         grid        Model parameters. (model_grid.grid)
#    keywords:
#         keys        Parameters to calculate expectation values for.
#         dir         Directory log likelihoods have been saved to
#
#    TODO: not sure this is useful as is. We need to define the output formats end use that directly.
#    I think not having so many files is better
#    TODO: if we keep this, we need to kill the m_test part (MC is not part a real data)
#    """
#    #make sure keys are real keys
#    for key in keys:
#        if not (key in grid.keys()):
#            raise KeyError('Key "%s" not recognized' % key)
#
#    n_stars, m_tests = get_nstars_mtests(dir)
#
#    means = np.zeros((n_stars, len(keys)))
#    stds = np.zeros((n_stars, len(keys)))
#    fakeinds = np.zeros(n_stars, dtype=np.int16)
#    temp_vals = np.zeros((m_tests, len(keys)))
#    with progressbar.PBar(n_stars, txt="E(key)") as pbar:
#        for i in range(n_stars):
#            for j in range(m_tests):
#                _fname = dir + 'fake_star_' + str(i) + '_' + str(j) + '.fits'
#                f = Table(_fname, silent=True)
#                prob = f.evalexpr('exp(lnp)')
#                indx = f['idx'].astype(int)
#                for e, key in enumerate(keys):
#                    temp_vals[j, e] = numexpr.evaluate('sum(prob * vals)', local_dict={'vals': grid.data[key][indx], 'prob': prob})
#                temp_vals[j, :] /= prob.sum()
#                if j == 0:
#                    fakeinds[i] = f.header.FAKE_IDX
#
#            means[i] = temp_vals.mean(axis=0)
#            stds[i] = temp_vals.std(axis=0)
#            pbar.update(i)
#
#    summary_tab = Table(name="Summary Table")
#    for e, key in enumerate(keys):
#        summary_tab.add_column(key + '_recovered', means[:, e])
#        summary_tab.add_column(key + '_std', stds[:, e])
#        summary_tab.add_column(key, grid.data[key][fakeinds])
#    return summary_tab
#
#
#def expectation_values_pytables(grid, keys=['Av', 'Rv', 'f_bump'], filepath='Tests/fake_many_0/fakestars_0.hf5', method='expectation'):
#    """
#    args:
#         grid        Model parameters. (model_grid.grid)
#    **kwargs:
#         keys        Parameters to calculate expectation values for.
#         filepath    File log likelihoods have been saved to.
#    """
#    #make sure keys are real keys
#    for key in keys:
#        if not (key in grid.keys()):
#            raise KeyError('Key "%s" not recognized' % key)
#    outfile = tables.openFile(filepath, 'r')
#    n_stars = outfile.root._v_nchildren - 1  # nchildren - 1 since wavelength also saved there
#    m_tests = outfile.root.fakeStar_0._v_nchildren
#
#    # Mean of m_tests realizations for each star
#    means = np.zeros((n_stars, len(keys)), dtype=np.float32)
#    #Standard deviation of m_tests realizations for each star
#    stds = np.zeros((n_stars, len(keys)), dtype=np.float32)
#    fakeinds = np.zeros(n_stars, dtype=np.int32)
#    temp_vals = np.zeros((m_tests, len(keys)), dtype=np.float32)
#
#    with progressbar.PBar(n_stars, txt="E(key)") as pbar:
#        for i in range(n_stars):
#            for j in range(m_tests):
#                if (method == 'expectation'):
#                    #the usage of numexpr here is not necessary. (Not to mention tables already has a better version built in
#                    prob = numexpr.evaluate('exp(lnp)', local_dict={'lnp': outfile.getNode('/fakeStar_%d/fake_%d/lnp' % (i, j)).read()})
#                    indx = outfile.getNode('/fakeStar_%d/fake_%d/idx' % (i, j)).read()
#                    for e, key in enumerate(keys):
#                        temp_vals[j, e] = numexpr.evaluate('sum(prob*vals)', local_dict={'vals': grid.getCol(key)[indx], 'prob': prob})
#                    temp_vals[j, :] /= prob.sum()
#                elif (method == 'maxprob'):
#                    lnp = outfile.getNode('/fakeStar_%d/fake_%d/lnp' % (i, j)).read()
#                    indx = outfile.getNode('/fakeStar_%d/fake_%d/idx' % (i, j)).read()
#                    sel = grid.getRow(indx[lnp.argmax()])
#                    for e, key in enumerate(keys):
#                        temp_vals[j, e] = sel[key]
#
#            fakeinds[i] = outfile.getNode('/fakeStar_%d' % i)._v_attrs.fakein
#            means[i] = temp_vals.mean(axis=0)
#            stds[i] = temp_vals.std(axis=0)
#            pbar.update(i)
#    outfile.close()
#    summary_tab = Table(name="Summary Table")
#    for e, key in enumerate(keys):
#        summary_tab.add_column(key + '_recovered', means[:, e])
#        summary_tab.add_column(key + '_std', stds[:, e])
#        summary_tab.add_column(key, grid.data[key][fakeinds])
#    return summary_tab
