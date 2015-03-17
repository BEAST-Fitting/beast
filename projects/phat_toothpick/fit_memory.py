"""
Fitting Pipeline for PHAT
BEAST Toothpick version (v1, 5 Feb 2015)
based on code by Morgan Fouesneau
major modifications by Karl Gordon
  - added a number of additional parameters
  - uses a fast 1D PDF generator
  - combines best,expectation,percentiles for speed
=============

This code uses `ezpipe`, a pipeline package written by Morgan Fouesneau to provide a clean
and flexible interface. In particular this package simplifies the management of
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

import os

import sys
import time
import numpy as np
import tables
import string

import numexpr

from astropy.coordinates import ICRS as ap_ICRS
from astropy import units as ap_units
from astropy.io import fits

from beast.core import grid
from beast.core.odict import odict
from beast.proba.likelihood import *
from beast.proba import expectation, percentile, getNorm_lnP
from beast.tools.pbar import Pbar
#from beast.external.eztables import Table
from astropy.table import Table
from beast.external.ezpipe.helpers import RequiredFile, task_decorator

from beast.core.pdf1d import pdf1d

__all__ = ['summary_table_memory']


def Q_all_memory(obs, sedgrid, ast, qnames, p=[16., 50., 84.], gridbackend='cache', max_nbins=50,
          pdf1d_outname=None, threshold=-40):
    """ Get the best, expectation, and percentile values of all the given grid property
      (done in once function for speed)

    keywords
    --------

    sedgrid: str or grid.SEDgrid instance
        model grid

    qnames: list of quantities or expresions

    p: array-like
        list of percentile values

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    max_nbins: maxiumum number of bins to use for the 1D likelihood calculations

    pdf1d_outname: set to output the 1D PDFs into a FITS file with extensions

    returns
    -------
    e_dict: dict with a (qname, ndarray) pair
    """

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid, backend=gridbackend)
    else:
        g0 = sedgrid

    g0_indxs, = np.where(g0['weight'] > 0.0)
    g0_weights = np.log(g0['weight'][g0_indxs])
    g0_weights_sum = np.log(g0['weight'][g0_indxs].sum())
    g0_weights = numexpr.evaluate("g0_weights - g0_weights_sum")

    ast_error = ast.root.error[:]
    ast_bias = ast.root.bias[:]

    nobs = len(obs)

    # setup the arrays to temp sore the results
    n_qnames = len(qnames)
    n_pers = len(p)
    best_vals = np.zeros((nobs, n_qnames))
    exp_vals = np.zeros((nobs, n_qnames))
    per_vals = np.zeros((nobs, n_qnames, n_pers))
    chi2_vals = np.zeros(nobs)
    chi2_indx = np.zeros(nobs)
    lnp_vals = np.zeros(nobs)
    lnp_indx = np.zeros(nobs)
    best_specgrid_indx = np.zeros(nobs)
    
    # setup the mapping for the 1D PDFs
    fast_pdf1d_objs = []
    save_pdf1d_vals = []
    for qname in qnames:
        q = g0[qname]
        
        n_uniq = len(np.unique(q))
        if len(np.unique(q)) > max_nbins: 
            nbins = max_nbins  # limit the number of bins in the 1D likelihood for speed
        else:
            nbins = n_uniq

        # setup the fast 1d pdf
        ignorebelow = None  # need to know so 'zeros' (defined at -100) are ignored
        if (string.find(qname,'_wd') > 0) | (string.find(qname,'_wd') > 0):
            ignorebelow = -99.99
        _tpdf1d = pdf1d(q, nbins, ignorebelow=ignorebelow)
        fast_pdf1d_objs.append(_tpdf1d)
        
        # setup the arrays to save the 1d PDFs
        save_pdf1d_vals.append(np.zeros((nobs+1, nbins)))
        save_pdf1d_vals[-1][nobs,:] = _tpdf1d.bin_vals

    # loop over the objects and get all the requested quantities
    g0_specgrid_indx = g0['specgrid_indx']
    _p = np.asarray(p, dtype=float)

    #loop over the obs and do the work
    if hasattr(g0.seds, 'read'):
        _seds = g0.seds.read()
    else:
        _seds = g0.seds

    for e, obj in Pbar(len(obs), desc='Calculating Lnp/Stats').iterover(obs.enumobs()):
    #with Pbar(nobs, desc='Best/Exp/Per') as pb:
    #    for e, obj in pb.iterover(enumerate(obs)):
            # get the full nD posterior
            (sed) = obj
            (lnp,chi2) = N_logLikelihood_NM(sed,_seds,ast_error,ast_bias,mask=None, lnp_threshold=abs(threshold) )
            
            lnp = lnp[g0_indxs]
            chi2 = chi2[g0_indxs]
            lnp = numexpr.evaluate('lnp + g0_weights')
            #lnp +=  g0_weights  # multiply by the prior weights (sum in log space)

            indx, = np.where((lnp - max(lnp[np.isfinite(lnp)])) > threshold)

            # now generate the sparse likelihood (remove later if this works by updating code below)
            lnps = lnp[indx]
            chi2 = chi2[indx]

            #log_norm = np.log(getNorm_lnP(lnps))
            #if not np.isfinite(log_norm):
            #    log_norm = lnps.max()
            log_norm = lnps.max()
            weights = np.exp(lnps - log_norm)
                
            # index to the full model grid for the best fit values
            best_full_indx = indx[weights.argmax()]

            # index to the spectral grid 
            best_specgrid_indx[e] = g0_specgrid_indx[best_full_indx]
            
            # goodness of fit quantities
            chi2_vals[e] = chi2.min()
            chi2_indx[e] = indx[chi2.argmin()]
            lnp_vals[e] = lnps.max()
            lnp_indx[e] = best_full_indx

            for k, qname in enumerate(qnames):
                q = g0[qname]

                # best value
                best_vals[e,k] = q[best_full_indx]

                # expectration value
                exp_vals[e,k] = expectation(q[indx], weights=weights)

                # percentile values
                pdf1d_bins, pdf1d_vals = fast_pdf1d_objs[k].gen1d(indx, np.exp(lnps))
                save_pdf1d_vals[k][e,:] = pdf1d_vals
                if pdf1d_vals.max() > 0:
                    pdf1d_vals /= pdf1d_vals.max()
                    per_vals[e,k,:] = percentile(pdf1d_bins, _p, weights=pdf1d_vals)
                else:
                    per_vals[e,k,:] = [0.0,0.0,0.0]

    # populate the dict array
    r = odict()
    for k, qname in enumerate(qnames):
        r['{0:s}_Best'.format(qname)] = best_vals[:,k]
        r['{0:s}_Exp'.format(qname)] = exp_vals[:,k]
        for i, pval in enumerate(p):
            r['{0:s}_p{1:d}'.format(qname, int(pval))] = per_vals[:,k,i]

    r['chi2min'] = chi2_vals
    r['chi2min_indx'] = chi2_indx.astype(int)
    r['Pmax'] = lnp_vals
    r['Pmax_indx'] = lnp_indx.astype(int)
    r['specgrid_indx'] = best_specgrid_indx

    # save the 1D PDFs
    if pdf1d_outname is not None:
        if os.path.isfile(pdf1d_outname):
            os.remove(pdf1d_outname)

        # write a small primary header
        fits.append(pdf1d_outname, np.zeros((2,2)))

        # write the 1D PDFs for all the objects, 1 set per extension
        for k, qname in enumerate(qnames):
            hdu = fits.PrimaryHDU(save_pdf1d_vals[k])
            pheader = hdu.header
            pheader.set('XTENSION','IMAGE') 
            pheader.set('EXTNAME',qname) 
            fits.append(pdf1d_outname, save_pdf1d_vals[k], header=pheader)

    return r

def IAU_names_and_extra_info(obsdata):
    """
    generates IAU approved names for the PHAT data using RA & DEC
      and extra information about the sources (ra, dec, photometry, etc.)

    keywords
    --------
    obs: Observations

    returns
    -------
    e_dict: dict 
        returns a dict with a (name, ndarray) pair
    """
    r = odict()

    # generate the IAU names
    _tnames = []
    for i in range(len(obsdata)):
        c = ap_ICRS(ra=obsdata.data['ra'][i], dec=obsdata.data['dec'][i],
                    unit=(ap_units.degree, ap_units.degree))
        _tnames.append('PHAT J' + 
                       c.ra.to_string(sep="",precision=2,alwayssign=False,pad=True) + 
                       c.dec.to_string(sep="",precision=2,alwayssign=True,pad=True))

    r['Name'] = _tnames

    # other useful information
    r['RA'] = obsdata.data['ra']
    r['DEC'] = obsdata.data['dec']
    r['field'] = obsdata.data['field']
    r['inside_brick'] = obsdata.data['inside_brick']
    r['inside_chipgap'] = obsdata.data['inside_chipgap']

    # include the observed filter fluxes
    for k, filtername in enumerate(obsdata.filters): 
        r[filtername] = (obsdata.data[filtername]*obsdata.vega_flux[k]).astype(float) 

    return r

def summary_table_memory(obs, noisemodel, sedgrid, keys=None, method=None, outname=None, gridbackend='cache'):
    """
    keywords
    --------

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

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid, backend=gridbackend)
    else:
        g0 = sedgrid

    if keys is None:
        keys = g0.keys()

    #make sure keys are real keys
    skip_keys = 'osl keep weight fullgrid_idx stage specgrid_indx'.split()
    keys = [k for k in keys if k not in skip_keys]

    if method is None:
        method = 'all'.split()

    for key in keys:
        if not (key in g0.keys()):
            raise KeyError('Key "{0}" not recognized'.format(key))

    r = {}
    #r = odict()

    # generate an IAU complient name for each source and add other inform
    #r.update( IAU_names_and_extra_info(obs) )
    res = IAU_names_and_extra_info(obs)
    for key in res.keys():
        r[key] = res[key]


    if ('all' in method):
        #r.update(Q_all_memory(obs, g0, noisemodel, keys, p=[16., 50., 84.],
        #                     pdf1d_outname=string.replace(outname,'stats.fits','pdf1d.fits')))
        res = Q_all_memory(obs, g0, noisemodel, keys, p=[16., 50., 84.],threshold=-10.,
                           pdf1d_outname=string.replace(outname,'stats.fits','pdf1d.fits'))
        for key in res.keys():
            r[key] = res[key]

    #summary_tab = Table(r, name="Summary Table")
    summary_tab = Table(r)

    if outname is not None:
        summary_tab.write(outname)

    return summary_tab


#---------------------------------------------------------
# Pipeline interface                        [sec:pipeline]
#---------------------------------------------------------

# needs updating or deleting due to change of doing all in memory (no lnp file)

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
