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

from itertools import islice

import numexpr

from astropy.coordinates import ICRS as ap_ICRS
from astropy import units as ap_units
from astropy.io import fits

from beast.core import grid
from beast.proba.likelihood import *
from beast.proba import expectation, percentile, getNorm_lnP
from beast.tools.pbar import Pbar
from astropy.table import Table

from beast.core.pdf1d import pdf1d

def save_stats(stats_outname, stats_dict_in, best_vals, exp_vals, per_vals, chi2_vals, chi2_indx, 
               lnp_vals, lnp_indx, best_specgrid_indx, qnames, p):

    stats_dict = stats_dict_in.copy()

    # populate the dict array
    for k, qname in enumerate(qnames):
        stats_dict['{0:s}_Best'.format(qname)] = best_vals[:,k]
        stats_dict['{0:s}_Exp'.format(qname)] = exp_vals[:,k]
        for i, pval in enumerate(p):
            stats_dict['{0:s}_p{1:d}'.format(qname, int(pval))] = per_vals[:,k,i]

    stats_dict['chi2min'] = chi2_vals
    stats_dict['chi2min_indx'] = chi2_indx.astype(int)
    stats_dict['Pmax'] = lnp_vals
    stats_dict['Pmax_indx'] = lnp_indx.astype(int)
    stats_dict['specgrid_indx'] = best_specgrid_indx

    summary_tab = Table(stats_dict)

    if stats_outname is not None:
        summary_tab.write(stats_outname, overwrite=True)

def save_pdf1d(pdf1d_outname, save_pdf1d_vals, qnames):

    # write a small primary header
    fits.writeto(pdf1d_outname, np.zeros((2,2)), clobber=True)

    # write the 1D PDFs for all the objects, 1 set per extension
    for k, qname in enumerate(qnames):
        hdu = fits.PrimaryHDU(save_pdf1d_vals[k])
        pheader = hdu.header
        pheader.set('XTENSION','IMAGE') 
        pheader.set('EXTNAME',qname) 
        fits.append(pdf1d_outname, save_pdf1d_vals[k], header=pheader)

def Q_all_memory(prev_result, obs, sedgrid, ast, qnames, p=[16., 50., 84.], gridbackend='cache', max_nbins=50,
                 stats_outname=None, pdf1d_outname=None, lnp_outname=None, lnp_npts=None, save_every_npts=None,
                 threshold=-40, resume=False):
    """ Get the best, expectation, and percentile values of all the given grid property
      (done in once function for speed)

    keywords
    --------
    prev_result: dict
        previous results to include in the output summary table
        usually basic data on each source

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

    # if this is a resume job, read in the already computed stats and fill the variables
    # also - find the start position for the resumed run
    if resume:
        stats_table = Table.read(stats_outname)
        
        for k, qname in enumerate(qnames):
            best_vals[:,k] = stats_table['{0:s}_Best'.format(qname)]
            exp_vals[:,k] = stats_table['{0:s}_Exp'.format(qname)]
            for i, pval in enumerate(p):
                per_vals[:,k,i] = stats_table['{0:s}_p{1:d}'.format(qname, int(pval))]
                
        chi2_vals = stats_table['chi2min']
        chi2_indx = stats_table['chi2min_indx']
        lnp_vals = stats_table['Pmax']
        lnp_indx = stats_table['Pmax_indx']
        best_specgrid_indx = stats_table['specgrid_indx']

        indxs, = np.where(stats_table['Pmax'] != 0.0)
        start_pos = max(indxs) + 1
        print('resuming run with start indx = ' + str(start_pos) + ' out of ' + str(len(stats_table['Pmax'])))

        # read in the already computed 1D PDFs
        if pdf1d_outname != None:
            print('restoring the already computed 1D PDFs from ' + pdf1d_outname)
            hdulist = fits.open(pdf1d_outname)
            for k in range(len(qnames)):
                save_pdf1d_vals[k] = hdulist[k+1].data
            hdulist.close()

        # setup the lnp file 
        if lnp_outname is not None:

#            outfile = tables.openFile(lnp_outname, 'r')
            #for group in outfile.walk_groups():
            #    print(group)
            outfile = tables.openFile(lnp_outname, 'a')
    else:
        start_pos = 0

        # setup a new lnp file
        if lnp_outname is not None:
            outfile = tables.openFile(lnp_outname, 'w')
            #Save wavelengths in root, remember #n_stars = root._v_nchildren -1
            outfile.createArray(outfile.root, 'grid_waves', g0.lamb[:])
            filters = obs.getFilters()
            outfile.createArray(outfile.root, 'obs_filters', filters[:])

    # loop over the objects and get all the requested quantities
    g0_specgrid_indx = g0['specgrid_indx']
    _p = np.asarray(p, dtype=float)

    #loop over the obs and do the work
    if hasattr(g0.seds, 'read'):
        _seds = g0.seds.read()
    else:
        _seds = g0.seds

    for e, obj in Pbar(len(obs)-start_pos, desc='Calculating Lnp/Stats').iterover(islice(obs.enumobs(),start_pos,None)):
        # calculate the full nD posterior
        (sed) = obj
        (lnp,chi2) = N_logLikelihood_NM(sed,_seds,ast_error,ast_bias,mask=None, lnp_threshold=abs(threshold) )
            
        lnp = lnp[g0_indxs]
        chi2 = chi2[g0_indxs]
        lnp = numexpr.evaluate('lnp + g0_weights')
        #lnp +=  g0_weights  # multiply by the prior weights (sum in log space)

        indx, = np.where((lnp - max(lnp[np.isfinite(lnp)])) > threshold)

        if lnp_outname is not None:
            try:
                star_group = outfile.createGroup('/', 'star_%d'  % e, title="star %d" % e)
            except tables.exceptions.NodeError:
                print('lnp for star ' + str(e) + ' already in file')
            else:
                outfile.createArray(star_group, 'input', np.array([sed]).T)
                if lnp_npts is not None:
                    rindx = np.random.choice(indx,size=lnp_npts)
                    outfile.createArray(star_group, 'idx', np.array(g0_indxs[rindx], dtype=np.int64))
                    outfile.createArray(star_group, 'lnp', np.array(lnp[rindx], dtype=np.float32))
                    outfile.createArray(star_group, 'chi2', np.array(chi2[rindx], dtype=np.float32))
                else:
                    outfile.createArray(star_group, 'idx', np.array(g0_indxs[indx], dtype=np.int64))
                    outfile.createArray(star_group, 'lnp', np.array(lnp[indx], dtype=np.float32))
                    outfile.createArray(star_group, 'chi2', np.array(chi2[indx], dtype=np.float32))
                    outfile.flush() #commit changes

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

        # incremental save (useful if job dies early to recover most of the computations)
        if save_every_npts is not None:
            if (e > 0) & (e%save_every_npts == 0):
                # save the 1D PDFs
                if pdf1d_outname is not None:
                    save_pdf1d(pdf1d_outname,save_pdf1d_vals, qnames)
    
                # save the stats/catalog file
                if stats_outname is not None:
                    save_stats(stats_outname, prev_result, best_vals, exp_vals, per_vals, chi2_vals, chi2_indx,
                               lnp_vals, lnp_indx, best_specgrid_indx, qnames, p)
    # save the 1D PDFs
    if pdf1d_outname is not None:
        save_pdf1d(pdf1d_outname,save_pdf1d_vals, qnames)
    
    # save the stats/catalog file
    if stats_outname is not None:
        save_stats(stats_outname, prev_result, best_vals, exp_vals, per_vals, chi2_vals, chi2_indx,
                   lnp_vals, lnp_indx, best_specgrid_indx, qnames, p)

def IAU_names_and_extra_info(obsdata):
    """
    generates IAU approved names for the PHAT data using RA & DEC
      and extra information about the sources (ra, dec, photometry, etc.)

    keywords
    --------
    obs: Observations

    returns
    -------
    r: dict 
        returns a dict with a (name, ndarray) pair
    """
    r = {}

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

def summary_table_memory(obs, noisemodel, sedgrid, keys=None, method=None, gridbackend='cache',
                         threshold=-10, save_every_npts=None, lnp_npts=None, resume=False,
                         stats_outname=None, pdf1d_outname=None, lnp_outname=None):
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

    for key in keys:
        if not (key in g0.keys()):
            raise KeyError('Key "{0}" not recognized'.format(key))

    # generate an IAU complient name for each source and add other inform
    res = IAU_names_and_extra_info(obs)

    Q_all_memory(res, obs, g0, noisemodel, keys, p=[16., 50., 84.], resume=resume,
                 threshold=-10.,save_every_npts=save_every_npts, lnp_npts=lnp_npts,
                 stats_outname=stats_outname,
                 pdf1d_outname=pdf1d_outname,
                 lnp_outname=lnp_outname)



