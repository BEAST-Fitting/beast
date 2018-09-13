"""
Fitting Pipeline for PHAT
BEAST Toothpick version (v1, 5 Feb 2015)
based on code by Morgan Fouesneau
major modifications by Karl Gordon (Feb-Mar 2015)
  - added a number of additional parameters
  - uses a fast 1D PDF generator
  - combines best,expectation,percentiles for speed
  - updated to compute the lnp right before it is needed
    this allows the option not to save the sparse lnp to disk
  - added the option to save a random sampling of the lnp to disk
  - code now can save the results every n stars if requested
    things allows a partially completed run to be recovered and continue
  - removed the use of ezpipe as everything is now packaged into a
    single routine and this routine often needs to be run even if
    results file already exist
  - switched from eztables to astropy.table to (potentially) avoid
    bus/memory errors
updated to allow for trunchen noise model by Karl Gordon (Dec 2015/Jan 2016)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

import sys
import time
import numpy as np
import math
import tables
import string
from itertools import islice

import numexpr

from astropy import units as ap_units
from astropy.coordinates import SkyCoord as ap_SkyCoord
#from astropy.coordinates import ICRS as ap_ICRS

from astropy.io import fits
from astropy.table import Table

from ..physicsmodel import grid
from ..tools.pbar import Pbar

from .fit_metrics.likelihood import (N_covar_logLikelihood,
                                     N_logLikelihood_NM,
                                     getNorm_lnP)
from .fit_metrics import expectation, percentile

from .pdf1d import pdf1d

__all__ = ['summary_table_memory',
           'Q_all_memory',
           'IAU_names_and_extra_info',
           'save_stats',
           'save_pdf1d',
           'save_lnp']

def save_stats(stats_outname, stats_dict_in, best_vals, exp_vals,
               per_vals, chi2_vals, chi2_indx, lnp_vals, lnp_indx,
               best_specgrid_indx, total_log_norm, qnames, p):
    """ Saves the stats to a file

    Keywords
    ----------
    stats_outname(str) : output filename
    stats_dict_in(dict) : input dictonary with ancilliary info
    best_vals(2D nparray) : best fit parameters
    exp_vals(2D nparray) : expectation fit parameters
    per_vals(3D nparray) : percentile fit parameters
    chi2_vals(1D nparray) : chisqr values (does not include model weights)
    chi2_indx(1D nparray) : indx in model grid of chisqr values
    lnp_vals(1D nparray) : P(max) values (includes model weights)
    lnp_indx(1D nparray) : indx in model grid of P(max) values
    best_specgrid_indx(1D nparray) : indx in spectroscopic model grid of
                                     P(max) values
    total_log_norm(1D nparray) : log of the total grid weight
    qnames(1D nparray) : list of the parameter names
    p(1D nparray) : list of percentiles use to create the per_vals

    Returns
    -------
    N/A
    """

    stats_dict = stats_dict_in.copy()

    # populate the dict array
    for k, qname in enumerate(qnames):
        stats_dict['{0:s}_Best'.format(qname)] = best_vals[:,k]
        stats_dict['{0:s}_Exp'.format(qname)] = exp_vals[:,k]
        for i, pval in enumerate(p):
            stats_dict['{0:s}_p{1:d}'.format(qname, int(pval))] = \
                                                    per_vals[:,k,i]

    stats_dict['chi2min'] = chi2_vals
    stats_dict['chi2min_indx'] = chi2_indx.astype(int)
    stats_dict['Pmax'] = lnp_vals
    stats_dict['Pmax_indx'] = lnp_indx.astype(int)
    stats_dict['specgrid_indx'] = best_specgrid_indx.astype(int)
    stats_dict['total_log_norm'] = total_log_norm

    summary_tab = Table(stats_dict)

    if stats_outname is not None:
        summary_tab.write(stats_outname, overwrite=True)

def save_pdf1d(pdf1d_outname, save_pdf1d_vals, qnames):
    """ Saves the 1D PDFs to a file

    Keywords
    ----------
    pdf1d_outname(str) : output filename
    save_pdf1d_vals(list) : list of 2D nparrays giving the 1D PDFs for
                            each parameter/variable
    qnames(1D nparray) : list of the parameter names

    Returns
    -------
    N/A
    """

    # write a small primary header
    fits.writeto(pdf1d_outname, np.zeros((2,2)), overwrite=True)

    # write the 1D PDFs for all the objects, 1 set per extension
    for k, qname in enumerate(qnames):
        hdu = fits.PrimaryHDU(save_pdf1d_vals[k])
        pheader = hdu.header
        pheader.set('XTENSION','IMAGE') 
        pheader.set('EXTNAME',qname) 
        fits.append(pdf1d_outname, save_pdf1d_vals[k], header=pheader)

def save_lnp(lnp_outname, save_lnp_vals, resume):
    """ Saves the nD lnps to a file

    Keywords
    ----------
    lnp_outname(str) : output filename
    save_lnp_vals(list) : list of 5 parameter lists giving the lnp/chisqr
                          info for each star
    resume(boolean) : **not used** remove later

    Returns
    -------
    N/A
    """

    # code needed if hdf5 is corrupted - usually due to job ending in the
    #    middle of the writing of the lnp file
    #  should be rare (not originally as the lnp file was open and
    #    written to continuously - 
    #                  should be fixed with the new code where the lnp
    #                  is saved every n stars instead)
    try:
        outfile = tables.open_file(lnp_outname, 'a')
    except Exception as error:
        print('partial run lnp file is corrupted - saving new lnp values in '
              + string.replace(lnp_outname,'lnp','lnp_partial'))
        outfile = tables.open_file(string.replace(lnp_outname,'lnp',
                                                 'lnp_partial'), 'a')
            
    for lnp_val in save_lnp_vals:
        e = lnp_val[0]
        try:
            star_group = outfile.create_group('/', 'star_%d'  % e,
                                             title="star %d" % e)
        except tables.exceptions.NodeError:
            #print('lnp for star ' + str(e) + ' already in file')
            pass
        else:
            outfile.create_array(star_group, 'input', lnp_val[4])
            outfile.create_array(star_group, 'idx', lnp_val[1])
            outfile.create_array(star_group, 'lnp', lnp_val[2])
            outfile.create_array(star_group, 'chi2', lnp_val[3])
    outfile.close()

def Q_all_memory(prev_result, obs, sedgrid, ast, qnames_in, p=[16., 50., 84.],
                 gridbackend='cache', max_nbins=50,
                 stats_outname=None, pdf1d_outname=None, grid_info_dict=None,
                 lnp_outname=None, lnp_npts=None, save_every_npts=None,
                 threshold=-40, resume=False,
                 use_full_cov_matrix=True, do_not_normalize=False):
    """ Fit each star, calculate various fit statistics, and output them
        to files
      (done in one function for speed and ability to resume partially
      completed runs)

    keywords
    --------
    prev_result: dict
        previous results to include in the output summary table
        usually basic data on each source

    obs: Observation object instance
        observation catalog

    sedgrid: str or grid.SEDgrid instance
        model grid

    ast: beast noisemodel instance
        noise model data

    qnames: list of quantities or expressions

    p: array-like
        list of percentile values

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    max_nbins: maxiumum number of bins to use for the 1D likelihood
               calculations

    save_every_npts: integer
        set to save the files below (if set) every n stars
        a requirement for recovering from partially complete runs

    resume: boolean
        set to designate this run is resuming a partially complete run

    use_full_cov_matrix: boolean
        set to use the full covariance matrix if it is present in the
        noise model file

    stats_outname: set to output the stats file into a FITS file with
                   extensions

    pdf1d_outname: set to output the 1D PDFs into a FITS file with extensions

    grid_info_dict: dict: {'qname': {'min': float,
                                     'max': float,
                                     'num_unique': int},
                           ...}
        Set to override the mins/maxes of the 1dpdfs, and the number of
        unique values.

    lnp_outname: set to output the sparse likelihoods into a (usually HDF5)
                 file

    threshold: value above which to use/save for the lnps
               (defines the sparse likelihood)

    lnp_npts: set to a number to output a random sampling of the
              lnp points above the threshold
              otherwise, the full sparse likelihood is output

    do_not_normalize: bool
        Do not normalize the prior weights before applying them. This
        should have no effect on the final outcome when using only a
        single grid, but is essential when using the subgridding
        approach.

    returns
    -------
    N/A
    """

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid, backend=gridbackend)
    else:
        g0 = sedgrid

    # remove weights that are less than zero
    g0_indxs, = np.where(g0['weight'] > 0.0)

    g0_weights = np.log(g0['weight'][g0_indxs])
    if not do_not_normalize:
        g0_weights_sum = np.log(g0['weight'][g0_indxs].sum())
        g0_weights = numexpr.evaluate("g0_weights - g0_weights_sum")

    if len(g0['weight']) != len(g0_indxs):
        print('some zero weight models exist')
        print('orig/g0_indxs', len(g0['weight']),len(g0_indxs))

    # get the model SEDs
    if hasattr(g0.seds, 'read'):
        _seds = g0.seds.read()
    else:
        _seds = g0.seds

    # get the names of all the children in the ast structure
    ast_children = []
    for label, node in list(ast.root._v_children.items()): 
        ast_children.append(label)

    # links to errors and biases
    ast_error = ast.root.error[:]
    ast_bias = ast.root.bias[:]

    # if the ast file includes the full covariance matrices, make links
    full_cov_mat = False
    if (use_full_cov_matrix &
        ('q_norm' in ast_children) &
        ('icov_diag' in ast_children) &
        ('icov_offdiag' in ast_children)):
        full_cov_mat = True
        ast_q_norm = ast.root.q_norm[:]
        ast_icov_diag = ast.root.icov_diag[:]
        ast_icov_offdiag = ast.root.icov_offdiag[:]

    if full_cov_mat:
        print('using full covariance matrix')
    else:
        print('not using full covariance matrix')

    # number of observed SEDs to fit
    nobs = len(obs)

    # augment the qnames to include the *full* model SED
    #  by this it means the physical model flux plus the noise model bias term
    qnames = qnames_in
    filters = sedgrid.filters
    for i, cfilter in enumerate(filters):
        qnames.append('symlog'+cfilter+'_wd_bias')

    # create the full model fluxes for later use
    #   save as symmetric log, since the fluxes can be negative
    full_model_flux = _seds + ast_bias
    logtempseds = np.array(full_model_flux)
    #full_model_flux = np.sign(logtempseds) * np.log10(1 + np.abs(logtempseds * math.log(10)))
    full_model_flux = np.sign(logtempseds) * np.log1p(np.abs(logtempseds * math.log(10)))/math.log(10)

    # setup the arrays to temp store the results
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
    total_log_norm = np.zeros(nobs)

    # variable to save the lnp files
    save_lnp_vals = []
    
    # setup the mapping for the 1D PDFs
    fast_pdf1d_objs = []
    save_pdf1d_vals = []

    for qname in qnames:
        #q = g0[qname][g0_indxs]
        if '_bias' in qname:
            fname = (qname.replace('_wd_bias','')).replace('symlog','')
            q = full_model_flux[:,filters.index(fname)]
        else:
            q = g0[qname]

        if grid_info_dict is not None and qname in grid_info_dict:
            # When processing a subgrid, we actuall need the number of
            # unique values across all the subgrids to make the 1dpdfs
            # compatible
            n_uniq = grid_info_dict[qname]['num_unique']
        else:
            n_uniq = len(np.unique(q))
        
        if n_uniq > max_nbins:
            # limit the number of bins in the 1D likelihood for speed
            nbins = max_nbins  
        else:
            nbins = n_uniq

        # temp code for BEAST paper figure
        if qname == 'Z':
            nbins = nbins + 1

        # setup the fast 1d pdf

        # needed for mass parameters as they are stored as linear values
        # computationally, less bins needed if 1D PDFs done as log spacing
        if qname in set(['M_ini', 'M_act','radius']):
            logspacing = True
        else:
            logspacing = False

        if grid_info_dict is not None and qname in grid_info_dict:
            minval = grid_info_dict[qname]['min']
            maxval = grid_info_dict[qname]['max']
        else:
            minval = None
            maxval = None

        # generate the fast 1d pdf mapping
        _tpdf1d = pdf1d(q, nbins,
                        logspacing=logspacing, minval=minval,
                        maxval=maxval)
        fast_pdf1d_objs.append(_tpdf1d)
        
        # setup the arrays to save the 1d PDFs
        save_pdf1d_vals.append(np.zeros((nobs+1, nbins)))
        save_pdf1d_vals[-1][nobs,:] = _tpdf1d.bin_vals

    # if this is a resume job, read in the already computed stats and
    #     fill the variables
    # also - find the start position for the resumed run
    if resume:
        stats_table = Table.read(stats_outname)
        
        for k, qname in enumerate(qnames):
            best_vals[:,k] = stats_table['{0:s}_Best'.format(qname)]
            exp_vals[:,k] = stats_table['{0:s}_Exp'.format(qname)]
            for i, pval in enumerate(p):
                per_vals[:,k,i] = stats_table['{0:s}_p{1:d}'.format(qname,
                                                                    int(pval))]
                
        chi2_vals = stats_table['chi2min']
        chi2_indx = stats_table['chi2min_indx']
        lnp_vals = stats_table['Pmax']
        lnp_indx = stats_table['Pmax_indx']
        best_specgrid_indx = stats_table['specgrid_indx']

        indxs, = np.where(stats_table['Pmax'] != 0.0)
        start_pos = max(indxs) + 1
        print('resuming run with start indx = ' + str(start_pos) +
              ' out of ' + str(len(stats_table['Pmax'])))

        # read in the already computed 1D PDFs
        if pdf1d_outname != None:
            print('restoring the already computed 1D PDFs from ' +
                  pdf1d_outname)
            hdulist = fits.open(pdf1d_outname)
            for k in range(len(qnames)):
                save_pdf1d_vals[k] = hdulist[k+1].data
            hdulist.close()
    else:
        start_pos = 0

        # setup a new lnp file
        if lnp_outname is not None:
            outfile = tables.open_file(lnp_outname, 'w')
            #Save wavelengths in root, remember #n_stars = root._v_nchildren -1
            outfile.create_array(outfile.root, 'grid_waves', g0.lamb[:])
            filters = obs.getFilters()
            outfile.create_array(outfile.root, 'obs_filters', filters[:])
            outfile.close()

    # loop over the objects and get all the requested quantities
    g0_specgrid_indx = g0['specgrid_indx']
    _p = np.asarray(p, dtype=float)

    it = Pbar(len(obs)-start_pos,
              desc='Calculating Lnp/Stats').iterover(islice(obs.enumobs(),
                                                            int(start_pos),None))
    for e, obj in it:
        # calculate the full nD posterior
        (sed) = obj

        cur_mask = (sed == 0)
        # need an alternate way to generate the mask as zeros can be
        # valid values in the observed SED (KDG 29 Jan 2016)
        # currently, set mask to False always
        cur_mask[:] = False
        
        if full_cov_mat:
            (lnp, chi2) = N_covar_logLikelihood(sed, _seds,
                                                ast_bias,
                                                ast_q_norm,
                                                ast_icov_diag,
                                                ast_icov_offdiag,
                                                lnp_threshold=abs(threshold))
        else:
            (lnp,chi2) = N_logLikelihood_NM(sed,_seds,ast_error,ast_bias,
                                            mask=cur_mask,
                                            lnp_threshold=abs(threshold) )
            
        lnp = lnp[g0_indxs]
        chi2 = chi2[g0_indxs]
        #lnp = numexpr.evaluate('lnp + g0_weights')
        lnp +=  g0_weights  # multiply by the prior weights (sum in log space)

        indx, = np.where((lnp - max(lnp[np.isfinite(lnp)])) > threshold)

        # now generate the sparse likelihood (remove later if this works
        #       by updating code below)
        #   checked if changing to the full likelihood speeds things up
        #       - the answer is no
        #   and is likely related to the switch here to the sparse
        #       likelihood for the weight calculation
        lnps = lnp[indx]
        chi2s = chi2[indx]

        #log_norm = np.log(getNorm_lnP(lnps))
        #if not np.isfinite(log_norm):
        #    log_norm = lnps.max()
        log_norm = lnps.max()
        weights = np.exp(lnps - log_norm)

        # normalize the weights make sure they sum to one
        #   needed for np.random.choice
        weight_sum = np.sum(weights)
        weights /= weight_sum

        # save the current set of lnps
        if lnp_outname is not None:
            if lnp_npts is not None:
                if lnp_npts < len(indx):
                    rindx = np.random.choice(indx, size=lnp_npts, replace=False)
                if lnp_npts >= len(indx):
                    rindx = indx
            else:
                rindx = indx
            save_lnp_vals.append([e,
                                  np.array(g0_indxs[rindx], dtype=np.int64),
                                  np.array(lnp[rindx], dtype=np.float32),
                                  np.array(chi2[rindx], dtype=np.float32),
                                  np.array([sed]).T])

        # To merge the stats for different subgrids, we need the total
        # weight of a grid, which is sum(exp(lnps)). Since sum(exp(lnps
        # - log_norm - log(weight_sum))) = 1, the relative weight of
        # each subgrid will be exp(log_norm + log(weight_sum)).
        # Therefore, we also store the following quantity:
        total_log_norm[e] = log_norm + np.log(weight_sum)

        # index to the full model grid for the best fit values
        best_full_indx = g0_indxs[indx[weights.argmax()]]

        # index to the spectral grid 
        best_specgrid_indx[e] = g0_specgrid_indx[best_full_indx]
            
        # goodness of fit quantities
        chi2_vals[e] = chi2s.min()
        chi2_indx[e] = g0_indxs[indx[chi2s.argmin()]]
        lnp_vals[e] = lnps.max()
        lnp_indx[e] = best_full_indx

        for k, qname in enumerate(qnames):
            if '_bias' in qname:
                fname = (qname.replace('_wd_bias','')).replace('symlog','')
                q = full_model_flux[:,filters.index(fname)]
            else:
                q = g0[qname]

            # best value
            best_vals[e,k] = q[best_full_indx]

            # expectation value
            exp_vals[e,k] = expectation(q[g0_indxs[indx]], weights=weights)

            # percentile values
            pdf1d_bins, pdf1d_vals = fast_pdf1d_objs[k].gen1d(g0_indxs[indx],
                                                              weights)

            save_pdf1d_vals[k][e,:] = pdf1d_vals
            if pdf1d_vals.max() > 0:
                # remove normalization to allow for post processing with 
                #   different distance runs (needed for the SMIDGE-SMC)
                #pdf1d_vals /= pdf1d_vals.max()
                per_vals[e,k,:] = percentile(pdf1d_bins, _p,
                                             weights=pdf1d_vals)
            else:
                per_vals[e,k,:] = [0.0,0.0,0.0]

        # incremental save (useful if job dies early to recover most
        #    of the computations)
        if save_every_npts is not None:
            if (e > 0) & (e%save_every_npts == 0):
                # save the 1D PDFs
                if pdf1d_outname is not None:
                    save_pdf1d(pdf1d_outname,save_pdf1d_vals, qnames)
    
                # save the stats/catalog
                if stats_outname is not None:
                    save_stats(stats_outname, prev_result, best_vals,
                               exp_vals, per_vals, chi2_vals, chi2_indx,
                               lnp_vals, lnp_indx, best_specgrid_indx,
                               total_log_norm, qnames, p)

                # save the lnps
                if lnp_outname is not None:
                    save_lnp(lnp_outname, save_lnp_vals, resume)
                    save_lnp_vals = []

    ## do the final save of everything (or the last set for the lnp values)

    # save the 1D PDFs
    if pdf1d_outname is not None:
        save_pdf1d(pdf1d_outname,save_pdf1d_vals, qnames)
    
    # save the stats/catalog
    if stats_outname is not None:
        save_stats(stats_outname, prev_result, best_vals, exp_vals,
                   per_vals, chi2_vals, chi2_indx, lnp_vals, lnp_indx,
                   best_specgrid_indx, total_log_norm, qnames, p)

    # save the lnps
    if lnp_outname is not None:
        save_lnp(lnp_outname, save_lnp_vals, resume)

def IAU_names_and_extra_info(obsdata, surveyname='PHAT',extraInfo=False):
    """
    generates IAU approved names for the data using RA & DEC
      and extra information about the sources (ra, dec, photometry, etc.)

    keywords
    --------
    obsdata: Observations
    surveyname: string
          name of survey [default = 'PHAT']
    extraInfo: bool
          set to get the HST specific PHAT software reduced survey information

    returns
    -------
    r: dict 
        returns a dict with a (name, ndarray) pair
    """
    r = {}

    go_name = False
    if 'ra' in list(obsdata.data.keys()):
        go_name = True
        ra_str = 'ra'
        dec_str = 'dec'

    if 'RA' in list(obsdata.data.keys()):
        go_name = True
        ra_str = 'RA'
        dec_str = 'DEC'

    if go_name:
        # generate the IAU names
        _tnames = []
        for i in range(len(obsdata)):
            c = ap_SkyCoord(ra=obsdata.data[ra_str][i]*ap_units.degree,
                            dec=obsdata.data[dec_str][i]*ap_units.degree,
                            frame='icrs')
            _tnames.append(surveyname + ' J' + 
                           c.ra.to_string(unit=ap_units.hourangle, 
                                          sep="",precision=2,
                                          alwayssign=False,pad=True) + 
                           c.dec.to_string(sep="",precision=2,
                                           alwayssign=True,pad=True))
            r['Name'] = _tnames

            # other useful information
            r['RA'] = obsdata.data[ra_str]
            r['DEC'] = obsdata.data[dec_str]
            if extraInfo:
                r['field'] = obsdata.data['field']
                r['inside_brick'] = obsdata.data['inside_brick']
                r['inside_chipgap'] = obsdata.data['inside_chipgap']
    else:
        r['Name'] = ["noname" for x in range(len(obsdata))]

    # include the observed filter fluxes
    for k, filtername in enumerate(obsdata.filters): 
        r[filtername] = (obsdata.data[filtername]*
                         obsdata.vega_flux[k]).astype(float) 

    return r

def summary_table_memory(obs, noisemodel, sedgrid, keys=None,
                         gridbackend='cache', threshold=-10,
                         save_every_npts=None, lnp_npts=None,
                         resume=False, stats_outname=None,
                         pdf1d_outname=None, grid_info_dict=None,
                         lnp_outname=None, use_full_cov_matrix=True,
                         surveyname='PHAT', extraInfo=False,
                         do_not_normalize=False):
    """
    keywords
    --------

    obs: Observation object instance
        observation catalog

    noisemodel: beast noisemodel instance
        noise model data

    sedgrid: str or grid.SEDgrid instance
        model grid

    keys: str or list of str
        if str:  name of the quantity or expression to evaluate from the
                 grid table
        if list: list of quantities or expresions

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    save_every_npts: integer
        set to save the files below (if set) every n stars
        a requirement for recovering from partially complete runs

    resume: boolean
        set to designate this run is resuming a partially complete run

    use_full_cov_matrix: boolean
        set to use the full covariance matrix if it is present in the
        noise model file

    stats_outname: set to output the stats file into a FITS file with
                   extensions

    pdf1d_outname: set to output the 1D PDFs into a FITS file with extensions

    grid_info_dict: dict: {'qname': {'min': float,
                                     'max': float,
                                     'num_unique': int},
                           ...}
        Set to override the mins/maxes of the 1dpdfs, and the number of
        unique values.

    lnp_outname: set to output the sparse likelihoods into a (usually HDF5)
                 file

    threshold: value above which to use/save for the lnps (defines the
               sparse likelihood)

    lnp_npts: set to a number to output a random sampling of the
              lnp points above the threshold
              otherwise, the full sparse likelihood is output

    surveyname: string
          name of survey [default = 'PHAT']

    extraInfo: bool
        set to get extra information, such as IAU name, brick, field, etc.

    do_not_normalize: bool
        Do not normalize the prior weights before applying them. This
        should have no effect on the final outcome when using only a
        single grid, but is essential when using the subgridding
        approach.

    returns
    -------
    N/A
    """

    if type(sedgrid) == str:
        g0 = grid.FileSEDGrid(sedgrid, backend=gridbackend)
    else:
        g0 = sedgrid

    if keys is None:
        keys = list(g0.keys())

    #make sure keys are real keys
    skip_keys = 'osl keep weight grid_weight prior_weight fullgrid_idx stage specgrid_indx'.split()
    keys = [k for k in keys if k not in skip_keys]

    for key in keys:
        if not (key in list(g0.keys())):
            raise KeyError('Key "{0}" not recognized'.format(key))

    # generate an IAU complient name for each source and add other inform
    res = IAU_names_and_extra_info(obs, 
                                   surveyname=surveyname,
                                   extraInfo=False)

    Q_all_memory(res, obs, g0, noisemodel, keys, p=[16., 50., 84.],
                 resume=resume,
                 threshold=threshold,save_every_npts=save_every_npts,
                 lnp_npts=lnp_npts,
                 stats_outname=stats_outname,
                 pdf1d_outname=pdf1d_outname,
                 grid_info_dict=grid_info_dict,
                 lnp_outname=lnp_outname,
                 use_full_cov_matrix=use_full_cov_matrix,
                 do_not_normalize=do_not_normalize)
