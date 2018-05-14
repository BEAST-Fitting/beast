#!/usr/bin/env python
"""
Script to run the BEAST on the PHAT-like data.
Assumes that the datamodel.py file exists in the same directory as this script.
  And it must be called datamodel.py
     used in make_model.py code in physicsmodel, more recoding needed to remove
     this dependency
"""

# system imports
from __future__ import (absolute_import, division, print_function)
import argparse
import time
import numpy as np
from multiprocessing import Pool

# BEAST imports
from beast.physicsmodel.create_project_dir import create_project_dir
from beast.physicsmodel.model_grid import (make_iso_table,
                                           make_spectral_grid,
                                           add_stellar_priors,
                                           make_extinguished_sed_grid)

import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel
from beast.observationmodel.ast import (make_ast_input_list,
                                        make_ast_xy_list)
from beast.fitting import fit
from beast.fitting import trim_grid
from beast.physicsmodel.grid import FileSEDGrid
from beast.tools import verify_params
from beast.tools import subgridding_tools

# I'm using this to import the datamodel module, so that this script can
# be run via a symlink
# https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path
import importlib.util
import os
datamodelfile = os.path.join(os.getcwd(), 'datamodel.py')
spec = importlib.util.spec_from_file_location('datamodel', datamodelfile)
datamodel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(datamodel)


if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--physicsmodel',
                        help='Generate the physics model grid',
                        action='store_true')
    parser.add_argument('-a', '--ast', help='Generate an input AST file',
                        action='store_true')
    parser.add_argument('-o', '--observationmodel',
                        help='Calculate the observation model (bias and noise)',
                        action='store_true')
    parser.add_argument('-t', '--trim',
                        help='Trim the physics and observation model grids',
                        action='store_true')
    parser.add_argument('-f', '--fit', help='Fit the observed data',
                        action='store_true')
    parser.add_argument('-r', '--resume', help='Resume a fitting run',
                        action='store_true')
    parser.add_argument('--nprocs', type=int, default=1,
                        help='Number of processes to use to process the subgrids')
    parser.add_argument('--nsubs', type=int, default=1,
                        help='Number of subgrids to split the physics model into')
    parser.add_argument('--subset', type=int, nargs=2, default=[None, None],
                        metavar=('START', 'STOP'),
                        help='Only process subgrids in the range [START, STOP[ '
                             '(o, t and f steps only)')
    parser.add_argument('-m', '--merge', action='store_true',
                        help='Merge the subgrid results (pdf1d and stats)')

    args = parser.parse_args()

    # check input parameters, print what is the problem, stop run_beast
    verify_params.verify_input_format(datamodel)

    def parallel_wrapper(function, argument):
        parallel = args.nprocs > 1
        if (parallel):
            p = Pool(args.nprocs)
            p.map(function, argument)
        else:
            for a in argument:
                function(a)

    outdir = os.path.join('.', datamodel.project)
    subgrid_names_file = os.path.join(outdir, 'subgrid_fnames.txt')

    def get_modelsubgridfiles():
        with open(subgrid_names_file, 'r') as f:
            modelsedgridfiles = f.read().split('\n')[:-1]
        return modelsedgridfiles

    subset_slice = slice(*args.subset)

    if args.physicsmodel:

        # make sure the project directory exists
        pdir = create_project_dir(datamodel.project)

        # download and load the isochrones
        (iso_fname, oiso) = make_iso_table(datamodel.project,
                                           oiso=datamodel.oiso,
                                           logtmin=datamodel.logt[0],
                                           logtmax=datamodel.logt[1],
                                           dlogt=datamodel.logt[2],
                                           z=datamodel.z)

        if hasattr(datamodel, 'add_spectral_properties_kwargs'):
            extra_kwargs = datamodel.add_spectral_properties_kwargs
        else:
            extra_kwargs = None

        # generate the spectral library (no dust extinction)
        (spec_fname, g_spec) = make_spectral_grid(
            datamodel.project,
            oiso,
            osl=datamodel.osl,
            distance=datamodel.distances,
            distance_unit=datamodel.distance_unit,
            add_spectral_properties_kwargs=extra_kwargs)

        # Work with the whole grid up to here (otherwise, priors need a
        # rework (they don't like having only a subset of the parameter
        # space, especially when there's only one age for example)
        (pspec_fname, g_pspec) = add_stellar_priors(datamodel.project, g_spec)

        # Make subgrids, by splitting the spectral grid into equal sized pieces
        custom_sub_pspec = subgridding_tools.split_grid(
            pspec_fname, args.nsubs)

        file_prefix = '{0}/{0}_'.format(datamodel.project)

        # process the subgrids individually
        def gen_subgrid(i, sub_name):
            sub_g_pspec = FileSEDGrid(sub_name)
            sub_seds_fname = '{}seds.gridsub{}.hd5'.format(file_prefix, i)

            (sub_seds_fname, sub_g_seds) = make_extinguished_sed_grid(
                datamodel.project,
                sub_g_pspec,
                datamodel.filters,
                extLaw=datamodel.extLaw,
                av=datamodel.avs,
                rv=datamodel.rvs,
                fA=datamodel.fAs,
                rv_prior_model=datamodel.rv_prior_model,
                av_prior_model=datamodel.av_prior_model,
                fA_prior_model=datamodel.fA_prior_model,
                add_spectral_properties_kwargs=extra_kwargs,
                seds_fname=sub_seds_fname)

            return sub_seds_fname

        par_tuples = [(i, sub_name)
                      for i, sub_name in enumerate(custom_sub_pspec)]

        parallel = args.nprocs > 1
        if (parallel):
            p = Pool(args.nprocs)
            final_sub_names = p.starmap(gen_subgrid, par_tuples)
        else:
            final_sub_names = [gen_subgrid(*pt) for pt in par_tuples]

        # Save a list of subgrid names
        with open(subgrid_names_file, 'w') as fname_file:
            for fname in final_sub_names:
                fname_file.write(fname + '\n')

        seds_fname = '{}seds.grid.hd5'.format(file_prefix)
        subgridding_tools.merge_grids(seds_fname, final_sub_names)

    if args.ast:
        # This ast section should be replaced by a more general implementation

        # Determine magnitude range for ASTs
        mag_cuts = datamodel.ast_maglimit
        if len(mag_cuts) == 1:
            tmp_cuts = mag_cuts
            obsdata = datamodel.get_obscat(datamodel.obsfile,
                                           datamodel.filters)

            min_mags = np.zeros(len(datamodel.filters))
            for k, filtername in enumerate(obsdata.filters):
                sfiltername = obsdata.data.resolve_alias(filtername)
                sfiltername = sfiltername.replace('rate', 'vega')
                sfiltername = sfiltername.replace('RATE', 'VEGA')
                keep, = np.where(obsdata[sfiltername] < 99.)
                min_mags[k] = np.percentile(obsdata[keep][sfiltername], 90.)

            # max. mags from the gst observation cat.
            mag_cuts = min_mags + tmp_cuts

        # Choose seds for ASTs
        modelsedgridfile = os.path.join(
            outdir, datamodel.project + '_seds.grid.hd5')
        outfile_chosenSEDs = os.path.join(
            outdir, datamodel.project + '_chosenSEDs.txt')
        N_per_age = datamodel.ast_models_selected_per_age
        Nfilters = datamodel.ast_bands_above_maglimit
        Nrealize = datamodel.ast_realization_per_model

        toothpick_style = True
        if toothpick_style:
            N_fluxes = 25
            min_N_per_flux = 50
            bins_outfile = os.path.join(
                outdir, datamodel.project + '_toothpick_style_bins.txt')
            chosen_seds = \
                make_ast_input_list.pick_models_toothpick_style(modelsedgridfile,
                                                                datamodel.filters,
                                                                mag_cuts,
                                                                Nfilters,
                                                                N_fluxes,
                                                                min_N_per_flux,
                                                                outfile_chosenSEDs,
                                                                bins_outfile,
                                                                mag_pad=0)
        else:
            chosen_seds = make_ast_input_list.pick_models(modelsedgridfile,
                                                          datamodel.filters,
                                                          mag_cuts,
                                                          Nfilters,
                                                          N_per_age,
                                                          Nrealize,
                                                          outfile_chosenSEDs)

        # Assign positions for ASTs
        outfile_inputAST = os.path.join(
            outdir, datamodel.project + '_inputAST.txt')
        make_ast_xy_list.pick_positions_per_background(chosen_seds,
                                                       bg_map=datamodel.ast_bg_map_file,
                                                       N_bg_bins=datamodel.ast_bg_nbins,
                                                       outfile=outfile_inputAST,
                                                       refimage=datamodel.ast_reference_image,
                                                       Nrealize=10)

    if args.observationmodel:
        print('Generating noise model from ASTs and absflux A matrix')

        modelsedgridfiles = get_modelsubgridfiles()[subset_slice]

        # Process the subgrids individually
        def gen_subobsmodel(modelsedgridfile):
            modelsedgrid = FileSEDGrid(modelsedgridfile)

            # generate the AST noise model
            noisefile = modelsedgridfile.replace('seds', 'noisemodel')
            noisemodel.make_toothpick_noise_model(
                noisefile,
                datamodel.astfile,
                modelsedgrid,
                absflux_a_matrix=datamodel.absflux_a_matrix)

        parallel_wrapper(gen_subobsmodel, modelsedgridfiles)

    if args.trim:
        print('Trimming the model and noise grids')

        # read in the observed data
        obsdata = datamodel.get_obscat(datamodel.obsfile,
                                       datamodel.filters)

        modelsedgridfiles = get_modelsubgridfiles()[subset_slice]

        # trim the models individually
        def trim_submodel(modelsedgridfile):
            modelsedgrid = FileSEDGrid(modelsedgridfile)

            # read in the noise model just created
            noisefile = modelsedgridfile.replace('seds', 'noisemodel')
            noisemodel_vals = noisemodel.get_noisemodelcat(noisefile)

            # trim the model sedgrid
            sed_trimname = modelsedgridfile.replace('seds', 'seds_trim')
            noisemodel_trimname = sed_trimname.replace('seds', 'noisemodel')

            trim_grid.trim_models(modelsedgrid, noisemodel_vals, obsdata,
                                  sed_trimname, noisemodel_trimname, sigma_fac=3.)

        parallel_wrapper(trim_submodel, modelsedgridfiles)

    if args.fit:
        start_time = time.clock()

        # read in the observed data
        obsdata = datamodel.get_obscat(datamodel.obsfile,
                                       datamodel.filters)

        modelsedgridfiles = get_modelsubgridfiles()[subset_slice]

        trimmed_modelsedgridfiles = [
            s.replace('seds', 'seds_trim') for s in modelsedgridfiles]
        trimmed_noisemodelfiles = [
            s.replace('seds', 'noisemodel') for s in trimmed_modelsedgridfiles]
        grid_info_dict = subgridding_tools.reduce_grid_info(
            trimmed_modelsedgridfiles,
            trimmed_noisemodelfiles, nprocs=args.nprocs)

        # perform fits for the subgrids individually
        def fit_submodel(modelsedgridfile):
            # input files
            trimmed_modelsedgridfile = modelsedgridfile.replace(
                'seds', 'seds_trim')
            trimmed_noisemodelfile = trimmed_modelsedgridfile.replace(
                'seds', 'noisemodel')

            # output files
            lnpfile = modelsedgridfile.replace('seds', 'lnp')
            statsfile = modelsedgridfile.replace('seds', 'stats')
            statsfile = statsfile.replace('.hd5', '.fits')
            pdf1dfile = statsfile.replace('stats', 'pdf1d')

            # load the subgrid seds and subgrid noisemodel
            modelsedgrid = FileSEDGrid(trimmed_modelsedgridfile)
            noisemodel_vals = noisemodel.get_noisemodelcat(
                trimmed_noisemodelfile)

            fit.summary_table_memory(obsdata, noisemodel_vals,
                                     modelsedgrid, resume=args.resume,
                                     threshold=-10.,
                                     save_every_npts=100, lnp_npts=60,
                                     stats_outname=statsfile,
                                     pdf1d_outname=pdf1dfile,
                                     grid_info_dict=grid_info_dict,
                                     lnp_outname=lnpfile,
                                     do_not_normalize=True)
            print('Done fitting on grid ' + trimmed_modelsedgridfile)

        parallel_wrapper(fit_submodel, modelsedgridfiles)

        new_time = time.clock()
        print('time to fit: ', (new_time - start_time) / 60., ' min')

    if args.merge:
        modelsedgridfiles = get_modelsubgridfiles()
        with_fits = [s.replace('.hd5', '.fits') for s in modelsedgridfiles]
        pdf1dfiles = [s.replace('seds', 'pdf1d') for s in with_fits]
        statsfiles = [s.replace('seds', 'stats') for s in with_fits]
        print('Merging')
        print(list(zip(pdf1dfiles, statsfiles)))
        subgridding_tools.merge_pdf1d_stats(pdf1dfiles, statsfiles)

    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()
