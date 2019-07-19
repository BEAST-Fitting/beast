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

import pickle


# import datamodel
# print("Import statement gets datamodel from {}.".format(datamodel.__file__))
# print("Distances are {}".format(datamodel.distances))

# import runpy, os

# print("runpy.run_path runs code from {}.".format(datamodelfile))
# datamodel_globals = runpy.run_path(datamodelfile)
# print("Distances are {}".format(datamodel_globals['distances']))

# https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path
import importlib.util
import os
datamodelfile = os.path.join(os.getcwd(), 'datamodel.py')
spec = importlib.util.spec_from_file_location('datamodel', datamodelfile)
datamodel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(datamodel)
print("Project name: {}".format(datamodel.project))

outdir = os.path.join('.', datamodel.project)
subgrid_names_file = os.path.join(outdir, 'subgrid_fnames.txt')

if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--physicsmodel",
                        help="Generate the physics model grid",
                        action="store_true")
    parser.add_argument("-a", "--ast", help="Generate an input AST file",
                        action="store_true")
    parser.add_argument("-o", "--observationmodel",
                        help="Calculate the observation model (bias and noise)",
                        action="store_true")
    parser.add_argument("-t", "--trim",
                        help="Trim the physics and observation model grids",
                        action="store_true")
    parser.add_argument("-f", "--fit", help="Fit the observed data",
                        action="store_true")
    parser.add_argument("-r", "--resume", help="Resume a fitting run",
                        action="store_true")
    parser.add_argument("--nprocs", type=int, default=1,
                        help='''Number of processes to use to process
                        the subgrids''')
    parser.add_argument("--nsubs", type=int, default=1,
                        help='''Number of subgrids to split the physics
                        model into''')
    parser.add_argument("--subset", type=int, nargs=2, default=[None, None],
                        help='''Only process subgrids in the range
                        [start, stop[. Should work for the o, t and f
                        steps''')
    parser.add_argument("-m", "--merge", action="store_true",
                        help="Merge the subgrid results")
    parser.add_argument("--ignore-missing-subresults", action="store_true",
                        help='''In some cases, it might not be possible
                        to perform the fit step for some subgrids, for
                        example because the trimmed grids are empty or
                        only have zero/negative weight points. In that
                        case, this option can be used to perform the
                        merge anyway. Just make sure that any missing
                        grids are missing for the right reasons.''')
    parser.add_argument("--dens_bin", type=int, default=None,
                        help='''Run for a certain source/background
                        density bin. Use the split_catalog_using_map `
                        tool to generate the right files for this
                        option.''')

    args = parser.parse_args()

    if args.dens_bin is not None:
        bin_subfolder = 'bin{}'.format(args.dens_bin)
        pdir = create_project_dir(bin_subfolder)

    def subcatalog_fname(full_cat_fname, dens_bin):
        return full_cat_fname.replace('.fits', '_bin{}.fits'.format(dens_bin))

    # check input parameters, print what is the problem, stop run_beast
    verify_params.verify_input_format(datamodel)

    def parallel_wrapper(function, argument):
        parallel = args.nprocs > 1
        if (parallel):
            p = Pool(args.nprocs)
            for r in p.imap_unordered(function, argument, chunksize=1):
                print(r)

        else:
            for a in argument:
                r = function(a)
                print(r)

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
            extLaw=datamodel.extLaw,
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
                      for i, sub_name in enumerate(custom_sub_pspec)][subset_slice]

        parallel = args.nprocs > 1
        if (parallel):
            p = Pool(args.nprocs)
            p.starmap(gen_subgrid, par_tuples)
        else:
            for pt in par_tuples:
                gen_subgrid(*pt)

        # Save a list of subgrid names that we expect to see
        required_names = ['{}seds.gridsub{}.hd5'.format(file_prefix, i)
                          for i in range(args.nsubs)]
        with open(subgrid_names_file, 'w') as fname_file:
            for fname in required_names:
                fname_file.write(fname + '\n')

        # seds_fname = '{}seds.grid.hd5'.format(file_prefix)
        # subgridding_tools.merge_grids(seds_fname, final_sub_names)

    if args.ast:
        # Determine magnitude range for ASTs
        mag_cuts = datamodel.ast_maglimit
        bright_cuts = None
        if len(mag_cuts) == 1:
            tmp_cuts = mag_cuts
            obsdata = datamodel.get_obscat(datamodel.obsfile,
                                           datamodel.filters)

            faintest_mags = np.zeros(len(datamodel.filters))
            brightest_mags = np.zeros(len(datamodel.filters))
            for k, filtername in enumerate(obsdata.filters):
                sfiltername = obsdata.data.resolve_alias(filtername)
                sfiltername = sfiltername.replace('rate', 'vega')
                sfiltername = sfiltername.replace('RATE', 'VEGA')
                keep, = np.where(obsdata[sfiltername] < 99.)
                faintest_mags[k] = np.percentile(obsdata[keep][sfiltername], 90.)
                brightest_mags[k] = np.amin(obsdata[keep][sfiltername])

            # max. mags from the gst observation cat.
            mag_cuts = faintest_mags + tmp_cuts # this many mags fainter thant the 90th percentile
            bright_cuts = brightest_mags - tmp_cuts # this many mags brighter than the brightest source

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
                                                                bright_cut=bright_cuts)
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
        make_ast_xy_list.pick_positions_from_map(chosen_seds,
                                                 input_map=datamodel.ast_bg_map_file,
                                                 N_bins=datamodel.ast_bg_nbins,
                                                 Npermodel=10,
                                                 outfile=outfile_inputAST,
                                                 refimage=datamodel.ast_reference_image,
                                                 Nrealize=1)

    if args.observationmodel:
        print('Generating noise model from ASTs and absflux A matrix')

        modelsedgridfiles = get_modelsubgridfiles()[subset_slice]

        # Process the subgrids individually
        def gen_subobsmodel(modelsedgridfile):
            modelsedgrid = FileSEDGrid(modelsedgridfile)

            # generate the AST noise model
            noisefile = modelsedgridfile.replace('seds', 'noisemodel')
            astfile = datamodel.astfile

            # If we are treating regions with different
            # backgrounds/source densities separately, pick one of the
            # split ast files, and put the results in a subfolder.
            if args.dens_bin is not None:
                noisefile = os.path.join(bin_subfolder, noisefile)
                create_project_dir(os.path.dirname(noisefile))
                astfile = subcatalog_fname(astfile, args.dens_bin)

            outname = noisemodel.make_toothpick_noise_model(
                noisefile,
                astfile,
                modelsedgrid,
                absflux_a_matrix=datamodel.absflux_a_matrix)

            return outname

        parallel_wrapper(gen_subobsmodel, modelsedgridfiles)

    if args.trim:
        print('Trimming the model and noise grids')

        # read in the observed data
        obsfile = datamodel.obsfile
        if args.dens_bin is not None:
            obsfile = subcatalog_fname(obsfile, args.dens_bin)
        obsdata = datamodel.get_obscat(obsfile, datamodel.filters)

        modelsedgridfiles = get_modelsubgridfiles()[subset_slice]

        # trim the models individually
        def trim_submodel(modelsedgridfile):
            modelsedgrid = FileSEDGrid(modelsedgridfile)


            noisefile = modelsedgridfile.replace('seds', 'noisemodel')
            sed_trimname = modelsedgridfile.replace('seds', 'seds_trim')
            noisemodel_trimname = sed_trimname.replace('seds', 'noisemodel')
            # When working with density bins, we nees to work in a subfolder
            if args.dens_bin is not None:
                noisefile = os.path.join(bin_subfolder, noisefile)
                sed_trimname = os.path.join(bin_subfolder, sed_trimname)
                noisemodel_trimname = os.path.join(bin_subfolder, noisemodel_trimname)

            # read in the noise model just created
            noisemodel_vals = noisemodel.get_noisemodelcat(noisefile)

            # trim the model sedgrid
            trim_grid.trim_models(modelsedgrid, noisemodel_vals, obsdata,
                                  sed_trimname, noisemodel_trimname, sigma_fac=3.)

        parallel_wrapper(trim_submodel, modelsedgridfiles)

    if args.fit:
        start_time = time.clock()

        # read in the observed data
        obsfile = datamodel.obsfile
        if args.dens_bin is not None:
            obsfile = subcatalog_fname(obsfile, args.dens_bin)

        obsdata = datamodel.get_obscat(obsfile, datamodel.filters)

        modelsedgridfiles = get_modelsubgridfiles()
        trimmed_modelsedgridfiles = [
            s.replace('seds', 'seds_trim') for s in modelsedgridfiles]
        trimmed_noisemodelfiles = [
            s.replace('seds', 'noisemodel') for s in trimmed_modelsedgridfiles]

        # File where the ranges and number of unique values for the grid
        # will be stored (this can take a while to calculate)
        grid_info_pkl = 'grid_info_dict.pkl'

        if args.dens_bin is not None:
            # Use the right subfolder
            trimmed_modelsedgridfiles, trimmed_noisemodelfiles = [
                [os.path.join(bin_subfolder, f) for f in l]
                for l in [trimmed_modelsedgridfiles, trimmed_noisemodelfiles]
            ]
            grid_info_pkl = os.path.join(bin_subfolder, grid_info_pkl)

        if not os.path.isfile(grid_info_pkl):
            grid_info_dict = subgridding_tools.reduce_grid_info(
                trimmed_modelsedgridfiles,
                trimmed_noisemodelfiles, nprocs=4)

            with open(grid_info_pkl, 'wb') as p:
                pickle.dump(grid_info_dict, p)
            print('wrote grid_info_dict to ' + grid_info_pkl)
        else:
            print('loading grid_info_dict from ' + grid_info_pkl)
            with open(grid_info_pkl, 'rb') as p:
                grid_info_dict = pickle.loads(p.read())

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

            if args.dens_bin is not None:
                # Put everything in the right subfolder
                trimmed_modelsedgridfile, trimmed_noisemodelfile, lnpfile, statsfile, pdf1dfile = [
                    os.path.join(bin_subfolder, f)
                    for f in [trimmed_modelsedgridfile, trimmed_noisemodelfile, lnpfile, statsfile,
                              pdf1dfile]
                ]


            # load the subgrid seds and subgrid noisemodel
            modelsedgrid = FileSEDGrid(trimmed_modelsedgridfile)
            noisemodel_vals = noisemodel.get_noisemodelcat(
                trimmed_noisemodelfile)

            try:
                fit.summary_table_memory(obsdata, noisemodel_vals,
                                         modelsedgrid, resume=args.resume,
                                         threshold=-10.,
                                         save_every_npts=100, lnp_npts=60,
                                         stats_outname=statsfile,
                                         pdf1d_outname=pdf1dfile,
                                         grid_info_dict=grid_info_dict,
                                         lnp_outname=lnpfile,
                                         do_not_normalize=True)
                print("Done fitting on grid " + trimmed_modelsedgridfile)
            except Exception as e:
                if not args.ignore_missing_subresults:
                    raise e

        parallel_wrapper(fit_submodel, modelsedgridfiles[subset_slice])

        new_time = time.clock()
        print('time to fit: ', (new_time - start_time) / 60., ' min')

    if args.merge:
        modelsedgridfiles = get_modelsubgridfiles()
        with_fits = [s.replace('.hd5', '.fits') for s in modelsedgridfiles]
        pdf1dfiles = [s.replace('seds', 'pdf1d') for s in with_fits]
        statsfiles = [s.replace('seds', 'stats') for s in with_fits]
        output_fname_base = os.path.join(datamodel.project, 'combined')
        if args.dens_bin is not None:
            pdf1dfiles, statsfiles = [[os.path.join(bin_subfolder, f) for f in l]
                                      for l in [pdf1dfiles, statsfiles]]
            output_fname_base = os.path.join(bin_subfolder, output_fname_base)

        if args.ignore_missing_subresults:
            # remove any missing filenames from the lists, and hope for the best
            def only_existing_files(file_list):
                return [f for f in file_list if os.path.isfile(f)]
            pdf1dfiles = only_existing_files(pdf1dfiles)
            statsfiles = only_existing_files(statsfiles)

        print("Merging")
        print(list(zip(pdf1dfiles, statsfiles)))

        subgridding_tools.merge_pdf1d_stats(pdf1dfiles, statsfiles,
                                            output_fname_base=output_fname_base)

    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()
