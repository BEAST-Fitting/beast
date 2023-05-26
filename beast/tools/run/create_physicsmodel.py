#!/usr/bin/env python

# system imports
import os
import stat
import argparse

from astropy import constants as const

# BEAST imports
from beast.physicsmodel.create_project_dir import create_project_dir
from beast.physicsmodel.model_grid import (
    make_evoltrack_table,
    make_spectral_grid,
    add_stellar_priors,
    make_extinguished_sed_grid,
)
from beast.physicsmodel.grid import SpectralGrid
from beast.tools.run.helper_functions import parallel_wrapper

# from beast.physicsmodel.stars.isochrone import ezIsoch
from beast.tools import beast_settings, subgridding_tools


def create_physicsmodel(beast_settings_info, nsubs=1, nprocs=1, subset=[None, None]):
    """
    Create the physics model grid.  If nsubs > 1, this will make sub-grids.


    Parameters
    ----------
    beast_settings_info : string or beast.tools.beast_settings.beast_settings instance
        if string: file name with beast settings
        if class: beast.tools.beast_settings.beast_settings instance

    nsubs : int (default=1)
        number of subgrids to split the physics model into

    nprocs : int (default=1)
        Number of parallel processes to use
        (currently only implemented for subgrids)

    subset : list of two ints (default=[None,None])
        Only process subgrids in the range [start,stop].
        (only relevant if nsubs > 1)

    """

    # process beast settings info
    if isinstance(beast_settings_info, str):
        settings = beast_settings.beast_settings(beast_settings_info)
    elif isinstance(beast_settings_info, beast_settings.beast_settings):
        settings = beast_settings_info
    else:
        raise TypeError(
            "beast_settings_info must be string or beast.tools.beast_settings.beast_settings instance"
        )

    # filename for the SED grid
    modelsedgrid_filename = "%s/%s_seds.grid.hd5" % (
        settings.project,
        settings.project,
    )

    # grab the current subgrid slice
    subset_slice = slice(subset[0], subset[1])

    # make sure the project directory exists
    create_project_dir(settings.project)

    # download and load the isochrones
    (iso_fname, oiso) = make_evoltrack_table(
        settings.project,
        oet=settings.oiso,
        age_info=settings.logt,
        mass_info=settings.logmass,
        z=settings.z,
    )

    if hasattr(settings, "add_spectral_properties_kwargs"):
        extra_kwargs = settings.add_spectral_properties_kwargs
    else:
        extra_kwargs = None

    if hasattr(settings, "velocity"):
        redshift = (settings.velocity / const.c).decompose().value
    else:
        redshift = 0

    # generate the spectral library (no dust extinction)
    (spec_fname, g_spec) = make_spectral_grid(
        settings.project,
        oiso,
        osl=settings.osl,
        redshift=redshift,
        distance=settings.distances,
        distance_unit=settings.distance_unit,
        extLaw=settings.extLaw,
        add_spectral_properties_kwargs=extra_kwargs,
    )

    # add the stellar priors as weights
    #   also computes the grid weights for the stellar part
    (pspec_fname, g_pspec) = add_stellar_priors(
        settings.project,
        g_spec,
        age_prior_model=settings.age_prior_model,
        mass_prior_model=settings.mass_prior_model,
        met_prior_model=settings.met_prior_model,
        distance_prior_model=settings.distance_prior_model,
    )

    # --------------------
    # no subgrids
    # --------------------

    if nsubs == 1:
        # generate the SED grid by integrating the filter response functions
        #   effect of dust extinction applied before filter integration
        #   also computes the dust priors as weights
        make_extinguished_sed_grid(
            settings.project,
            g_pspec,
            settings.filters,
            extLaw=settings.extLaw,
            av=settings.avs,
            rv=settings.rvs,
            fA=settings.fAs,
            rv_prior_model=settings.rv_prior_model,
            av_prior_model=settings.av_prior_model,
            fA_prior_model=settings.fA_prior_model,
            spec_fname=modelsedgrid_filename,
            add_spectral_properties_kwargs=extra_kwargs,
        )

    # --------------------
    # use subgrids
    # --------------------

    if nsubs > 1:
        # Work with the whole grid up to there (otherwise, priors need a
        # rework - they don't like having only a subset of the parameter
        # space, especially when there's only one age for example)

        # Make subgrids, by splitting the spectral grid into equal sized pieces
        custom_sub_pspec = subgridding_tools.split_grid(pspec_fname, nsubs)

        file_prefix = "{0}/{0}_".format(settings.project)

        # function to process the subgrids individually
        def gen_subgrid(i, sub_name):
            sub_g_pspec = SpectralGrid(sub_name)
            sub_seds_fname = "{}seds.gridsub{}.hd5".format(file_prefix, i)

            # generate the SED grid by integrating the filter response functions
            #   effect of dust extinction applied before filter integration
            #   also computes the dust priors as weights
            (sub_seds_fname, sub_g_seds) = make_extinguished_sed_grid(
                settings.project,
                sub_g_pspec,
                settings.filters,
                extLaw=settings.extLaw,
                av=settings.avs,
                rv=settings.rvs,
                fA=settings.fAs,
                rv_prior_model=settings.rv_prior_model,
                av_prior_model=settings.av_prior_model,
                fA_prior_model=settings.fA_prior_model,
                add_spectral_properties_kwargs=extra_kwargs,
                seds_fname=sub_seds_fname,
            )

            return sub_seds_fname

        # run the above function
        par_tuples = [(i, sub_name) for i, sub_name in enumerate(custom_sub_pspec)][
            subset_slice
        ]

        parallel_wrapper(gen_subgrid, par_tuples, nprocs=nprocs)

        # Save a list of subgrid names that we expect to see
        required_names = [
            "{}seds.gridsub{}.hd5".format(file_prefix, i) for i in range(nsubs)
        ]

        outdir = os.path.join(".", settings.project)
        subgrid_names_file = os.path.join(outdir, "subgrid_fnames.txt")

        with open(subgrid_names_file, "w") as fname_file:
            for fname in required_names:
                fname_file.write(fname + "\n")


def split_create_physicsmodel(beast_settings_info, nsubs=1, nprocs=1):
    """
    Making the physics model grid takes a while for production runs.  This
    creates scripts to run each subgrid as a separate job.

    Parameters
    ----------
    beast_settings_info : string or beast.tools.beast_settings.beast_settings instance
        if string: file name with beast settings
        if class: beast.tools.beast_settings.beast_settings instance

    nsubs : int (default=1)
        number of subgrids to split the physics model into

    nprocs : int (default=1)
        Number of parallel processes to use
        (currently only implemented for subgrids)

    """

    # process beast settings info
    if isinstance(beast_settings_info, str):
        settings = beast_settings.beast_settings(beast_settings_info)
    elif isinstance(beast_settings_info, beast_settings.beast_settings):
        settings = beast_settings_info
    else:
        raise TypeError(
            "beast_settings_info must be string or beast.tools.beast_settings.beast_settings instance"
        )

    # make sure the project directory exists
    create_project_dir(settings.project)

    # directory for scripts
    job_path = "./{0}/model_batch_jobs/".format(settings.project)
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path + "logs/"
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    for i in range(nsubs):

        joblist_file = job_path + "create_physicsmodel_" + str(i) + ".job"
        with open(joblist_file, "w") as jf:

            jf.write(
                "python -m beast.tools.run.create_physicsmodel "
                + " {0} ".format(settings.settings_file)
                + " --nsubs "
                + str(nsubs)
                + " --nprocs "
                + str(nprocs)
                + " --subset "
                + str(i)
                + " "
                + str(i + 1)
                + " >> "
                + log_path
                + "create_physicsmodel_"
                + str(i)
                + ".log\n"
            )

        # slurm needs it to be executable
        os.chmod(joblist_file, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)


if __name__ == "__main__":  # pragma: no cover
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "beast_settings_file", type=str, help="file name with beast settings",
    )
    parser.add_argument(
        "--nsubs",
        type=int,
        default=1,
        help="""Number of subgrids to split the physics
                        model into""",
    )
    parser.add_argument(
        "--nprocs",
        type=int,
        default=1,
        help="""Number of processes to use to process
                        the subgrids""",
    )
    parser.add_argument(
        "--subset",
        type=int,
        nargs=2,
        default=[None, None],
        help="""Only process subgrids in the range
                        [start, stop].""",
    )

    args = parser.parse_args()

    create_physicsmodel(
        beast_settings_info=args.beast_settings_file,
        nsubs=args.nsubs,
        nprocs=args.nprocs,
        subset=args.subset,
    )
