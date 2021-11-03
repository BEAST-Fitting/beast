#!/usr/bin/env python

"""
Code to create the batch files for fitting
"""
import os
import stat
import argparse
import tables

#####
import numpy as np
from astropy.table import Table
from astropy.io import fits

#####

from beast.tools import beast_settings
from beast.tools.run import create_filenames


def setup_batch_beast_fit(
    beast_settings_info,
    num_percore=5,
    nice=None,
    overwrite_logfile=True,
    prefix=None,
    use_sd=True,
    pdf2d_param_list=['Av', 'Rv', 'f_A', 'M_ini', 'logA', 'Z', 'distance'],
    nsubs=1,
    nprocs=1,
):
    """
    Sets up batch files for submission to the 'at' queue on
    linux (or similar) systems

    Parameters
    ----------
    beast_settings_info : string or beast.tools.beast_settings.beast_settings instance
        if string: file name with beast settings
        if class: beast.tools.beast_settings.beast_settings instance

    num_percore : int (default = 5)
        number of fitting runs per core

    nice : int (default = None)
        set this to an integer (-20 to 20) to prepend a "nice" level
        to the fitting command

    overwrite_logfile : boolean (default = True)
        if True, will overwrite the log file; if False, will append to
        existing log file

    prefix : string (default=None)
        Set this to a string (such as 'source activate astroconda') to prepend
        to each batch file (use '\n's to make multiple lines)

    use_sd : boolean (default=True)
        If True, split runs based on source density (determined by finding
        matches to settings.astfile with SD info)

    pdf2d_param_list : list of strings or None
        If set, do 2D PDFs of these parameters.  If None, don't make 2D PDFs.

    nsubs : int (default=1)
        number of subgrids used for the physics model

    nprocs : int (default=1)
        Number of parallel processes to use when doing the fitting
        (currently only implemented for subgrids)


    Returns
    -------
    run_info_dict : dict
        Dictionary indicating which catalog files have complete modeling, and
        which job files need to be run

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

    # setup the subdirectory for the batch and log files
    job_path = settings.project + "/fit_batch_jobs/"
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path + "logs/"
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    # get file name lists (to check if they exist and/or need to be resumed)
    file_dict = create_filenames.create_filenames(settings, use_sd=use_sd, nsubs=nsubs)

    # - input files
    photometry_files = file_dict["photometry_files"]
    # modelsedgrid_files = file_dict['modelsedgrid_files']
    # noise_files = file_dict['noise_files']

    # - output files
    stats_files = file_dict["stats_files"]
    pdf_files = file_dict["pdf_files"]
    lnp_files = file_dict["lnp_files"]

    # - total number of files
    n_files = len(photometry_files)

    # - other useful info
    sd_sub_info = file_dict["sd_sub_info"]
    gridsub_info = file_dict["gridsub_info"]

    # names of output log files
    log_files = []

    # initialize a variable name (otherwise it got auto-added in the wrong
    # place and broke the code)
    pf = None

    for i in range(n_files):

        sd_piece = ""
        if use_sd is True:
            sd_piece = "_bin" + sd_sub_info[i][0] + "_sub" + sd_sub_info[i][1]

        gridsub_piece = ""
        if nsubs > 1:
            gridsub_piece = "_gridsub" + str(gridsub_info[i])

        log_files.append("beast_fit" + sd_piece + gridsub_piece + ".log")

    # start making the job files!

    pf_open = False
    cur_f = 0
    cur_total_size = 0.0
    j = -1

    # keep track of which files are done running
    run_info_dict = {
        "phot_file": photometry_files,
        "done": np.full(n_files, False),
        "files_to_run": [],
    }

    for i, phot_file in enumerate(photometry_files):

        print("")

        # check if this is a full run
        reg_run = False
        run_done = False
        if not os.path.isfile(stats_files[i]):
            reg_run = True
            print("no stats file")
        if not os.path.isfile(pdf_files[i]):
            reg_run = True
            print("no pdf1d file")
        if not os.path.isfile(lnp_files[i]):
            reg_run = True
            print("no lnp file")

        # first check if the pdf1d mass spacing is correct
        if not reg_run:
            hdulist = fits.open(pdf_files[i])
            delta1 = hdulist["M_ini"].data[-1, 1] - hdulist["M_ini"].data[-1, 0]
            if delta1 > 1.0:  # old linear spacing
                print("pdf1d lin mass spacing - full refitting needed")
                old_mass_spacing = True
            else:
                old_mass_spacing = False
                print("pdf1d log mass spacing - ok")

            if old_mass_spacing:
                run_done = False
                reg_run = True

        # now check if the number of results is the same as
        #    the number of observations
        if not reg_run:
            # get the observed catalog
            obs = Table.read(phot_file)

            # get the fit results catalog
            t = Table.read(stats_files[i], hdu=1)
            # get the number of stars that have been fit
            (indxs,) = np.where(t["Pmax"] != 0.0)

            # get the number of entries in the lnp file
            f = tables.open_file(lnp_files[i], "r")
            nlnp = f.root._v_nchildren - 2
            f.close()

            print("# obs, stats, lnp = ", len(obs), len(indxs), nlnp)
            if (len(indxs) == len(obs)) & (nlnp == len(obs)):

                # final check, is the pdf1d file correctly populated
                tot_prob = np.sum(hdulist["M_ini"].data, axis=1)
                (tindxs,) = np.where(tot_prob > 0.0)
                print("# good pdf1d = ", len(tindxs) - 1)
                if len(tindxs) == (len(obs) + 1):
                    run_done = True

        if run_done:
            print(stats_files[i] + " done")
            run_info_dict["done"][i] = True
        else:
            j += 1
            if j % num_percore == 0:
                cur_f += 1

                # close previous files
                if j != 0:
                    pf.close()
                    # slurm needs the job file to be executable
                    #   flake8/codestyle error ignored as this if statement only executed
                    #   for j > 0 and appropriate joblist_file defined in j - 1
                    os.chmod(joblist_file, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)  # noqa: F821

                    print(
                        "total sed_trim size [Gb] = ",
                        cur_total_size / (1024.0 * 1024.0 * 1024.0),
                    )
                    cur_total_size = 0.0

                # open the slurm and param files
                pf_open = True
                joblist_file = job_path + "beast_batch_fit_" + str(cur_f) + ".joblist"
                pf = open(joblist_file, "w")
                run_info_dict["files_to_run"].append(joblist_file)

                # write out anything at the beginning of the file
                if prefix is not None:
                    pf.write(prefix + "\n")

            # flag for resuming
            resume_str = ""
            if reg_run:
                print(
                    stats_files[i]
                    + " does not exist "
                    + "- adding job as a regular fit job (not resume job)"
                )
            else:
                print(
                    stats_files[i]
                    + " not done - adding to continue fitting list ("
                    + str(len(indxs))
                    + "/"
                    + str(len(t["Pmax"]))
                    + ")"
                )
                resume_str = "-r"

            # prepend a `nice` value
            nice_str = ""
            if nice is not None:
                nice_str = "nice -n" + str(int(nice)) + " "

            # choose whether to append or overwrite log file
            pipe_str = " > "
            if not overwrite_logfile:
                pipe_str = " >> "

            # set SD+sub option
            sd_str = ""
            if use_sd is True:
                sd_str = ' --choose_sd_sub "{0}" "{1}" '.format(
                    sd_sub_info[i][0], sd_sub_info[i][1]
                )

            # set gridsub option
            gs_str = ""
            if nsubs > 1:
                gs_str = " --choose_subgrid {0} ".format(gridsub_info[i])

            # set 2D PDF option
            if pdf2d_param_list is None:
                pdf2d_str = "None"
            else:
                pdf2d_str = " " + " ".join(pdf2d_param_list) + " "

            job_command = (
                nice_str
                + "python -m beast.tools.run.run_fitting "
                + " {0} ".format(settings.settings_file)
                + resume_str
                + sd_str
                + gs_str
                + " --nsubs "
                + str(nsubs)
                + " --nprocs "
                + str(nprocs)
                + " --pdf2d_param_list "
                + pdf2d_str
                + pipe_str
                + log_path
                + log_files[i]
            )

            pf.write(job_command + "\n")

    if pf_open:
        pf.close()

        # slurm needs the job file to be executable
        os.chmod(joblist_file, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)

    # return the info about completed modeling
    return run_info_dict


if __name__ == "__main__":  # pragma: no cover

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "beast_settings_file",
        type=str,
        help="file name with beast settings",
    )
    parser.add_argument(
        "--num_percore", default=5, type=int, help="number of fitting runs per core"
    )
    parser.add_argument(
        "--nice",
        default=None,
        type=int,
        help="set this to an integer (-20 to 20) to prepend a 'nice' level to the trimming command",
    )
    parser.add_argument(
        "--overwrite_logfile",
        default=0,
        type=int,
        help="""if True (1), will overwrite the log file; if False (0), will append
        to existing log file""",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        type=str,
        help="Set this to a string to prepend to each batch file",
    )
    parser.add_argument(
        "--use_sd",
        default=1,
        type=int,
        help="if True (1), remove sources with flag=99 in flag_filter",
    )
    parser.add_argument(
        "--pdf2d_param_list",
        type=str,
        nargs="+",
        default=['Av', 'Rv', 'f_A', 'M_ini', 'logA', 'Z', 'distance'],
        help="If set, do 2D PDFs of these parameters. If None, don't make 2D PDFs."
    )
    parser.add_argument(
        "--nsubs",
        default=1,
        type=int,
        help="number of subgrids used for the physics model",
    )
    parser.add_argument(
        "--nprocs",
        default=1,
        type=int,
        help="Number of parallel processes to use when doing the fitting",
    )

    args = parser.parse_args()

    if 'None' in args.pdf2d_param_list:
        args.pdf2d_param_list = None

    setup_batch_beast_fit(
        beast_settings_info=args.beast_settings_file,
        num_percore=args.num_percore,
        nice=args.nice,
        overwrite_logfile=bool(args.overwrite_logfile),
        prefix=args.prefix,
        use_sd=bool(args.use_sd),
        pdf2d_param_list=args.pdf2d_param_list,
        nsubs=args.nsubs,
        nprocs=args.nprocs,
    )
