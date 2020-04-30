# system imports
import os
import argparse
import time
import pickle


# BEAST imports
from beast.fitting import fit
from beast.physicsmodel.grid import FileSEDGrid
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel
from beast.tools import verify_params
from beast.tools.run.helper_functions import parallel_wrapper
from beast.tools import subgridding_tools
from beast.tools.run import create_filenames


import datamodel
import importlib


def run_fitting(
    use_sd=True,
    nsubs=1,
    nprocs=1,
    choose_sd_sub=None,
    choose_subgrid=None,
    pdf2d_param_list=['Av', 'Rv', 'f_A', 'M_ini', 'logA', 'Z', 'distance'],
    resume=False,
):
    """
    Run the fitting.  If nsubs > 1, this will find existing subgrids.
    If use_sd is True, will also incorporate source density info.

    The additional choose_* options are to make queue scripts usable,
    by specifying a given SD+sub and/or subgrid for the fitting run.


    Parameters
    ----------
    use_sd : boolean (default=True)
        If True, create source density dependent noise models (determined by
        finding matches to datamodel.astfile with SD info)

    nsubs : int (default=1)
        number of subgrids used for the physics model

    nprocs : int (default=1)
        Number of parallel processes to use
        (currently only implemented for subgrids)

    choose_sd_sub : list of two strings (default=None)
        If this is set, the fitting will just be for this combo of SD+sub,
        rather than all of them.  Overrides use_sd.
        format of the list: ['#','#']

    choose_subgrid : int (default=None)
        If this is set, the fitting with just be for this subgrid index.
        If nsubs=1, this is ignored.

    pdf2d_param_list : list of strings or None
        If set, do 2D PDFs of these parameters.  If None, don't make 2D PDFs.

    resume : boolean (default=False)
        choose whether to resume existing run or start over

    """

    # before doing ANYTHING, force datamodel to re-import (otherwise, any
    # changes within this python session will not be loaded!)
    importlib.reload(datamodel)
    # check input parameters
    verify_params.verify_input_format(datamodel)

    # keep track of time
    start_time = time.clock()

    # --------------------
    # make lists of file names
    # --------------------

    file_dict = create_filenames.create_filenames(
        use_sd=use_sd,
        nsubs=nsubs,
        choose_sd_sub=choose_sd_sub,
        choose_subgrid=choose_subgrid,
    )

    # input files
    photometry_files = file_dict["photometry_files"]
    # modelsedgrid_files = file_dict["modelsedgrid_files"]
    modelsedgrid_trim_files = file_dict["modelsedgrid_trim_files"]
    # noise_files = file_dict["noise_files"]
    noise_trim_files = file_dict["noise_trim_files"]

    # output files
    stats_files = file_dict["stats_files"]
    pdf_files = file_dict["pdf_files"]
    pdf2d_files = file_dict["pdf2d_files"]
    if pdf2d_param_list is None:
        pdf2d_files = [None for i in range(len(pdf2d_files))]
    lnp_files = file_dict["lnp_files"]

    # total number of files
    n_files = len(photometry_files)

    # other potentially useful info
    sd_sub_info = file_dict["sd_sub_info"]
    # gridsub_info = file_dict['gridsub_info']

    # if using subgrids, make the grid dictionary file:
    # File where the ranges and number of unique values for the grid
    # will be stored (this can take a while to calculate)

    if nsubs > 1:

        gridpickle_files = file_dict["gridpickle_files"]

        for i in range(len(gridpickle_files)):
            if not os.path.isfile(gridpickle_files[i]):

                # list of corresponding SED grids and noise models

                # - with SD+sub: get file list for ALL subgrids at current SD+sub
                if use_sd or (choose_sd_sub is not None):
                    temp = create_filenames.create_filenames(
                        nsubs=nsubs, choose_sd_sub=sd_sub_info[i], choose_subgrid=None
                    )
                    modelsedgrid_trim_list = temp["modelsedgrid_trim_files"]
                    noise_trim_list = temp["noise_trim_files"]

                # - no SD info: get file list for ALL subgrids
                else:
                    temp = create_filenames.create_filenames(
                        use_sd=False, nsubs=nsubs, choose_subgrid=None
                    )
                    modelsedgrid_trim_list = temp["modelsedgrid_trim_files"]
                    noise_trim_list = temp["noise_trim_files"]

                # create the grid info dictionary
                print("creating grid_info_dict for " + gridpickle_files[i])
                grid_info_dict = subgridding_tools.reduce_grid_info(
                    modelsedgrid_trim_list, noise_trim_list, nprocs=nprocs
                )
                # save it
                with open(gridpickle_files[i], "wb") as p:
                    pickle.dump(grid_info_dict, p)
                print("wrote grid_info_dict to " + gridpickle_files[i])

    # --------------------
    # do the fitting!
    # --------------------

    # set up function inputs

    if nsubs == 1:

        input_list = [
            (
                photometry_files[i],
                modelsedgrid_trim_files[i],
                noise_trim_files[i],
                stats_files[i],
                pdf_files[i],
                pdf2d_files[i],
                pdf2d_param_list,
                lnp_files[i],
                None,
                resume,
            )
            for i in range(n_files)
        ]

    if nsubs > 1:

        input_list = [
            (
                photometry_files[i],
                modelsedgrid_trim_files[i],
                noise_trim_files[i],
                stats_files[i],
                pdf_files[i],
                pdf2d_files[i],
                pdf2d_param_list,
                lnp_files[i],
                gridpickle_files[i],
                resume,
            )
            for i in range(n_files)
        ]

    # run the fitting (via parallel wrapper)

    parallel_wrapper(fit_submodel, input_list, nprocs=nprocs)

    # see how long it took!
    new_time = time.clock()
    print("time to fit: ", (new_time - start_time) / 60.0, " min")


def fit_submodel(
    photometry_file,
    modelsedgrid_file,
    noise_file,
    stats_file,
    pdf_file,
    pdf2d_file,
    pdf2d_param_list,
    lnp_file,
    grid_info_file=None,
    resume=False,
):
    """
    Code to run the SED fitting

    Parameters
    ----------
    photometry_file : string
        path+name of the photometry file

    modelsedgrid_file : string
        path+name of the physics model grid file

    noise_file : string
        path+name of the noise model file

    stats_file : string
        path+name of the file to contain stats output

    pdf_file : string
        path+name of the file to contain 1D PDF output

    pdf2d_file : string
        path+name of the file to contain 2D PDF output

    pdf2d_param_list: list of strings or None
        parameters for which to make 2D PDFs (or None)

    lnp_file : string
        path+name of the file to contain log likelihood output

    grid_info_file : string (default=None)
        path+name for pickle file that contains dictionary with subgrid
        min/max/n_unique (required for a run with subgrids)

    resume : boolean (default=False)
        choose whether to resume existing run or start over


    Returns
    -------
    noisefile : string
        name of the created noise file

    """

    # read in the photometry catalog
    obsdata = datamodel.get_obscat(photometry_file, datamodel.filters)

    # check if it's a subgrid run by looking in the file name
    if "gridsub" in modelsedgrid_file:
        subgrid_run = True
        print("loading grid_info_dict from " + grid_info_file)
        with open(grid_info_file, "rb") as p:
            grid_info_dict = pickle.loads(p.read())
    else:
        subgrid_run = False

    # load the SED grid and noise model
    modelsedgrid = FileSEDGrid(modelsedgrid_file)
    noisemodel_vals = noisemodel.get_noisemodelcat(noise_file)

    if subgrid_run:
        fit.summary_table_memory(
            obsdata,
            noisemodel_vals,
            modelsedgrid,
            resume=resume,
            threshold=-10.0,
            save_every_npts=100,
            lnp_npts=500,
            stats_outname=stats_file,
            pdf1d_outname=pdf_file,
            pdf2d_outname=pdf2d_file,
            pdf2d_param_list=pdf2d_param_list,
            grid_info_dict=grid_info_dict,
            lnp_outname=lnp_file,
            do_not_normalize=True,
            surveyname=datamodel.surveyname,
        )
        print("Done fitting on grid " + modelsedgrid_file)

    else:

        fit.summary_table_memory(
            obsdata,
            noisemodel_vals,
            modelsedgrid,
            resume=resume,
            threshold=-10.0,
            save_every_npts=100,
            lnp_npts=500,
            stats_outname=stats_file,
            pdf1d_outname=pdf_file,
            pdf2d_outname=pdf2d_file,
            pdf2d_param_list=pdf2d_param_list,
            lnp_outname=lnp_file,
            surveyname=datamodel.surveyname,
        )
        print("Done fitting on grid " + modelsedgrid_file)


if __name__ == "__main__":  # pragma: no cover
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--use_sd",
        help="create source density dependent noise models",
        action="store_true",
    )
    parser.add_argument(
        "--nsubs",
        type=int,
        default=1,
        help="""Number of subgrids that the physics model
                        was split into""",
    )
    parser.add_argument(
        "--nprocs",
        type=int,
        default=1,
        help="""Number of processes to use to process
                        the subgrids""",
    )
    parser.add_argument(
        "--choose_sd_sub",
        nargs=2,
        default=None,
        help="Fit just this combo of SD+sub. Format: ['#','#']",
    )
    parser.add_argument(
        "--choose_subgrid", type=int, default=None, help="Fit just this subgrid number"
    )
    parser.add_argument(
        "--pdf2d_param_list",
        type=str,
        nargs="+",
        default=['Av', 'Rv', 'f_A', 'M_ini', 'logA', 'Z', 'distance'],
        help="If set, do 2D PDFs of these parameters. If None, don't make 2D PDFs."
    )
    parser.add_argument(
        "-r", "--resume", help="resume a fitting run", action="store_true"
    )

    args = parser.parse_args()

    if 'None' in args.pdf2d_param_list:
        args.pdf2d_param_list = None

    run_fitting(
        use_sd=args.use_sd,
        nsubs=args.nsubs,
        nprocs=args.nprocs,
        choose_sd_sub=args.choose_sd_sub,
        choose_subgrid=args.choose_subgrid,
        pdf2d_param_list=args.pdf2d_param_list,
        resume=args.resume,
    )
