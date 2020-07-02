# system imports
import os
import numpy as np
import glob

# BEAST imports
from beast.tools import beast_settings
from beast.tools.run.helper_functions import get_modelsubgridfiles


def create_filenames(
    beast_settings_info, use_sd=True, nsubs=1, choose_sd_sub=None, choose_subgrid=None,
):
    """
    Helper function to make all of the filenames.  SED grid and noise model
    are trimmed versions.

    Parameters
    ----------
    beast_settings_info : string or beast.tools.beast_settings.beast_settings instance
        if string: file name with beast settings
        if class: beast.tools.beast_settings.beast_settings instance

    use_sd : boolean (default=True)
        If True, create source density dependent noise models (determined by
        finding matches to settings.astfile with SD info)

    nsubs : int (default=1)
        number of subgrids used for the physics model

    choose_sd_sub : list of two strings (default=None)
        If this is set, the fitting will just be for this combo of SD+sub,
        rather than all of them.  Overrides use_sd.
        format of the list: ['#','#']

    choose_subgrid : int (default=None)
        If this is set, the fitting with just be for this subgrid index.
        If nsubs=1, this is ignored.

    Returns
    -------
    dictionary with the lists of filenames, plus the corresponding SD+sub and
    gridsub values for easy referencing

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

    # input files
    photometry_files = []
    modelsedgrid_files = []
    modelsedgrid_trim_files = []
    noise_files = []
    noise_trim_files = []

    # output files
    stats_files = []
    pdf_files = []
    pdf2d_files = []
    lnp_files = []

    # other potentially useful things
    sd_sub_info = []
    gridsub_info = []

    # ** no subgrids **

    if nsubs == 1:

        # -- SD+sub specified
        if choose_sd_sub is not None:

            photometry_files.append(
                settings.obsfile.replace(
                    ".fits",
                    "_bin{0}_sub{1}.fits".format(choose_sd_sub[0], choose_sd_sub[1]),
                )
            )
            modelsedgrid_files.append("{0}/{0}_seds.grid.hd5".format(settings.project))
            modelsedgrid_trim_files.append(
                "{0}/{0}_bin{1}_sub{2}_seds_trim.grid.hd5".format(
                    settings.project, choose_sd_sub[0], choose_sd_sub[1]
                )
            )
            noise_files.append(
                "{0}/{0}_noisemodel_bin{1}.grid.hd5".format(
                    settings.project, choose_sd_sub[0]
                )
            )
            noise_trim_files.append(
                "{0}/{0}_bin{1}_sub{2}_noisemodel_trim.grid.hd5".format(
                    settings.project, choose_sd_sub[0], choose_sd_sub[1]
                )
            )

            stats_files.append(
                "{0}/{0}_bin{1}_sub{2}_stats.fits".format(
                    settings.project, choose_sd_sub[0], choose_sd_sub[1]
                )
            )
            pdf_files.append(
                "{0}/{0}_bin{1}_sub{2}_pdf1d.fits".format(
                    settings.project, choose_sd_sub[0], choose_sd_sub[1]
                )
            )
            pdf2d_files.append(
                "{0}/{0}_bin{1}_sub{2}_pdf2d.fits".format(
                    settings.project, choose_sd_sub[0], choose_sd_sub[1]
                )
            )
            lnp_files.append(
                "{0}/{0}_bin{1}_sub{2}_lnp.hd5".format(
                    settings.project, choose_sd_sub[0], choose_sd_sub[1]
                )
            )

            sd_sub_info.append([choose_sd_sub[0], choose_sd_sub[1]])

        # -- using source density info
        elif use_sd is True:

            photometry_files = sorted(
                glob.glob(settings.obsfile.replace(".fits", "_bin*_sub*.fits"))
            )

            for phot_file in photometry_files:
                # get the sd/sub number
                dpos = phot_file.rfind("_bin")
                spos = phot_file.rfind("sub")
                ppos = phot_file.rfind(".")
                curr_sd = phot_file[dpos + 4 : spos - 1]
                curr_sub = phot_file[spos + 3 : ppos]

                # construct other file names
                modelsedgrid_files.append(
                    "{0}/{0}_seds.grid.hd5".format(settings.project)
                )
                modelsedgrid_trim_files.append(
                    "{0}/{0}_bin{1}_sub{2}_seds_trim.grid.hd5".format(
                        settings.project, curr_sd, curr_sub
                    )
                )
                noise_files.append(
                    "{0}/{0}_noisemodel_bin{1}.grid.hd5".format(
                        settings.project, curr_sd
                    )
                )
                noise_trim_files.append(
                    "{0}/{0}_bin{1}_sub{2}_noisemodel_trim.grid.hd5".format(
                        settings.project, curr_sd, curr_sub
                    )
                )

                stats_files.append(
                    "{0}/{0}_bin{1}_sub{2}_stats.fits".format(
                        settings.project, curr_sd, curr_sub
                    )
                )
                pdf_files.append(
                    "{0}/{0}_bin{1}_sub{2}_pdf1d.fits".format(
                        settings.project, curr_sd, curr_sub
                    )
                )
                pdf2d_files.append(
                    "{0}/{0}_bin{1}_sub{2}_pdf2d.fits".format(
                        settings.project, curr_sd, curr_sub
                    )
                )
                lnp_files.append(
                    "{0}/{0}_bin{1}_sub{2}_lnp.hd5".format(
                        settings.project, curr_sd, curr_sub
                    )
                )

                sd_sub_info.append([curr_sd, curr_sub])

        # -- no source density splitting
        else:

            photometry_files.append(settings.obsfile)
            modelsedgrid_files.append("{0}/{0}_seds.grid.hd5".format(settings.project))
            modelsedgrid_trim_files.append(
                "{0}/{0}_seds_trim.grid.hd5".format(settings.project)
            )
            noise_files.append("{0}/{0}_noisemodel.grid.hd5".format(settings.project))
            noise_trim_files.append(
                "{0}/{0}_noisemodel_trim.grid.hd5".format(settings.project)
            )

            stats_files.append("{0}/{0}_stats.fits".format(settings.project))
            pdf_files.append("{0}/{0}_pdf1d.fits".format(settings.project))
            pdf2d_files.append("{0}/{0}_pdf2d.fits".format(settings.project))
            lnp_files.append("{0}/{0}_lnp.hd5".format(settings.project))

    # ** with subgrids **

    # subgrids require a pickle file with grid info
    gridpickle_files = []

    if nsubs > 1:

        # start with getting the model grid files (note these aren't trimmed ones)
        outdir = os.path.join(".", settings.project)
        subgrid_names_file = os.path.join(outdir, "subgrid_fnames.txt")
        temp = get_modelsubgridfiles(subgrid_names_file)
        # use that to get the number of subgrids and make a list of them
        gridsub_list = np.arange(len(temp))
        # or a subset if set
        if choose_subgrid is not None:
            gridsub_list = [choose_subgrid]

        # -- SD+sub specified
        if choose_sd_sub is not None:

            for gridsub in gridsub_list:

                photometry_files.append(
                    settings.obsfile.replace(
                        ".fits",
                        "_bin{0}_sub{1}.fits".format(
                            choose_sd_sub[0], choose_sd_sub[1]
                        ),
                    )
                )

                modelsedgrid_files.append(
                    "{0}/{0}_seds.gridsub{1}.hd5".format(settings.project, gridsub)
                )
                modelsedgrid_trim_files.append(
                    "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_seds_trim.grid.hd5".format(
                        settings.project, choose_sd_sub[0], choose_sd_sub[1], gridsub
                    )
                )
                noise_files.append(
                    "{0}/{0}_noisemodel_bin{1}.gridsub{2}.hd5".format(
                        settings.project, choose_sd_sub[0], gridsub
                    )
                )
                noise_trim_files.append(
                    "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_noisemodel_trim.grid.hd5".format(
                        settings.project, choose_sd_sub[0], choose_sd_sub[1], gridsub
                    )
                )

                stats_files.append(
                    "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_stats.fits".format(
                        settings.project, choose_sd_sub[0], choose_sd_sub[1], gridsub
                    )
                )
                pdf_files.append(
                    "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_pdf1d.fits".format(
                        settings.project, choose_sd_sub[0], choose_sd_sub[1], gridsub
                    )
                )
                pdf2d_files.append(
                    "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_pdf2d.fits".format(
                        settings.project, choose_sd_sub[0], choose_sd_sub[1], gridsub
                    )
                )
                lnp_files.append(
                    "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_lnp.hd5".format(
                        settings.project, choose_sd_sub[0], choose_sd_sub[1], gridsub
                    )
                )

                gridpickle_files.append(
                    "{0}/bin{1}_sub{2}/grid_info_dict.pkl".format(
                        settings.project, choose_sd_sub[0], choose_sd_sub[1]
                    )
                )

                sd_sub_info.append([choose_sd_sub[0], choose_sd_sub[1]])
                gridsub_info.append(gridsub)

        # -- using source density info
        elif use_sd is True:

            phot_file_list = sorted(
                glob.glob(settings.obsfile.replace(".fits", "_bin*_sub*.fits"))
            )

            for phot_file in phot_file_list:
                # get the sd/sub number
                dpos = phot_file.rfind("_bin")
                spos = phot_file.rfind("sub")
                ppos = phot_file.rfind(".")
                curr_sd = phot_file[dpos + 4 : spos - 1]
                curr_sub = phot_file[spos + 3 : ppos]

                # construct other file names
                for gridsub in gridsub_list:
                    photometry_files.append(phot_file)
                    modelsedgrid_files.append(
                        "{0}/{0}_seds.gridsub{1}.hd5".format(settings.project, gridsub)
                    )
                    modelsedgrid_trim_files.append(
                        "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_seds_trim.grid.hd5".format(
                            settings.project, curr_sd, curr_sub, gridsub
                        )
                    )
                    noise_files.append(
                        "{0}/{0}_noisemodel_bin{1}.gridsub{2}.hd5".format(
                            settings.project, curr_sd, gridsub
                        )
                    )
                    noise_trim_files.append(
                        "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_noisemodel_trim.grid.hd5".format(
                            settings.project, curr_sd, curr_sub, gridsub
                        )
                    )

                    stats_files.append(
                        "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_stats.fits".format(
                            settings.project, curr_sd, curr_sub, gridsub
                        )
                    )
                    pdf_files.append(
                        "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_pdf1d.fits".format(
                            settings.project, curr_sd, curr_sub, gridsub
                        )
                    )
                    pdf2d_files.append(
                        "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_pdf2d.fits".format(
                            settings.project, curr_sd, curr_sub, gridsub
                        )
                    )
                    lnp_files.append(
                        "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}_gridsub{3}_lnp.hd5".format(
                            settings.project, curr_sd, curr_sub, gridsub
                        )
                    )

                    gridpickle_files.append(
                        "{0}/bin{1}_sub{2}/grid_info_dict.pkl".format(
                            settings.project, curr_sd, curr_sub
                        )
                    )

                    sd_sub_info.append([curr_sd, curr_sub])
                    gridsub_info.append(gridsub)

        # -- no source density splitting
        else:

            for gridsub in gridsub_list:
                photometry_files.append(settings.obsfile)
                modelsedgrid_files.append(
                    "{0}/{0}_seds.gridsub{1}.hd5".format(settings.project, gridsub)
                )
                modelsedgrid_trim_files.append(
                    "{0}/{0}_gridsub{1}_seds_trim.grid.hd5".format(
                        settings.project, gridsub
                    )
                )
                noise_files.append(
                    "{0}/{0}_noisemodel.gridsub{1}.hd5".format(
                        settings.project, gridsub
                    )
                )
                noise_trim_files.append(
                    "{0}/{0}_gridsub{1}_noisemodel_trim.grid.hd5".format(
                        settings.project, gridsub
                    )
                )

                stats_files.append(
                    "{0}/{0}_gridsub{1}_stats.fits".format(settings.project, gridsub)
                )
                pdf_files.append(
                    "{0}/{0}_gridsub{1}_pdf1d.fits".format(settings.project, gridsub)
                )
                pdf2d_files.append(
                    "{0}/{0}_gridsub{1}_pdf2d.fits".format(settings.project, gridsub)
                )
                lnp_files.append(
                    "{0}/{0}_gridsub{1}_lnp.hd5".format(settings.project, gridsub)
                )

                gridpickle_files.append(
                    "{0}/grid_info_dict.pkl".format(settings.project)
                )

                gridsub_info.append(gridsub)

    # double check that all file lists are the same length
    n_file_list = [
        len(x)
        for x in [
            photometry_files,
            modelsedgrid_files,
            modelsedgrid_trim_files,
            noise_files,
            noise_trim_files,
            stats_files,
            pdf_files,
            pdf2d_files,
            lnp_files,
        ]
    ]
    if len(np.unique(n_file_list)) > 1:
        print("file list lengths don't match!")
        return None

    return {
        "photometry_files": photometry_files,
        "modelsedgrid_files": modelsedgrid_files,
        "modelsedgrid_trim_files": modelsedgrid_trim_files,
        "noise_files": noise_files,
        "noise_trim_files": noise_trim_files,
        "stats_files": stats_files,
        "pdf_files": pdf_files,
        "pdf2d_files": pdf2d_files,
        "lnp_files": lnp_files,
        "gridpickle_files": gridpickle_files,
        "sd_sub_info": sd_sub_info,
        "gridsub_info": gridsub_info,
    }
