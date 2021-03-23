import argparse

# BEAST imports
from beast.tools import beast_settings, subgridding_tools, merge_beast_stats
from beast.tools.run import create_filenames


def merge_files(beast_settings_info, use_sd=True, nsubs=1, partial=False):
    """
    Merge all of the results from the assorted fitting sub-files (divided by
    source density, subgrids, or both).

    If fitting is in progress but you want to check results of completed stars,
    set partial=True.  This is only relevant when using subgrids.

    Parameters
    ----------
    beast_settings_info : string or beast.tools.beast_settings.beast_settings instance
        if string: file name with beast settings
        if class: beast.tools.beast_settings.beast_settings instance

    use_sd : boolean (default=True)
        set to True if the fitting used source density bins

    nsubs : int (default=1)
        number of subgrids used for the physics model

    partial : boolean (default=False)
        If True, the output merged files will only have stars that have been
        run across all subgrids.  If stars have only been fit in some subgrids
        and not others, they will be discarded in the "partial" output files.
        Currently only implemented for 1D PDFs and stats (not lnP) files.

    """

    # if there's no SD and no subgridding, running this is unnecessary
    if (not use_sd) and (nsubs == 1):
        print("No merging necessary")
        return

    # process beast settings info
    if isinstance(beast_settings_info, str):
        settings = beast_settings.beast_settings(beast_settings_info)
    elif isinstance(beast_settings_info, beast_settings.beast_settings):
        settings = beast_settings_info
    else:
        raise TypeError(
            "beast_settings_info must be string or beast.tools.beast_settings.beast_settings instance"
        )

    # get file name lists (to check if they exist and/or need to be resumed)
    file_dict = create_filenames.create_filenames(settings, use_sd=use_sd, nsubs=nsubs)

    # - input files
    # photometry_files = file_dict['photometry_files']
    # modelsedgrid_files = file_dict['modelsedgrid_files']
    # noise_files = file_dict['noise_files']

    # - output files
    stats_files = file_dict["stats_files"]
    pdf_files = file_dict["pdf_files"]
    lnp_files = file_dict["lnp_files"]

    # - other useful info
    sd_sub_info = file_dict["sd_sub_info"]
    # gridsub_info = file_dict['gridsub_info']
    # the unique sets of gridsub
    unique_sd_sub = [x for i, x in enumerate(sd_sub_info) if i == sd_sub_info.index(x)]

    # --------------------
    # no subgrids
    # --------------------

    if nsubs == 1:

        out_filebase = "{0}/{0}".format(settings.project)
        reorder_tags = ["bin{0}_sub{1}".format(x[0], x[1]) for x in unique_sd_sub]
        merge_beast_stats.merge_stats_files(
            stats_files, out_filebase, reorder_tag_list=reorder_tags
        )

    # --------------------
    # use subgrids
    # --------------------

    if nsubs > 1:

        # runs were split by source density
        if use_sd:

            # lists to save the merged file names
            merged_pdf_files = []
            merged_stats_files = []
            merged_lnp_files = []

            for sd_sub in unique_sd_sub:

                # indices with the current sd_sub
                ind = [j for j, x in enumerate(sd_sub_info) if x == sd_sub]

                # merge the subgrid files for that SD+sub
                out_filebase = "{0}/bin{1}_sub{2}/{0}_bin{1}_sub{2}".format(
                    settings.project, sd_sub[0], sd_sub[1]
                )
                if partial:
                    out_filebase += "_partial"

                # - 1D PDFs and stats
                (
                    merged_pdf1d_fname,
                    merged_stats_fname,
                ) = subgridding_tools.merge_pdf1d_stats(
                    [pdf_files[j] for j in ind],
                    [stats_files[j] for j in ind],
                    re_run=False,
                    output_fname_base=out_filebase,
                    partial=partial,
                )

                merged_pdf_files.append(merged_pdf1d_fname)
                merged_stats_files.append(merged_stats_fname)

                # - lnP files
                if not partial:
                    merged_lnp_fname = subgridding_tools.merge_lnp(
                        [lnp_files[j] for j in ind],
                        re_run=False,
                        output_fname_base=out_filebase,
                        threshold=-10,
                    )
                    merged_lnp_files.append(merged_lnp_fname)

            # merge the merged stats files
            out_filebase = "{0}/{0}".format(settings.project)
            reorder_tags = ["bin{0}_sub{1}".format(x[0], x[1]) for x in unique_sd_sub]
            merge_beast_stats.merge_stats_files(
                merged_stats_files, out_filebase, reorder_tag_list=reorder_tags
            )

        # runs weren't split by source density
        else:

            out_filebase = "{0}/{0}".format(settings.project)

            # - 1D PDFs and stats
            subgridding_tools.merge_pdf1d_stats(
                pdf_files, stats_files, output_fname_base=out_filebase, partial=partial,
            )

            # - lnP files
            if not partial:
                subgridding_tools.merge_lnp(
                    lnp_files,
                    re_run=False,
                    output_fname_base=out_filebase,
                    threshold=-10,
                )


if __name__ == "__main__":  # pragma: no cover
    # commandline parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "beast_settings_file", type=str, help="file name with beast settings",
    )
    parser.add_argument(
        "--use_sd",
        type=int,
        default=1,
        help="set to True (1) if the fitting used source density bins",
    )
    parser.add_argument(
        "--nsubs",
        type=int,
        default=1,
        help="number of subgrids used for the physics model",
    )
    parser.add_argument(
        "--partial",
        action="store_true",
        help="Set this if fitting is partially complete and you want to merge stars that are done (relevant to subgrids only)",
    )

    args = parser.parse_args()

    merge_files(
        beast_settings_info=args.beast_settings_file,
        use_sd=bool(args.use_sd),
        nsubs=args.nsubs,
        partial=args.partial,
    )

    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()
