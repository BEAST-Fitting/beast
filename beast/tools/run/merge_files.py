# system imports
from __future__ import (absolute_import, division, print_function)


# BEAST imports
from beast.tools import (verify_params,
                         subgridding_tools,
                         merge_beast_stats)
from beast.tools.run import create_filenames


import datamodel
import importlib

#import pdb

def merge_files(use_sd=True, nsubs=1):
    """
    Merge all of the results from the assorted fitting sub-files (divided by 
    source density, subgrids, or both).


    Parameters
    ----------
    use_sd : boolean (default=True)
        If True, create source density dependent noise models (determined by
        finding matches to datamodel.astfile with SD info)

    nsubs : int (default=1)
        number of subgrids used for the physics model

    """

    # if there's no SD and no subgridding, running this is unnecessary
    if (use_sd == False) and (nsubs == 1):
        print('No merging necessary')
        return
    
    # before doing ANYTHING, force datamodel to re-import (otherwise, any
    # changes within this python session will not be loaded!)
    importlib.reload(datamodel)
    # check input parameters
    verify_params.verify_input_format(datamodel)


    # get file name lists (to check if they exist and/or need to be resumed)
    file_dict = create_filenames.create_filenames(use_sd=use_sd, nsubs=nsubs)    

    # - input files
    #photometry_files = file_dict['photometry_files']
    #modelsedgrid_files = file_dict['modelsedgrid_files']
    #noise_files = file_dict['noise_files']

    # - output files
    stats_files = file_dict['stats_files']
    pdf_files = file_dict['pdf_files']
    #lnp_files = file_dict['lnp_files']

    # - other useful info
    sd_sub_info = file_dict['sd_sub_info']
    #gridsub_info = file_dict['gridsub_info']
    # the unique sets of gridsub
    unique_sd_sub = [x for i, x in enumerate(sd_sub_info) if i == sd_sub_info.index(x)]


    # --------------------
    # no subgrids
    # --------------------
    
    if nsubs == 1:

        out_filebase = '{0}/{0}'.format(datamodel.project)
        reorder_tags = ['sd{0}_sub{1}'.format(x[0],x[1]) for x in unique_sd_sub]
        merge_beast_stats.merge_stats_files(stats_files, out_filebase,
                                            reorder_tag_list=reorder_tags)



    # --------------------
    # use subgrids
    # --------------------

    if nsubs > 1:
    
        # runs were split by source density
        if use_sd == True:

            # lists to save the merged file names
            merged_pdf_files = []
            merged_stats_files = []

            for i,sd_sub in enumerate(unique_sd_sub):

                # indices with the current sd_sub
                ind = [j for j,x in enumerate(sd_sub_info) if x == sd_sub]

                # merge the subgrid files for that SD+sub
                out_filebase = '{0}/SD{1}_sub{2}/{0}_SD{1}_sub{2}'.format(
                    datamodel.project, sd_sub[0], sd_sub[1])

                merged_pdf1d_fname, merged_stats_fname = \
                    subgridding_tools.merge_pdf1d_stats([pdf_files[j] for j in ind],
                                                        [stats_files[j] for j in ind],
                                                        re_run=False,
                                                        output_fname_base=out_filebase)

                merged_pdf_files.append(merged_pdf1d_fname)
                merged_stats_files.append(merged_stats_fname)
 
            # merge the merged stats files
            out_filebase = '{0}/{0}'.format(datamodel.project)
            reorder_tags = ['sd{0}_sub{1}'.format(x[0],x[1]) for x in unique_sd_sub]
            merge_beast_stats.merge_stats_files(merged_stats_files, out_filebase,
                                                reorder_tag_list=reorder_tags)


                

        # runs weren't split by source density
        if use_sd == False:

            out_filebase = '{0}/{0}'.format(datamodel.project)
 
            subgridding_tools.merge_pdf1d_stats(pdf_files,
                                                stats_files,
                                                output_fname_base=out_filebase)

