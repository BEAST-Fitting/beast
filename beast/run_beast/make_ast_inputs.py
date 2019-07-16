#!/usr/bin/env python
"""
Script to run the BEAST on the PHAT-like data.
Assumes that the datamodel.py file exists in the same directory as this script.
  And it must be called datamodel.py 
"""

# system imports
from __future__ import (absolute_import, division, print_function)
import argparse
import numpy as np

# BEAST imports
from beast.observationmodel.ast.make_ast_input_list import (pick_models,
                                                            pick_models_toothpick_style)
from beast.observationmodel.ast.make_ast_xy_list import pick_positions
from beast.tools import verify_params

import datamodel
import importlib

#import pdb


def make_ast_inputs(flux_bin_method=True):
    """
    Make the list of artificial stars to be run through the photometry pipeline

    Parameters
    ----------
    flux_bin_method : boolean (default=True)
        If True, use the flux bin method to select SEDs.  If False, randomly
        select SEDs from the model grid.

    """

    # before doing ANYTHING, force datamodel to re-import (otherwise, any
    # changes within this python session will not be loaded!)
    importlib.reload(datamodel)
    # check input parameters
    verify_params.verify_input_format(datamodel)
    

    # construct magnitude cuts
    
    mag_cuts = datamodel.ast_maglimit
    obsdata = datamodel.get_obscat(datamodel.obsfile, datamodel.filters)

    if len(mag_cuts) == 1:
        tmp_cuts = mag_cuts
        min_mags = np.zeros(len(datamodel.filters))
        for k, filtername in enumerate(obsdata.filters):
            sfiltername = obsdata.data.resolve_alias(filtername)
            sfiltername = sfiltername.replace('rate','vega')
            sfiltername = sfiltername.replace('RATE','VEGA')
            keep, = np.where(obsdata[sfiltername] < 99.)
            min_mags[k] = np.percentile(obsdata[keep][sfiltername],90.)

        # max. mags from the gst observation cat.
        mag_cuts = min_mags + tmp_cuts


        
    # --------------------
    # select SEDs
    # --------------------

    Nrealize = datamodel.ast_realization_per_model
    Nfilters = datamodel.ast_bands_above_maglimit

    # file names for stars and corresponding SED parameters
    outfile = './' + datamodel.project + '/' + datamodel.project + '_inputAST.txt'
    outfile_params = './' + datamodel.project + '/' + datamodel.project + '_ASTparams.fits'


    if flux_bin_method == True:
        
        N_fluxes = datamodel.ast_n_flux_bins
        min_N_per_flux = datamodel.ast_n_per_flux_bin
        bins_outfile = './' + datamodel.project + '/' + datamodel.project + '_ASTfluxbins.txt'
        
        chosen_seds = pick_models_toothpick_style(modelsedgrid_filename, datamodel.filters,
                            mag_cuts, Nfilters, N_fluxes, min_N_per_flux,
                            outfile=outfile, outfile_params=outfile_params,
                            bins_outfile=bins_outfile)


    if flux_bin_method == False:

        N_models = datamodel.ast_models_selected_per_age
        
        chosen_seds = pick_models(modelsedgrid_filename, datamodel.filters,
                                  mag_cuts, Nfilter=Nfilters, N_stars=N_models, Nrealize=Nrealize,
                                  outfile=outfile, outfile_params=outfile_params)


    # --------------------
    # assign positions
    # --------------------

        
    if datamodel.ast_with_positions == True:
        separation = datamodel.ast_pixel_distribution
        filename = datamodel.project + '/' + datamodel.project + '_inputAST.txt'

        if datamodel.ast_reference_image is not None:
            # With reference image, use one of these options
            if datamodel.ast_source_density_table is not None:
                pick_positions_from_map(obsdata,
                                        chosen_seds,
                                        datamodel.ast_source_density_table,
                                        datamodel.ast_N_bins,
                                        datamodel.ast_realization_per_model,
                                        outfile=filename,
                                        refimage=datamodel.ast_reference_image,
                                        refimage_hdu=0,
                                        Nrealize=1,
                                        set_coord_boundary=datamodel.ast_coord_boundary)

            elif datamodel.ast_background_table is not None:
                pick_positions_from_map(obsdata,
                                        chosen_seds,
                                        datamodel.ast_background_table,
                                        datamodel.ast_N_bins,
                                        datamodel.ast_realization_per_model,
                                        outfile=filename,
                                        refimage=datamodel.ast_reference_image,
                                        refimage_hdu=0,
                                        Nrealize=1,
                                        set_coord_boundary=datamodel.ast_coord_boundary)
            else:
                pick_positions(obsdata, filename, separation,
                               	refimage=datamodel.ast_reference_image)

        else:
            # Without reference image, we can only use this function
            if (datamodel.ast_source_density_table is None and
                datamodel.ast_background_table is None):
                pick_positions(obsdata, filename, separation)
            else:
                print("To use ast_source_density_table or ast_background_table, ast_reference_image must be specified.")





if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--flux_bin_method",
                        help="Choose SEDs using flux bin method",
                        action="store_true")

    args = parser.parse_args()


    make_ast_inputs(flux_bin_method=args.flux_bin_method)

    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()

