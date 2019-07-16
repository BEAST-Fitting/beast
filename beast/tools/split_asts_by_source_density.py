from __future__ import print_function

import numpy as np
from astropy.table import Table

def split_asts(ast_file, sd_map_file, bin_width=1.):
    """
    Split the ASTs into sub-files for each source density bin.

    Parameters
    ----------
    ast_file : string
        Name of the file that contains the AST results (mag_in and mag_out)

    sd_map_file : string
        Name of the fits file that contains the source densities

    bin_width : float
        Width of source density bin in star/arcsec

    """

    # read in the AST file
    ast_data = Table.read(ast_file)

    # read in the source density file
    sd_data = Table.read(sd_map_file)

    # define SD bins
    sd_bins = np.arange(np.floor(np.min(sd_data['value'])),
                            np.ceil(np.max(sd_data['value']))+1,
                            bin_width,
                            dtype=int)


    # go through each source density bin
    for i in range(len(sd_bins)-1):

        print('getting ASTs in SD bin '+str(sd_bins[i])+'-'+str(sd_bins[i+1]) )

        # list to hold all of the indices
        ast_all_ind = []

        # indices for this bin
        sd_ind = np.where((sd_data['value'] >= sd_bins[i]) &
                              (sd_data['value'] < sd_bins[i+1]))[0]

        # for each index, grab the ASTs within that RA/Dec box
        for j,ind in enumerate(sd_ind):

            ast_ind = np.where((ast_data['RA_J2000'] > sd_data['min_ra'][ind]) &
                                   (ast_data['RA_J2000'] < sd_data['max_ra'][ind]) &
                                   (ast_data['DEC_J2000'] > sd_data['min_dec'][ind]) & 
                                   (ast_data['DEC_J2000'] < sd_data['max_dec'][ind]) )[0]

            # append indices to master list
            if len(ast_ind) > 0:
                ast_all_ind += ast_ind.tolist()


        # make a new table out of the indices
        print('   found '+str(len(ast_all_ind))+' ASTs')
        if len(ast_all_ind) > 0:
            print('   writing new table')
            ast_table = ast_data[ast_all_ind]
            new_filename = ast_file.replace('.fits','_SD'+str(sd_bins[i])+'-'+str(sd_bins[i+1])+'.fits')
            ast_table.write(new_filename, overwrite=True)
        

