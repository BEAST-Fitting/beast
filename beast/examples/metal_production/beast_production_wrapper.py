import numpy as np
import glob
import subprocess
import sys
import os

from run_beast_production import run_beast_production
from beast.tools import create_source_density_map
from beast.tools import create_background_density_map
from beast.tools import subdivide_obscat_by_source_density
from beast.tools import merge_beast_stats
from beast.tools import setup_batch_beast_trim
from beast.tools import setup_batch_beast_fit
from beast.tools import reorder_beast_results_spatial
from beast.tools import condense_beast_results_spatial

from astropy.table import Table


def beast_production_wrapper():
    """
    This does all of the steps for a full production run, and can be used as
    a wrapper to automatically do most steps for multiple fields.
    * make datamodel.py file
    * make source density map
    * split catalog by source density
    * make physics model (SED grid)
    * make input list for ASTs
    * make noise model
    * generate batch script to trim models
    * generate batch script to fit models
    * merge stats files back together
    * spatially reorder the results

    Places for user to manually do things:
    * editing code before use
        - datamodel_template.py: setting up the file with desired parameters
        - here: list the catalog filter names with the corresponding BEAST names
        - here: choose settings (pixel size, filter, mag range) for the source density map
        - here: choose settings (pixel size, reference image) for the background map
        - here: choose settings (filter, number per file) for dividing catalog by source density
        - here: choose settings (# files, nice level) for the trimming/fitting batch scripts
    * process the ASTs, as described in BEAST documentation
    * run the trimming scripts
    * run the fitting scripts

    BEWARE: When running the trimming/fitting scripts, ensure that the correct
    datamodel.py file is in use.  Since it gets updated every time this code is
    run, you may be unexpectedly be using one from another field.
    """

    # list of field names:
    #    These should be the base of the GST file names, and will be used
    #    to create the project name
    field_names = ['14675_LMC-5665ne-12232']
    n_field = len(field_names)

    
    # Need to know what the correspondence is between filter names in the
    # catalog and the BEAST filter names.
    #
    # These will be used to automatically determine the filters present in
    # each GST file and fill in the datamodel.py file.  The order doesn't
    # matter, as long as the order in one list matches the order in the other
    # list.
    #
    gst_filter_names = ['F225W','F275W','F336W','F475W','F814W','F110W','F160W']
    beast_filter_names = ['HST_WFC3_F225W','HST_WFC3_F275W','HST_WFC3_F336W',
                              'HST_ACS_WFC_F475W','HST_ACS_WFC_F814W',
                              'HST_WFC3_F110W','HST_WFC3_F160W']
    
    

    for b in range(n_field):

        print('********')
        print('field ' + field_names[b])
        print('********')

        
        # -----------------
        # data file names
        # -----------------

        # paths for the data/AST files
        gst_file = './data/' + field_names[b]+'.gst.fits'
        ast_file = './data/' + field_names[b]+'.gst.fake.fits'
        # path for the reference image (if using for the background map)
        im_file = './data/'+field_names[b]+'_F475W.fits.gz'
        
        # -----------------
        # make datamodel file
        # -----------------

        # need to do this first, because otherwise any old version that exists
        # will be imported, and changes made here won't get imported again

        print('')
        print('creating datamodel file')
        print('')

        create_datamodel(gst_file, gst_filter_names, beast_filter_names)


        # -----------------
        # make a source density map
        # -----------------

        print('')
        print('making source density map')
        print('')
        
        if not os.path.isfile(gst_file.replace('.fits','_source_den_image.fits')):
            # - pixel size of 10 arcsec
            # - use F475W between vega mags of 17 and 27
            create_source_density_map.make_source_dens_map(gst_file, pix_size=10,
                                                            mag_name='F475W_VEGA', mag_cut=[17,27])

        # new file name with the source density column
        gst_file_new = gst_file.replace('.fits', '_with_sourceden.fits')

        
        # -----------------
        # make a background map
        # -----------------

        print('')
        print('making background map')
        print('')
        
        if not os.path.isfile(gst_file_new.replace('.fits','_F475W_bg_map.hd5')):
            # - pixel dimensions: 15x15
            create_background_density_map.create_background_density_map(gst_file_new, npix=15,
                                                      reference=im_file, filebase='F475W')

        # new file name with the background column
        #gst_file_new = gst_file_new.replace('.fits', '_with_bg.fits')


        # -----------------
        # split observations by source density
        # -----------------

        print('')
        print('splitting observations by source density')
        print('')

        if len(glob.glob(gst_file_new.replace('.fits','*sub*fits') )) == 0:

            # a smaller value for Ns_file will mean more individual files/runs,
            # but each run will take a shorter amount of time
            
            subdivide_obscat_by_source_density.split_obs_by_source_density(gst_file_new, bin_width=1,
                                                                            sort_col='F475W_RATE', Ns_file=6250)

        # figure out how many files there are
        tot_sub_files = len(glob.glob(gst_file_new.replace('.fits','*sub*fits') ))
        print('** total subfiles: '+str(tot_sub_files))

        
        
        # -----------------
        # make physics model
        # -----------------

        print('')
        print('making physics model')
        print('')

        # the file name for the model grid
        physics_model_file = './' + field_names[b] + '_beast/' + field_names[b] + '_beast_seds.grid.hd5'

        # only make the physics model if it doesn't already exist
        if not os.path.isfile(physics_model_file):
            run_beast_production(gst_file, physicsmodel=True)


        # -----------------
        # make ASTs
        # -----------------

        # only create an AST input list if the ASTs don't already exist
        ast_input_file = './' + field_names[b] + '_beast/' + field_names[b] + '_inputAST.txt'
            
        if not os.path.isfile(ast_file):
            if not os.path.isfile(ast_input_file):
                print('')
                print('creating artificial stars')
                print('')
                run_beast_production(gst_file, ast=True)
                
            print('\n**** go run ASTs for '+field_names[b]+'! ****\n')
            continue
       
            
        # -----------------
        # make noise model
        # -----------------
      
        print('')
        print('making noise model')
        print('')

        # eventually, this may be divided into files based on source density,
        # but for now it'll be one giant file

        # the file name for the model grid
        noise_model_file = './' + field_names[b] + '_beast/' + field_names[b] + '_beast_noisemodel.hd5'

        if not os.path.isfile(noise_model_file):
            run_beast_production(gst_file, observationmodel=True)#, source_density='0', sub_source_density='0')


        # -----------------
        # make script to trim models
        # -----------------

        print('')
        print('setting up script to trim models')
        print('')
        
        # check if the trimmed grids exist before moving on
        trim_files = glob.glob('./' + field_names[b] + '_beast/' + field_names[b] + '_beast_*_sed_trim.grid.hd5')
        
        if len(trim_files) < tot_sub_files:

            # choose num_subtrim to be the number of CPUs you'll run it on
            # (if it's more than 1, you'll need to split the joblist file manually)

            setup_batch_beast_trim.setup_batch_beast_trim(field_names[b] + '_beast',
                                                            gst_file, ast_file,
                                                            num_subtrim=1, nice=19)

            print('\n**** go run trimming code for '+field_names[b]+'! ****')
            print('Here is the command to run:')
            print('at -f ./'+field_names[b]+'_beast/trim_batch_jobs/beast_batch_trim.joblist now \n')
            continue
        
        else:
            print('all files are trimmed for '+field_names[b])


        # -----------------
        # make script to fit models
        # -----------------

        print('')
        print('setting up script to fit models')
        print('')
        
        fit_run_info = setup_batch_beast_fit.setup_batch_beast_fit(field_names[b] + '_beast',
                                                                       gst_file,
                                                                       num_percore=1, nice=19,
                                                                       overwrite_logfile=False)

        # check if the fits exist before moving on
        tot_remaining = len(fit_run_info['done']) - np.sum(fit_run_info['done'])
        if tot_remaining > 0:
            print('\n**** go run fitting code for '+field_names[b]+'! ****')
            print('Here are the '+str(len(fit_run_info['files_to_run']))+' commands to run:')
            for job_file in fit_run_info['files_to_run']:
                print('at -f ./'+job_file+' now')
            continue
        else:
            print('all fits are complete for '+field_names[b])


        
        # -----------------
        # merge stats files from each fit
        # -----------------

        print('')
        print('merging stats files')
        print('')

        stats_filebase = './' + field_names[b] + '_beast/' + field_names[b] + '_beast'
        
        merge_beast_stats.merge_stats_files(glob.glob(stats_filebase+'*sub*_stats.fits'), stats_filebase)
        


        # -----------------
        # reorganize results into spatial regions
        # -----------------

        print('')
        print('doing spatial reorganizing')
        print('')

        region_filebase = './' + field_names[b] + '_beast/' + field_names[b] + '_beast_sd'
        output_filebase = './' + field_names[b] + '_beast/spatial/' + field_names[b]

        reorder_beast_results_spatial.reorder_beast_results_spatial(stats_filename=stats_filebase + '_stats.fits',
                                                                        region_filebase=region_filebase,
                                                                        output_filebase=output_filebase)

        condense_beast_results_spatial.condense_files(filedir='./' + field_names[b] + '_beast/spatial/')
        
        

    
def create_datamodel(basename, gst_filter_label, beast_filter_label):
    """
    Create a datamodel.py file for the given field.  This will open the file to
    determine the filters present - the `*_filter_label` inputs are references
    to properly interpret the file's information.

    Parameters
    ----------
    basename : string
        the path+name of the GST file

    gst_filter_label : list of strings
        Labels used to represent each filter in the photometry catalog

    beast_filter_label : list of strings
        The corresponding full labels used by the BEAST

    Returns
    -------
    nothing

    """

    # read in the catalog
    cat = Table.read(basename)

    # get the list of filters
    filter_list_base = []
    filter_list_long = []
    for f in range(len(gst_filter_label)):
        filt_exist = [gst_filter_label[f] in c for c in cat.colnames]
        if np.sum(filt_exist) > 0:
            filter_list_base.append(gst_filter_label[f])
            filter_list_long.append(beast_filter_label[f])


    # read in the template datamodel file
    orig_file = open('datamodel_template.py','r')
    datamodel_lines = np.array( orig_file.readlines() )
    orig_file.close()

    # write out an edited datamodel
    new_file = open('datamodel.py','w')

    for i in range(len(datamodel_lines)):

        # replace project name with the field ID
        if datamodel_lines[i][0:10] == 'project = ':
            new_file.write('project = "' + basename.split('/')[-1].split('.')[0] + '_beast"\n')
        # obsfile (just a placeholder - will be updated in run_beast_production.py)
        elif datamodel_lines[i][0:10] == 'obsfile = ':
            new_file.write('obsfile = "' + basename + '"\n')
        # AST file name
        elif datamodel_lines[i][0:10] == 'astfile = ':
            new_file.write('astfile = "'+ basename.replace('.gst.fits','.gst.fake.fits') +'"\n')
        # BEAST filter names
        elif datamodel_lines[i][0:10] == 'filters = ':
            new_file.write("filters = ['" + "','".join(filter_list_long) + "'] \n")
        # catalog filter names
        elif datamodel_lines[i][0:14] == 'basefilters = ':
            new_file.write("basefilters = ['" + "','".join(filter_list_base) + "'] \n")
        else:
            new_file.write(datamodel_lines[i])
        
    new_file.close()

    

    
            
        
if __name__ == '__main__':

    beast_production_wrapper()
