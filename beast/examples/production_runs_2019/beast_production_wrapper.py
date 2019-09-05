import numpy as np
import glob
import os
import types

from beast.tools.run import (create_physicsmodel,
                             make_ast_inputs,
                             create_obsmodel,
                             run_fitting,
                             merge_files,
                             create_filenames)

from beast.plotting import plot_mag_hist
from beast.tools import (create_background_density_map,
                             split_ast_input_file,
                             subdivide_obscat_by_source_density,
                             cut_catalogs,
                             split_asts_by_source_density,
                             setup_batch_beast_trim,
                             setup_batch_beast_fit)

import importlib
importlib.reload(create_physicsmodel)
importlib.reload(make_ast_inputs)
importlib.reload(create_obsmodel)
importlib.reload(run_fitting)
importlib.reload(setup_batch_beast_fit)
importlib.reload(merge_files)
importlib.reload(create_filenames)

from astropy.table import Table
from astropy.coordinates import Angle
from astropy import units as u



def beast_production_wrapper():
    """
    This does all of the steps for a full production run, and can be used as
    a wrapper to automatically do most steps for multiple fields.
    * make datamodel.py file
    * make source density map
    * make background density map
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

    # the list of fields
    field_names = ['15275_IC1613']

    # distance moduli and velocities
    # http://adsabs.harvard.edu/abs/2013AJ....146...86T
    dist_mod = [24.36]
    velocity = [-236]

    # the path+file for a reference image
    im_path = ['../beast_dwarfs/images/15275_IC1613_F555W_drz.fits.gz']
    ref_filter = ['F555W']

    # choose a filter to use for removing artifacts
    # (remove catalog sources with filter_FLAG > 99)
    flag_filter = ['F555W']

    # number of fields
    n_field = len(field_names)

    # Need to know what the correspondence is between filter names in the
    # catalog and the BEAST filter names.
    #
    # These will be used to automatically determine the filters present in
    # each GST file and fill in the datamodel.py file.  The order doesn't
    # matter, as long as the order in one list matches the order in the other
    # list.
    #
    gst_filter_names = ['F275W','F336W','F390M','F555W',
                            'F814W','F110W','F160W']
    beast_filter_names = ['HST_WFC3_F275W','HST_WFC3_F336W','HST_WFC3_F390M','HST_WFC3_F555W',
                              'HST_WFC3_F814W', 'HST_WFC3_F110W','HST_WFC3_F160W']



    for b in range(n_field):
    #for b in [0]:

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
        im_file = im_path[b]



        # region file with catalog stars
        #make_region_file(gst_file)

        # -----------------
        # 0. make datamodel file
        # -----------------

        # need to do this first, because otherwise any old version that exists
        # will be imported, and changes made here won't get imported again

        print('')
        print('creating datamodel file')
        print('')

        create_datamodel(gst_file, ast_file, gst_filter_names, beast_filter_names,
                             dist_mod[b], velocity[b], ref_image=im_file)

        # load in datamodel to get number of subgrids
        import datamodel; importlib.reload(datamodel)

        # -----------------
        # 1a. make magnitude histograms
        # -----------------

        print('')
        print('making magnitude histograms')
        print('')

        #if not os.path.isfile('./data/'+field_names[b]+'.gst_maghist.pdf'):
        peak_mags = plot_mag_hist.plot_mag_hist(gst_file, stars_per_bin=70, max_bins=75)
        #test = plot_mag_hist.plot_mag_hist(ast_file, stars_per_bin=200, max_bins=30)


        # -----------------
        # 1b. make a source density map
        # -----------------

        print('')
        print('making source density map')
        print('')


        # not currently doing background density bins
        #use_bg_info = True
        use_bg_info = False
        if use_bg_info:
            background_args = types.SimpleNamespace(subcommand='background',
                                                    catfile=gst_file, pixsize=5, npix=None,
                                                    reference=im_file,
                                                    mask_radius=10, ann_width=20,
                                                    cat_filter=[ref_filter[b],'90'])
            create_background_density_map.main_make_map(background_args)

        # but we are doing source density bins!
        if not os.path.isfile(gst_file.replace('.fits','_source_den_image.fits')):
            # - pixel size of 10 arcsec
            # - use ref_filter[b] between vega mags of 17 and peak_mags[ref_filter[b]]-0.5
            sourceden_args = types.SimpleNamespace(subcommand='sourceden',
                                                   catfile=gst_file, pixsize=8, npix=None,
                                                   mag_name=ref_filter[b]+'_VEGA',
                                                   mag_cut=[15,peak_mags[ref_filter[b]]-0.5])
            create_background_density_map.main_make_map(sourceden_args)

        # new file name with the source density column
        gst_file_sd = gst_file.replace('.fits', '_with_sourceden.fits')



        # -----------------
        # 2. make physics model
        # -----------------

        print('')
        print('making physics model')
        print('')

        # see which subgrid files already exist
        gs_str = ''
        if datamodel.n_subgrid > 1:
            gs_str = 'sub*'

        spec_files = glob.glob('./' + field_names[b] + '_beast/' + field_names[b]
                                + '_beast_spec_w_priors.grid'+gs_str+'.hd5')

        # only make the physics model they don't already exist
        if len(spec_files) < datamodel.n_subgrid:
            create_physicsmodel.create_physicsmodel(nprocs=1, nsubs=datamodel.n_subgrid)

        # list of SED files
        model_grid_files = sorted(glob.glob('./' + field_names[b] + '_beast/' + field_names[b]
                                            + '_beast_seds.grid'+gs_str+'.hd5'))


        # -----------------
        # 3. make ASTs
        # -----------------

        # only create an AST input list if the ASTs don't already exist
        ast_input_file = './' + field_names[b] + '_beast/' + field_names[b] + '_beast_inputAST.txt'

        if not os.path.isfile(ast_file):
            if not os.path.isfile(ast_input_file):
                print('')
                print('creating artificial stars')
                print('')
                make_ast_inputs.make_ast_inputs(flux_bin_method=True)

            split_ast_input_file.split_asts(
                field_names[b] + '_beast', ast_input_file, 2000)

            print('\n**** go run ASTs for '+field_names[b]+'! ****\n')
            continue

        # -----------------
        # 4/5. edit photometry/AST catalogs
        # -----------------

        # remove sources that are
        # - in regions without full imaging coverage,
        # - flagged in flag_filter

        print('')
        print('editing photometry/AST catalogs')
        print('')

        # - photometry
        gst_file_cut = gst_file.replace('.fits', '_with_sourceden_cut.fits')
        cut_catalogs.cut_catalogs(gst_file_sd, gst_file_cut, partial_overlap=True,
                                      flagged=True, flag_filter=flag_filter[b],
                                      region_file=True)

        # - ASTs
        ast_file_cut = ast_file.replace('.fits', '_cut.fits')
        cut_catalogs.cut_catalogs(ast_file, ast_file_cut, partial_overlap=True,
                                      flagged=True, flag_filter=flag_filter[b],
                                      region_file=True)
        #test = plot_mag_hist.plot_mag_hist(ast_file_cut, stars_per_bin=200, max_bins=30)

        # edit the datamodel.py file to have the correct photometry file name
        # (AST file name is already automatically the cut version)
        create_datamodel(gst_file_cut, ast_file_cut, gst_filter_names, beast_filter_names,
                             dist_mod[b], velocity[b], ref_image=im_file)


        # -----------------
        # 6. split observations by source density
        # -----------------

        print('')
        print('splitting observations by source density')
        print('')

        # - photometry

        if len(glob.glob(gst_file_cut.replace('.fits','*sub*fits') )) == 0:

            # a smaller value for Ns_file will mean more individual files/runs,
            # but each run will take a shorter amount of time

            subdivide_obscat_by_source_density.split_obs_by_source_density(gst_file_cut, bin_width=1,
                                                                           sort_col=ref_filter[b]+'_RATE',
                                                                           Ns_file=6250)


        # - ASTs

        # check if any files exist already
        ast_files = sorted(glob.glob(ast_file_cut.replace('.fits','_SD_*.fits')))

        if len(ast_files) == 0:
            split_asts_by_source_density.split_asts(ast_file_cut,
                                                    gst_file.replace('.fits','_sourceden_map.hd5'))


        # -- at this point, we can run the code to create lists of filenames
        file_dict = create_filenames.create_filenames(use_sd=True,
                                                      nsubs=datamodel.n_subgrid)

        # figure out how many files there are
        sd_sub_info = file_dict['sd_sub_info']
        # - number of SD bins
        temp = set([i[0] for i in sd_sub_info])
        print('** total SD bins: '+str(len(temp)))
        # - the unique sets of SD+sub
        unique_sd_sub = [x for i, x in enumerate(sd_sub_info) if i == sd_sub_info.index(x)]
        print('** total SD subfiles: '+str(len(unique_sd_sub)) )


        # -----------------
        # 7. make noise models
        # -----------------

        print('')
        print('making noise models')
        print('')

        # create the noise model (this code will check if it exists)
        create_obsmodel.create_obsmodel(use_sd=True,
                                        nsubs=datamodel.n_subgrid,
                                        nprocs=1)


        # -----------------
        # 8. make script to trim models
        # -----------------

        print('')
        print('setting up script to trim models')
        print('')


        # save any at-queue commands
        at_list = []

        # iterate through each model grid
        for i in range(datamodel.n_subgrid):

            # gst list
            temp = file_dict['photometry_files']
            gst_input_list = [x for i, x in enumerate(temp) if i == temp.index(x)]

            # create corresponding files for each of those
            ast_input_list = []
            noise_files = []
            trim_prefix = []

            for j in range(len(gst_input_list)):
                # get the sd/sub number
                curr_sd = unique_sd_sub[j][0]
                curr_sub = unique_sd_sub[j][1]
                subfolder = 'SD{0}_sub{1}'.format(curr_sd, curr_sub)


                # create file names
                ast_input_list.append(ast_file_cut.replace('.fits','_SD'+curr_sd+'.fits'))
                if datamodel.n_subgrid > 1:
                    noise_files.append('./' + field_names[b] + '_beast/' + field_names[b]
                                       + '_beast_noisemodel_SD'+curr_sd+'.gridsub'+str(i)+'.hd5')
                    trim_prefix.append('./'+field_names[b]+'_beast/'+subfolder+'/'+field_names[b]+'_beast_'
                                       +subfolder+ '_gridsub'+str(i))
                if datamodel.n_subgrid == 1:
                    noise_files.append('./' + field_names[b] + '_beast/' + field_names[b]
                                       + '_beast_noisemodel_SD'+curr_sd+'.hd5')
                    trim_prefix.append('./'+field_names[b]+'_beast/'+field_names[b]+'_beast_'
                                       +subfolder)



            # check if the trimmed grids exist before moving on
            if datamodel.n_subgrid > 1:
                trim_files = sorted(glob.glob('./' + field_names[b] + '_beast/SD*_sub*/'+field_names[b] +
                                       '_beast_*_gridsub'+str(i)+'_sed_trim.grid.hd5') )
            if datamodel.n_subgrid == 1:
                trim_files = sorted(glob.glob('./' + field_names[b] + '_beast/'+field_names[b] +
                                       '_beast_*_sub*_sed_trim.grid.hd5') )


            if len(trim_files) < len(gst_input_list):

                job_path = './' + field_names[b] + '_beast/trim_batch_jobs/'
                if datamodel.n_subgrid > 1:
                    file_prefix='BEAST_gridsub'+str(i)
                if datamodel.n_subgrid == 1:
                    file_prefix='BEAST'

                # generate trimming at-queue commands
                setup_batch_beast_trim.generic_batch_trim(model_grid_files[i],
                                                          noise_files,
                                                          gst_input_list,
                                                          ast_input_list,
                                                          trim_prefix,
                                                          job_path=job_path,
                                                          file_prefix=file_prefix,
                                                          num_subtrim=1, nice=19,
                                                          prefix='source activate b13')


                at_list.append('at -f '+job_path+file_prefix+'_batch_trim.joblist now')


        if len(at_list) > 0:
            print('\n**** go run trimming code for '+field_names[b]+'! ****')
            print('Here are the command(s) to run:')
            for cmd in at_list:
                print(cmd)
            return
        else:
            print('all files are trimmed for '+field_names[b])



        # -----------------
        # 9. make script to fit models
        # -----------------

        print('')
        print('setting up script to fit models')
        print('')

        fit_run_info = setup_batch_beast_fit.setup_batch_beast_fit(
                num_percore=1, nice=19,
                overwrite_logfile=False,
                prefix='source activate b13',
                use_sd=True, nsubs=datamodel.n_subgrid, nprocs=1)


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
        # 10. merge stats files from each fit
        # -----------------

        print('')
        print('merging stats files')
        print('')

        merge_files.merge_files(use_sd=True, nsubs=datamodel.n_subgrid)





def create_datamodel(gst_file, ast_file, gst_filter_label, beast_filter_label,
                         dist_mod, velocity, ref_image='None'):
    """
    Create a datamodel.py file for the given field.  This will open the file to
    determine the filters present - the `*_filter_label` inputs are references
    to properly interpret the file's information.

    Parameters
    ----------
    gst_file : string
        the path+name of the GST file

    ast_file : string
        the path+name of the AST file

    gst_filter_label : list of strings
        Labels used to represent each filter in the photometry catalog

    beast_filter_label : list of strings
        The corresponding full labels used by the BEAST

    dist_mod : float
        distance modulus of the source (mag)

    velocity : float
        velocity to use for the source (km/s)

    ref_image : string (default='None')
        path+name of image to use as reference for ASTs

    Returns
    -------
    nothing

    """

    # read in the catalog
    cat = Table.read(gst_file)

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
            new_file.write('project = "' + gst_file.split('/')[-1].split('.')[0] + '_beast"\n')
        # obsfile
        elif datamodel_lines[i][0:10] == 'obsfile = ':
            new_file.write('obsfile = "' + gst_file + '"\n')
        # AST file name
        elif datamodel_lines[i][0:10] == 'astfile = ':
            new_file.write('astfile = "'+ ast_file +'"\n')
        # BEAST filter names
        elif datamodel_lines[i][0:10] == 'filters = ':
            new_file.write("filters = ['" + "','".join(filter_list_long) + "'] \n")
        # catalog filter names
        elif datamodel_lines[i][0:14] == 'basefilters = ':
            new_file.write("basefilters = ['" + "','".join(filter_list_base) + "'] \n")
        # distance modulus
        elif datamodel_lines[i][0:12] == 'distances = ':
            new_file.write('distances = ['+str(dist_mod)+'] \n')
        # velocity
        elif datamodel_lines[i][0:11] == 'velocity = ':
            new_file.write('velocity = '+str(velocity)+' * units.km / units.s \n')
        # AST stuff
        #elif datamodel_lines[i][0:27] == 'ast_source_density_table = ':
        #    new_file.write('ast_source_density_table = "' + gst_file.replace('.fits','_sourcedens_map.hd5')+'" \n')
        elif datamodel_lines[i][0:22] == 'ast_reference_image = ':
            new_file.write('ast_reference_image = "' + ref_image+'" \n')
        # none of those -> write line as-is
        else:
            new_file.write(datamodel_lines[i])

    new_file.close()




def make_region_file(gst_file):
    """
    Make a region file out of the catalog file
    """

    with open(gst_file.replace('.fits','.reg'),'w') as ds9_file:
        ds9_file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        ds9_file.write('fk5\n')

        cat = Table.read(gst_file)
        for i in range(len(cat)):
            #if cat['F555W_RATE'][i] != 0:
            #if cat[flag_filter[b]+'_FLAG'][i] < 99:
            if cat[ref_filter[b]+'_VEGA'][i] < 26.65:
                ds9_file.write('circle(' +
                                   Angle(cat['RA'][i], u.deg).to_string(unit=u.hour, sep=':')
                                   + ',' +
                                   Angle(cat['DEC'][i], u.deg).to_string(unit=u.deg, sep=':')
                                   + ',0.1")\n' )
            else:
                ds9_file.write('circle(' +
                                   Angle(cat['RA'][i], u.deg).to_string(unit=u.hour, sep=':')
                                   + ',' +
                                   Angle(cat['DEC'][i], u.deg).to_string(unit=u.deg, sep=':')
                                   + ',0.1") # color=magenta \n' )





if __name__ == '__main__':

    beast_production_wrapper()
