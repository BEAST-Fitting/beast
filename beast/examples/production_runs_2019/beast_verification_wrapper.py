import numpy as np
import glob
import os

from beast.tools.run import (
    run_fitting,
    merge_files,
    create_filenames,
    make_trim_scripts,
)

from beast.tools import (
    simulate_obs,
    setup_batch_beast_fit,
)
from beast.plotting import (
    plot_param_recovery,
    plot_param_err,
    plot_triangle,
)

import importlib

importlib.reload(make_trim_scripts)
importlib.reload(run_fitting)
importlib.reload(setup_batch_beast_fit)
importlib.reload(merge_files)
importlib.reload(create_filenames)

from astropy.table import Table, vstack


def beast_verification_wrapper():
    """
    This wrapper does the processing for BEAST verification

    Parameter recovery
    * create simulated data for a given model grid + noise model
    * generate batch script to trim models
    * generate batch script to fit models
    * merge stats files back together

    Places for user to manually do things:
    * editing code before use
        - datamodel_template.py: setting up the file with desired parameters
        - here: list the catalog filter names with the corresponding BEAST names
        - here: number of simulated stars to generate
        - here: choose settings (# files, nice level) for the trimming/fitting batch scripts
    * process the ASTs, as described in BEAST documentation
    * run the trimming scripts
    * run the fitting scripts

    BEWARE: When running the trimming/fitting scripts, ensure that the correct
    datamodel.py file is in use.  Since it gets updated every time this code is
    run, you may be unexpectedly be using one from another field.
    """

    # the list of fields
    field_names = ["15275_IC1613"]

    # distance moduli and velocities
    # http://adsabs.harvard.edu/abs/2013AJ....146...86T
    dist_mod = [24.36]
    velocity = [-236]

    # the path+file for a reference image
    im_path = ["../beast_dwarfs/images/15275_IC1613_F555W_drz.fits.gz"]
    ref_filter = ["F555W"]

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
    gst_filter_names = ["F275W", "F336W", "F390M", "F555W", "F814W", "F110W", "F160W"]
    beast_filter_names = [
        "HST_WFC3_F275W",
        "HST_WFC3_F336W",
        "HST_WFC3_F390M",
        "HST_WFC3_F555W",
        "HST_WFC3_F814W",
        "HST_WFC3_F110W",
        "HST_WFC3_F160W",
    ]

    for b in range(n_field):
        # for b in [0]:

        print("********")
        print("field " + field_names[b])
        print("********")

        # -----------------
        # 0. get original file names
        # -----------------

        print('')
        print('retrieving original file names')
        print('')

        # paths for the data/AST files
        gst_file_orig = './data/' + field_names[b]+'.gst_with_sourceden_cut.fits'
        ast_file_orig = './data/' + field_names[b]+'.gst.fake_cut.fits'
        # path for the reference image (if using for the background map)
        im_file = im_path[b]

        create_datamodel(
            gst_file_orig,
            ast_file_orig,
            gst_filter_names,
            beast_filter_names,
            dist_mod[b],
            velocity[b],
            ref_image=im_file,
            proj_type='beast',
        )

        # load in datamodel to get number of subgrids
        import datamodel; importlib.reload(datamodel)

        # grab relevant file names
        file_dict = create_filenames.create_filenames(
            use_sd=True,
            nsubs=datamodel.n_subgrid,
        )
        modelsedgrid_files = file_dict['modelsedgrid_files']
        noise_files = file_dict['noise_files']
        sd_sub_info = file_dict["sd_sub_info"]
        tot_files = len(sd_sub_info)


        # -----------------
        # 1. create simulated data
        # -----------------

        print('')
        print('simulating data')
        print('')

        gst_file = gst_file_orig.replace('.gst','.sim.gst')
        gst_subfile_form = gst_file.replace(
            '.fits','_bin*_sub0.fits')

        # loop through all files
        # only grab files and simulate data when we're at sd_sub = [N,0]
        for i in range(tot_files):

            # find matches to sd_sub = [i,0]
            inds_to_use = [ind for ind in range(tot_files)
                if sd_sub_info[ind]==[str(i),'0']]

            # if there are matches, use those corresponding files
            if len(inds_to_use) > 0:

                output_catalog = gst_subfile_form.replace('_bin*','_bin'+str(i))

                if not os.path.isfile(output_catalog):

                    print('generating simulated observations for bin='+str(i))

                    grid_sublist = [modelsedgrid_files[x] for x in inds_to_use]
                    noise_sublist = [noise_files[x] for x in inds_to_use]

                    simulate_obs.simulate_obs(
                        grid_sublist,
                        noise_sublist,
                        output_catalog,
                        nsim=5000,
                        compl_filter=ref_filter[b],
                    )

                else:
                    print('simulated observations already exist for bin='+str(i))


        # combine them all into one catalog
        if not os.path.isfile(gst_file):
            table_list = []
            for cat_file in glob.glob(gst_subfile_form):
                table_list.append(Table.read(cat_file))
            vstack(table_list).write(gst_file, overwrite=True)

        # -----------------
        # 2. make new datamodel file
        # -----------------

        print('')
        print('creating datamodel file')
        print('')

        create_datamodel(
            gst_file,
            ast_file_orig,
            gst_filter_names,
            beast_filter_names,
            dist_mod[b],
            velocity[b],
            ref_image=im_file,
            proj_type='sim',
        )

        # load in datamodel again
        importlib.reload(datamodel)


        # -----------------
        # 3. make symbolic links to model grids and noise models
        # -----------------

        # make new directory
        if not os.path.isdir('./'+datamodel.project):
            os.mkdir('./'+datamodel.project)

        # symlink the physics/noise models
        orig_phys = list(set(modelsedgrid_files))
        for grid in orig_phys:
            source = os.path.abspath(grid)
            dest = os.path.abspath(grid.replace('_beast','_sim'))
            if not os.path.islink(dest):
                os.symlink(source, dest)
        orig_noise = list(set(noise_files))
        for grid in orig_noise:
            source = os.path.abspath(grid)
            dest = os.path.abspath(grid.replace('_beast','_sim'))
            if not os.path.islink(dest):
                os.symlink(source, dest)

        # -----------------
        # 4. make script to trim models
        # -----------------

        print("")
        print("setting up script to trim models")
        print("")

        job_file_list = make_trim_scripts.make_trim_scripts(
            num_subtrim=1, prefix='source activate b13'
        )

        if len(job_file_list) > 0:
            print('\n**** go run trimming code for '+field_names[b]+'! ****')
            print('Here are the command(s) to run:')
            for job in job_file_list:
                print('at -f '+job+' now')
            return
        else:
            print('all files are trimmed for '+field_names[b])



        # -----------------
        # 5. make script to fit models
        # -----------------

        print("")
        print("setting up script to fit models")
        print("")

        fit_run_info = setup_batch_beast_fit.setup_batch_beast_fit(
            num_percore=1,
            nice=19,
            overwrite_logfile=False,
            prefix="source activate b13",
            use_sd=True,
            nsubs=datamodel.n_subgrid,
            nprocs=1,
        )

        # check if the fits exist before moving on
        tot_remaining = len(fit_run_info["done"]) - np.sum(fit_run_info["done"])
        if tot_remaining > 0:
            print("\n**** go run fitting code for " + field_names[b] + "! ****")
            print(
                "Here are the "
                + str(len(fit_run_info["files_to_run"]))
                + " commands to run:"
            )
            for job_file in fit_run_info["files_to_run"]:
                print("at -f ./" + job_file + " now")
            continue
        else:
            print("all fits are complete for " + field_names[b])


        # -----------------
        # 6. plots
        # -----------------

        print('')
        print('making plots')
        print('')

        # grab relevant file names
        file_dict = create_filenames.create_filenames(
            use_sd=True,
            nsubs=datamodel.n_subgrid,
        )

        plot_param_recovery.plot_param_recovery(
            file_dict['photometry_files'],
            file_dict['stats_files'],
            field_names[b]+'_param_recovery.pdf',
            max_nbins=20,
        )

        for stats_file in file_dict['stats_files']:
            plot_param_err.plot(stats_file, n_bins=10)
            plot_triangle.plot_triangle(stats_file)


def create_datamodel(
    gst_file,
    ast_file,
    gst_filter_label,
    beast_filter_label,
    dist_mod,
    velocity,
    ref_image="None",
    proj_type="beast",
):
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

    proj_type : string (default='beast')
        type of project: BEAST run ('beast') or simulation ('sim')

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
    orig_file = open("datamodel_template.py", "r")
    datamodel_lines = np.array(orig_file.readlines())
    orig_file.close()

    # write out an edited datamodel
    new_file = open("datamodel.py", "w")

    for i in range(len(datamodel_lines)):

        # replace project name with the field ID
        if datamodel_lines[i][0:10] == 'project = ':
            new_file.write(
                'project = "' + gst_file.split('/')[-1].split('.')[0]+ '_'+proj_type+'"\n'
            )
        # obsfile
        elif datamodel_lines[i][0:10] == "obsfile = ":
            new_file.write('obsfile = "' + gst_file + '"\n')
        # AST file name
        elif datamodel_lines[i][0:10] == "astfile = ":
            new_file.write('astfile = "' + ast_file + '"\n')
        # BEAST filter names
        elif datamodel_lines[i][0:10] == "filters = ":
            new_file.write("filters = ['" + "','".join(filter_list_long) + "'] \n")
        # catalog filter names
        elif datamodel_lines[i][0:14] == "basefilters = ":
            new_file.write("basefilters = ['" + "','".join(filter_list_base) + "'] \n")
        # distance modulus
        elif datamodel_lines[i][0:12] == "distances = ":
            new_file.write("distances = [" + str(dist_mod) + "] \n")
        # velocity
        elif datamodel_lines[i][0:11] == "velocity = ":
            new_file.write("velocity = " + str(velocity) + " * units.km / units.s \n")
        # AST stuff
        # elif datamodel_lines[i][0:27] == 'ast_source_density_table = ':
        #    new_file.write('ast_source_density_table = "' + gst_file.replace('.fits','_sourcedens_map.hd5')+'" \n')
        elif datamodel_lines[i][0:22] == "ast_reference_image = ":
            new_file.write('ast_reference_image = "' + ref_image + '" \n')
        # none of those -> write line as-is
        else:
            new_file.write(datamodel_lines[i])

    new_file.close()



if __name__ == "__main__":

    beast_verification_wrapper()
