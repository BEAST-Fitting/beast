"""
Example of running the BEAST on M31 PHAT data.
Karl G. - 5 Nov 2014

------ previous comments by Morgan
Everything I need to generate a grid, make fake data and fit

I use the pipeline package I wrote in order to clean the syntax and allow more
flexibilities. In particular it will simplifies the managment of intermediate
results or broken jobs.

> python run [--models] [-?]

    --models         generates the models only (needed once)
    -?, --help       display this message

Make models
-----------
:func:`make_models` is equivalent to using individual tasks as follow:

project, noisefile, grid = project | t_project_dir
                                   | t_isochrones(**iso_kwargs)
                                   | t_spectra(**spec_kwargs)
                                   | t_seds(filters, **seds_kwargs)
                                   | t_gen_noise_model(astfile, **noise_kwargs)

Running the fit
---------------
:func:`run_fit`  is equivalent to using individual tasks as follow:

project, stat, obs, sedgrid = project | t_project_dir
                                      | t_get_obscat(**obscat_kwargs)
                                      | t_fit(g, noise, **fit_kwargs)
                                      | t_summary_table(g, **stat_kwargs)
"""

# system imports
import sys

# BEAST imports
from pipeline import run_fit, make_models
import datamodel
import noisemodel 
from merge_phat_asts import merge_phat_asts
from beast.core import prior_weights

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if ('-?' in sys.argv[1]) or ('--help' in sys.argv[1]):
            print(__doc__)
        elif '-models' in sys.argv[1]:
            make_models()
        elif '-priors' in sys.argv[1]:
            modelsedgrid = '{project:s}/{project:s}_seds.grid.hd5'.format(project=datamodel.project)
            modelsedgridwpriors = '{project:s}/{project:s}_seds_w_priors.grid.hd5'.format(project=datamodel.project)
            prior_weights.add_priors_sedgrid(modelsedgrid, modelsedgridwpriors, datamodel.filters)
        elif '-noise' in sys.argv[1]:
            print 'generating noise model from ASTs'
            # merge the single camera ASTs into a single file
            merge_phat_asts(datamodel.uvastfile,datamodel.optastfile,datamodel.irastfile,datamodel.astfile)

            # get the modesedgrid on which to generate the noisemodel
            from beast.core.grid import FileSEDGrid
            modelsedgrid = FileSEDGrid('{project:s}/{project:s}_seds_w_priors.grid.hd5'.format(project=datamodel.project))
                        
            # generate the AST noise model
            noisemodel.make_toothpick_noise_model(datamodel.noisefile, datamodel.astfile, modelsedgrid, datamodel.absflux_a_matrix)
        else:
            print sys.argv[1] + ' option is not supported'

    else:
        # define the file in which to store the grid of model SED
        modelsedgrid = '{project:s}/{project:s}_seds_w_priors.grid.hd5'.format(project=datamodel.project)

        # read in the the AST noise model
        noisemodel_vals = noisemodel.get_noisemodelcat(datamodel.noisefile)

        run_fit(datamodel.project, modelsedgrid, noise=noisemodel_vals, obsfile=datamodel.obsfile)
