"""
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

from pipeline import make_models
import datamodel
import noisemodel
import sys, argparse
from beast.core import prior_weights, trim_grid
from beast.core.grid import FileSEDGrid
from beast.fitting import fit


if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--models", help="Generate the model grid",
                        action="store_true")
    parser.add_argument("-n", "--noise", help="Calculate the noise model",
                        action="store_true")
    parser.add_argument("-t", "--trim", help="Trim the model and noise grids",
                        action="store_true")
    parser.add_argument("-f", "--fit", help="Fit the observed data",
                        action="store_true")
    parser.add_argument("-r", "--resume", help="Resume a run",
                        action="store_true")
    args = parser.parse_args()

    
    if args.models:
        make_models()

        
    elif args.noise:
        print('Generating noise model from ASTs and absflux A matrix')
        # get the modelsedgrid on which to generate the noisemodel
        modelsedgrid = FileSEDGrid('{project:s}/{project:s}_seds.grid.hd5'.format(project=datamodel.project))
        # generate the AST noise model  
        noisemodel.make_trunchen_noise_model(datamodel.noisefile, datamodel.astfile, 
                                             datamodel.ast_colnames,modelsedgrid, 
                                             generic_absflux_a_matrix= datamodel.generic_absflux_a_matrix)
    
    elif args.trim:
        print('Trimming the model and noise grids')

        # get the modesedgrid on which to generate the noisemodel  
        modelsedgrid = FileSEDGrid('{project:s}/{project:s}_seds.grid.hd5'.format(project=datamodel.project))
        # read in the noise model just created
        noisemodel_vals = noisemodel.get_noisemodelcat(datamodel.noisefile)
        # read in the observed data
        obsdata = datamodel.get_obscat(datamodel.obsfile, datamodel.distanceModulus, datamodel.filters)
        # trim the model sedgrid
        sed_trimname = '{project:s}/{project:s}_seds.grid_trim.hd5'.format(project=datamodel.project)
        noisemodel_trimname = '{project:s}/{project:s}_noisemodel_trim.hd5'.format(project=datamodel.project)

        trim_grid.trim_models(modelsedgrid, noisemodel_vals, obsdata, sed_trimname, noisemodel_trimname, 
                              sigma_fac=3., trunchen=True)

    if args.fit:
        print('Fitting the SEDs')
        
        # the files for the trimmed model grid and noisemodel grid
        g = '{project:s}/{project:s}_seds.grid_trim.hd5'.format(project=datamodel.project)
        noisefile = '{project:s}/{project:s}_noisemodel_trim.hd5'.format(project=datamodel.project)
        statsfile = '{project:s}/{project:s}_memory_stats.fits'.format(project=datamodel.project) 
       
        noise = noisemodel.get_noisemodelcat(noisefile)
        obsdata = datamodel.get_obscat(datamodel.obsfile,datamodel.distanceModulus,datamodel.filters)

        fit.summary_table_memory(obsdata, noise, g, resume=args.resume,
                                 threshold=-10, save_every_npts=100, lnp_npts=60,
                                 stats_outname=statsfile,
                                 pdf1d_outname=statsfile.replace('stats.fits','pdf1d.fits'),
                                 lnp_outname=statsfile.replace('stats.fits','lnp.hd5'),
                                 use_full_cov_matrix=datamodel.use_full_cov_matrix)
