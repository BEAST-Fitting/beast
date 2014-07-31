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

from pipeline import run_fit, make_models
import datamodel
import noisemodel
import sys

if __name__ == '__main__':
    if ('-?' in sys.argv[1]) or ('--help' in sys.argv[1]):
        print(__doc__)
    elif '-models' in sys.argv[1]:
        make_models()
    else:
        g = '{project:s}/{project:s}_seds.grid.hd5'.format(project=datamodel.project)
        noisefile = '{project:s}/{project:s}_noisemodel.hd5'.format(project=datamodel.project)
        noise = noisemodel.get_noisemodelcat(noisefile)

        run_fit(datamodel.project, g, noise=noise, obsfile=datamodel.obsfile)
