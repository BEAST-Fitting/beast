"""
Everything I need to generate a grid, make fake data and fit

I use the pipeline package I wrote in order to clean the syntax and allow more
flexibilities. In particular it will simplifies the managment of intermediate
results or broken jobs.

make models
-----------
the pipeline is sequence of tasks
tasks = ( t_isochrones(**iso_kwargs),  t_spectra(**spec_kwargs), t_seds(filters, **seds_kwargs) )

models = Pipeline('make_models', tasks)
job, (project, grid) = models(project)

The pipeline is equivalent to using individual tasks as follow:
project, grid = project | t_isochrones(**iso_kwargs) | t_spectra(**spec_kwargs) | t_seds(filters, **seds_kwargs)

However it allows more management options, including logs, memoizations etc.

TODO: make better documentation
TODO: make all static code go into a different module
"""

from pipeline import run_fit
import datamodel
import noisemodel
import sys

if __name__ == '__main__':
    if '-?' in sys.argv[1]:
        print(__doc__)
    else:
        g = '{project:s}/{project:s}_seds.grid.hd5'.format(project=datamodel.project)
        ast = noisemodel.get_noisemodelcat(datamodel.astfile)
        run_fit(datamodel.project, g, ast=ast, obsfile=datamodel.obsfile)
