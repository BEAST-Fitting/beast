Beast TODO list
===============

External libs
-------------
 * [x] add -- ezpadova
 		added frozen version to beast.core

 * [x] add -- ezpipe
                added submodule to beast.external
 
 * [ ] rm -- mytables need to be removed completely (use eztables instead)

Stellib
-------
 * [x] bug/feature correction -- swap boundary path axis from (g vs T) to (T vs g)
  		added keyword to get_boundaries, plot_boundary

HDFStore
--------
 * [ ] add -- create a HDF Table like class which offers proxy to common eztables.Table functions: 
           e.g., `aliases`, `selectWhere`, `evalexpr`

PosteriorProxy
--------------
 * [ ] add -- Make an easy way to prior definition from the models (set `self.pvalues`)

PosteriorResults
----------------
 * [ ] add -- 1d pdf plots

Documentation
-------------
 * [ ] create -- create an API documentation. Idea: use Sphinx for auto nice website.


Project/morgan
--------------

 models.py:
 * [x] rm -- models.make_spectra: remove extLaw extra arg.
 * [x] update -- update proper documentation (should be ready for sphinx)

 fit.py:
 * [x] update -- update the code to include grid backend feature
                    made cache the default. All functions can use a SpectraGrid object
 * [x] update -- documentation

 fake.py:
 * [ ] update -- FakeData class should return SN errors (currently N only)
 
