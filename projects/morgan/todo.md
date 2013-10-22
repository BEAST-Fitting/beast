Beast TODO list
===============

External libs
-------------
 * [x] add -- ezpadova
 		added frozen version to beast.core

 * [x] add -- ezpipe
                added submodule to beast.external
 
 * [ ] rm -- mytables need to be removed completely (use eztables instead)

 * [x] add -- matplotlib nxutils removed in version > 1.2
              Path.contains_points do not exists version < 1.2
	      added core.future which uses normal Path if v > 1.2 and a derived
	      version otherwise.

ModelGrid
---------

 * [ ] bug -- copy method works only for MemoryBackend
              move copy to backends: 
	          * memory --> deepcopy
		  * cache --> deepcopy (whatever is loaded only)
		  * hdf --> copy node/file?

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
 
