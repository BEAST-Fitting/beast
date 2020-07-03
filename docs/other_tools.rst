#########################
BEAST Miscellaneous Tools
#########################

The following miscellaneous tools are useful for handling data formats and
doing assorted data characterization.

.. _other_beast_tools:

Convert hdf5 to FITS
---------------------

Although the BEAST can read data in various formats, it may be of interest to
convert a photometric catalog from hdf5 to FITS format. Below, we show an example for
to use the BEAST tool convert_hdf5_to_fits to convert a catalog from hdf5 to FITS
and save it to disk:

.. code-block:: python

  >>> from beast.tools import convert_hdf5_to_fits #doctest: +SKIP
  >>>
  >>> # Specify the HDF5 file name
  >>> phot_file = 'my_awesome_catalog.hdf5' #doctest: +SKIP
  >>>
  >>> # Convert the HDF5 file to a FITS file (yay!) and save it to disk
  >>> convert_hdf5_to_fits.st_file(file_name = phot_file) #doctest: +SKIP


Observation Depth
-----------------

The noise model contains completeness information for each filter.  The
`calculate_depth` tool uses that to find the Vega magnitude (or flux in
erg/s/cm^2/A) at which a given completeness is reached.  This is useful for
evaluating the depth of your observations.

.. code-block:: python

  >>> from beast.tools import calculate_depth #doctest: +SKIP
  >>>
  >>> # Find the 50% and 75% completeness for phat_small example
  >>> depth = calculate_depth.calculate_depth(  #doctest: +SKIP
          'beast_example_phat_seds.grid.hd5',
          'beast_example_phat_noisemodel.grid.hd5',
          completeness_value=[0.5, 0.75],
          vega_mag=True
      )
  >>> # Depth in F275W (Vega mag)
  >>> depth['HST_WFC3_F275W'] #doctest: +SKIP
  [25.000309202589012, 24.80610510139205]
  >>>
  >>> # Depth in F814W (Vega mag)
  >>> # NaNs show that ASTs don't go deep enough to evaluate 50% completeness
  >>> depth['HST_ACS_WFC_F814W'] #doctest: +SKIP
  [nan, 24.368742437736692]
