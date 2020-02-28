####################
BEAST Miscellaneous Tools
####################

The following miscellaneous tools are useful for handling data formats and

Convert HDF5 to FITS
---------------------

Although the BEAST can read data in various formats, it may be of interest to
convert a photometric catalog from HDF5 to FITS. Below, we show an example for
to use the BEAST tool convert_hdf5_to_fits to do this:

.. code-block:: python

  >>> from beast.tools import convert_hdf5_to_fits #doctest: +SKIP
  >>>
  >>> # Specify the HDF5 file name
  >>> phot_file = 'my_awesome_catalog.hdf5'
  >>>
  >>> # Convert the HDF5 file to a FITS file (yay!) and save it to disk
  >>> convert_hdf5_to_fits.st_file(file_name = phot_file)
