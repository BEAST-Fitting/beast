####################
BEAST Analysis Tools
####################

The BEAST has several tools for analyzing the results of the BEAST runs.

Star type probability
---------------------

One use case of the BEAST results is to find a set of stars of a particular
type with some probability.  For instance, the user may want stars that are 90%
likely to be high mass stars.

There are two types of stars currently implemented, though more can easily be
added.  When there are cuts along two parameters, the 2D PDFs are utilized.  If
2D PDFs of the required parameters don't exist, the star type will be skipped.

* Extinguished high mass stars (`star_type_probability.ext_O_star`): Stars above
  a certain mass (default `M_ini` > 10 solar masses) with some minimum
  foreground dust column (default `Av` > 0.5 magnitudes).  These are candidates
  for follow-up UV spectroscopy for extinction curves.
* Dusty AGB stars (`star_type_probability.dusty_agb`): As described in
  :doc:`Known Issues <beast_issues>`, the BEAST does not include models for
  dusty AGB stars, and many of them are incorrectly fit as heavily extinguished
  hot stars.  This finds stars using their extinction (default `Av` > 7) and
  effective temperature (default `logT` between 3.7 and 4.2).

Below, we show an example for extinction curve candidates in phat_small.

.. code-block:: python

  >>> from beast.tools import star_type_probability #doctest: +SKIP
  >>> from astropy.io import fits #doctest: +SKIP
  >>> # calculate probabilities
  >>> star_prob = star_type_probability.star_type_probability( #doctest: +SKIP
          'beast_example_phat_pdf1d.fits',
          'beast_example_phat_pdf2d.fits',
          output_filebase=None,
          ext_O_star_params={'min_M_ini':10, 'min_Av':0.5}
      )
  >>> # stars with >80% likelihood of being extinguished massive stars
  >>> np.where(star_prob['ext_O_star'] > 0.8)[0] #doctest: +SKIP
  array([29, 54])
  >>> # confirm their best fit masses and Av
  >>> with fits.open('beast_example_phat_stats.fits') as hdu: #doctest: +SKIP
          print('Masses:', hdu[1].data['M_ini_p50'][[29,54]])
          print('Av:', hdu[1].data['Av_p50'][[29,54]])
  Masses: [16.84200042 15.23882141]
  Av: [3.98209536 3.99047602]


Spectral type comparison
------------------------

When there are spectrally typed stars in our catalogs, we would like to compare
the BEAST parameters to those inferred from the spectral types.  The code in
`beast/tools/compare_spec_type.py` simplifies this comparison.  It ingests the
photometry catalog, stats catalog, spectral types and coordinates of comparison
stars, and several settings.  Here is a summary of its procedure:

1. Convert the spectral types into an effective temperature (Teff) and surface
   gravity (logg).  This is done by interpolating on tables of Teff (from
   Stellar Spectral Classification; R. Gray & C. Corbally) and logg (from
   Allen's Astrophysical Quantities; A. Cox).
2. For each spectrally-typed star, find all matches in the photometry catalog
   within a user-defined radius (default 1").  Since spectrally-typed stars are
   generally bright, choose the brightest source in the user-defined filter
   (default='F475W').  If no match is found, all output values for this star
   (see below) will be set to `None`.
3. Save the p16, p50, and p84 values of logT and logg for the matched BEAST
   star.  Also calculate the number of standard deviations between the BEAST
   fits and star values (assuming no uncertainty on the star values).
4. If the keyword `output_filebase` is set, the results will be saved into a
   file.  Otherwise, the results will be returned in a dictionary.

Below is an example for phat_small.  The spectral types are not from any real
catalog, and are for illustrative purposes only.

.. code-block:: python

  >>> from beast.tools import compare_spec_type #doctest: +SKIP
  >>> spec_match = compare_spec_type.compare_spec_type( #doctest: +SKIP
          'data/b15_4band_det_27_A.fits',
          'beast_example_phat/beast_example_phat_stats.fits',
          []
      )
