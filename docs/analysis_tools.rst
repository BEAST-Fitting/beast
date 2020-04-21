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
2D PDFs of the required parameters don't exist, the probability of that star
type will be `NaN`.

* Extinguished high mass stars (`star_type_probability.ext_O_star`): Stars above
  a certain mass (default `M_ini` > 10 solar masses) with some minimum
  foreground dust column (default `Av` > 0.5 magnitudes).  You may also wish to
  set a maximum `Av` to avoid artifacts (such as the dusty AGB stars).  These
  stars are candidates for follow-up UV spectroscopy for extinction curves.
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
          ext_O_star_params={'min_M_ini':10, 'min_Av':0.5, 'max_Av':5}
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

1. Convert the spectral types into an effective temperature (T_eff) and surface
   gravity (logg).  This is done by interpolating on tables of T_eff (from
   Stellar Spectral Classification; R. Gray & C. Corbally) and logg (from
   Allen's Astrophysical Quantities; A. Cox).
2. For each spectrally-typed star, find all matches in the photometry catalog
   within a user-defined radius (default 1").  Since spectrally-typed stars are
   generally bright, choose the brightest source in the user-defined filter
   (default='F475W').  If no match is found, all output values for this star
   will be set to `None`.
3. Save the p16, p50, and p84 values of logT and logg for the matched BEAST
   star.  Also calculate the number of standard deviations between the BEAST
   fits and star values (assuming no uncertainty on the star values).
4. If the keyword `output_filebase` is set, the results will be saved into a
   file.  Otherwise, the results will be returned in a dictionary.

Below is an example for phat_small.  The spectral types are not from any real
catalog, and are for illustrative purposes only.  The output dictionary shows
the information about the spectrally-typed star, the indices of the matched
BEAST star (in both the photometry and stats catalogs), the BEAST fits, and the
number of standard deviations (sigma) apart they are.  The first star (A2II) is
a good match to BEAST star 27.  The second star (G7II) is a close, but
imperfect, match to BEAST star 8.  The third star (B4V) is outside the catalog
footprint and therefore has no match.

.. code-block:: python

  >>> from beast.tools import compare_spec_type #doctest: +SKIP
  >>> compare_spec_type.compare_spec_type( #doctest: +SKIP
          'data/b15_4band_det_27_A.fits',  # Photometry catalog
          'beast_example_phat/beast_example_phat_stats.fits', # Stats catalog
          [11.2335881, 11.23342557, 1.0],  # RA
          [41.9001895, 41.90006316, 1.0],  # Dec
          ['A', 'G', 'B'],                 # Spectral type
          [2, 7, 4]                        # Subtype
          ['II', 'II', 'V'],               # Luminosity class
          match_radius=0.2                 # Match radius (arcsec)
      )
  {'spec_ra': [11.2335881, 11.23342557, 1.0],
   'spec_dec': [41.9001895, 41.90006316, 1.0],
   'spec_type': ['A 2 II', 'G 7 II', 'B 4 V'],
   'spec_teff': [9000.0, 4916.666666666667, None],
   'spec_logg': [2.7164474106543732, 1.7184474106543735, None],
   'phot_cat_ind': [27, 8, None],
   'stats_cat_ind': [27, 8, None],
   'beast_teff_p50': [9046.250020338754, 4528.230977991138, None],
   'beast_teff_p16': [8643.670633196869, 4335.617282355577, None],
   'beast_teff_p84': [9536.391362054928, 4729.401710221546, None],
   'beast_logg_p50': [2.714286917261312, 1.7684285714285717, None],
   'beast_logg_p16': [2.636272525730954, 1.7014832653061227, None],
   'beast_logg_p84': [2.799534708811963, 1.8353738775510207, None],
   'teff_sigma': [-0.11488422362383206, 1.9308757510045778, None],
   'logg_sigma': [0.025343687546173433, -0.7465969411324851, None]}
