##################
BEAST Output Files
##################

Below are details regarding the output files produced by the BEAST:

* `*_stats.fits`: Statistics for each of the fitted and derived parameters,
  including the 16th/50th/84th percentiles, mean, and expectation value
* `*_pdf1d.fits`: Marginalized 1D PDFs for each of the fitted and derived
  parameters
* `*_pdf2d.fits`: Marginalized 2D PDFs for pairs of parameters
* `*_lnp.hd5`: Sparsely sampled log likelihoods

Several of the BEAST output files are saved in the hdf5 format, which can be
more challenging to access than fits files.  There are functions in
`tools/read_beast_data.py` to facilitate reading those files.


Statistics file
===============

Data Parameters
---------------

* `ID Name`: IAU suggested naming scheme used (example: PHAT J113759.63+421022.03)
* `RA`: right ascension from photometry catalog
* `DEC`: declination from photometry catalog
* `field`: field in the brick
* `inside_brick`: inside the brick boundaries
* `inside_chipgap`: in ACS chip gap
* Photometry: listed as *flux* (not mag), units are normalized Vega fluxes
  (e.g., flux/flux_vega)

  * `HST_WFC3_F275W`
  * `HST_WFC3_F336W`
  * `HST_ACS_WFC_F475W`
  * `HST_ACS_WFC_F814W`
  * `HST_WFC3_F110W`
  * `HST_WFC3_F160W`

Goodness-of-fit metrics
-----------------------

* `Pmax`: maximum probability of nD PDF
* `Pmax_indx`: index in BEAST model grid corresponding to `Pmax`
* `specgrid_indx`: index in spectroscopic grid corresponding to `Pmax`
* `chi2min`: minimum value of chisqr
* `chi2min_indx`: index in BEAST model grid corresponding to `chi2min`

Fitted and derived parameters
-----------------------------

Each parameter (listed below) has five values associated with it:

* `X_Best`: best fit value ["traditional" values]
* `X_Exp`: expectation value (average weighted by 1D PDF) [best when not using
  uncertainties]
* `X_p50`: 50th percentile from 1D PDF
* `X_p16`: 16th percentile from 1D PDF (p50-p16 is proxy for -1 sigma)
* `X_p84`: 84th percentile from 1D PDF (p84-p50 is proxy for +1 sigma)

Dust Parameters
"""""""""""""""

First 3 primary, others derived

* `Av = A(V)`: visual extinction in magnitudes
* `Rv`: R(V) = A(V)/E(B-V) = ratio of total to selective extinction
* `f_A`: fraction in extinction curve from A component (MW)
* `Rv_A`: R(V)_A = R(V) of A component of BEAST R(V)-f_A model of extinction curves

Stellar Parameters
""""""""""""""""""

First 3 primary, others derived

* `M_ini`: initial stellar mass (in solar masses)
* `logA`: log10 of the stellar age (in years)
* `Z`: stellar metallicity
* `M_act`: current stellar mass (in solar masses)
* `logL`: log10 of the stellar luminosity
* `logT`: log10 of the stellar effective temperature
* `logg`: log10 of the stellar surface gravity
* `mbol`: bolometric magnitude
* `radius`: stellar radius

Predicted Fluxes
""""""""""""""""

The fitting process also predicts fluxes, both in the observed bands and in
other bands of interest.

* `logHST_WFC3_F275W_nd`: log10 of the unextinguished WFC3 F275W flux
* `logHST_WFC3_F275W_wd`: log10 of the extinguished WFC3 F275W flux
* `logHST_WFC3_F336W_nd`: log10 of the unextinguished WFC3 F336W flux
* `logHST_WFC3_F336W_wd`: log10 of the extinguished WFC3 F336W flux
* `logHST_ACS_WFC_F475W_nd`: log10 of the unextinguished ACS F475W flux
* `logHST_ACS_WFC_F475W_wd`: log10 of the extinguished ACS F475W flux
* `logHST_ACS_WFC_F814W_nd`: log10 of the unextinguished ACS F814W flux
* `logHST_ACS_WFC_F814W_wd`: log10 of the extinguished ACS F814W flux
* `logHST_WFC3_F110W_nd`: log10 of the unextinguished WFC3 F110W flux
* `logHST_WFC3_F110W_wd`: log10 of the extinguished WFC3 F110W flux
* `logHST_WFC3_F160W_nd`: log10 of the unextinguished WFC3 F160W flux
* `logHST_WFC3_F160W_wd`: log10 of the extinguished WFC3 F160W flux
* `logGALEX_FUV_nd`: log10 of the unextinguished GALEX FUV flux
* `logGALEX_FUV_wd`: log10 of the extinguished GALEX FUV flux
* `logGALEX_NUV_nd`: log10 of the unextinguished GALEX FUV flux
* `logGALEX_NUV_wd`: log10 of the extinguished GALEX FUV flux
* `logF_UV_6_13e_nd`: log10 of the unextinguished flux between 6 and 13 eV
* `logF_UV_6_13e_wd`: log10 of the extinguished flux between 6 and 13 eV
* `logF_QION_nd`: log10 of the unextinguished ionizing flux (***do not use for
  PHAT results - incorrect***)
* `logF_QION_wd`: log10 of the extinguished ionizing flux (***do not use for
  PHAT results - incorrect***)


1D PDF file
===========

Each extension in the fits file is for one of the parameters listed above.  It
contains an array with dimensions `(N_obs+1, N_bin)`, where `N_obs` is the
number of stars and `N_bin` is the number of bins for that parameter.  Each
entry in the array is the probability (NOT logarithmic) in each bin.  The bin
values are listed in the last line of the array.

Below is an example for `Rv` in the `phat_small` example.

.. code-block:: python

  >>> from astropy.io import fits #doctest: +SKIP
  >>> hdu = fits.open('beast_example_phat_pdf1d.fits') #doctest: +SKIP
  >>> hdu.info() #doctest: +SKIP
  Filename: beast_example_phat_pdf1d.fits
  No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU       6   (2, 2)   float64
  1  Av            1 ImageHDU         8   (11, 270)   float64
  2  M_act         1 ImageHDU         8   (50, 270)   float64
  3  M_ini         1 ImageHDU         8   (50, 270)   float64
  4  Rv            1 ImageHDU         8   (5, 270)   float64
  5  Rv_A          1 ImageHDU         8   (9, 270)   float64
  6  Z             1 ImageHDU         8   (5, 270)   float64
  ...
  >>> hdu['Rv'].data[0,:]  # 1D PDF for star 0 #doctest: +SKIP
  array([0.00000000e+00, 9.99753477e-01, 2.46523236e-04, 0.00000000e+00,
       0.00000000e+00])
  >>> hdu['Rv'].data[-1,:]  # corresponding bin values #doctest: +SKIP
  array([2., 3., 4., 5., 6.])


2D PDF file
===========

Each extension in the fits file is for one of the pairs of fitting parameters
(the default is the 7 main parameters, but the user may have selected a
different set).  The saved arrays have dimensions `(N_obs+2, N_bin_1, N_bin_2)`,
where `N_obs` is the number of stars, `N_bin_1` is the number of bins for the
first parameter, and `N_bin_2` is the number of bins for the second parameter.
The last two slices contain the bin values.

Below is an example of the `Rv` and `f_A` 2D PDF in the `phat_small` example.

.. code-block:: python

  >>> from astropy.io import fits #doctest: +SKIP
  >>> hdu = fits.open('beast_example_phat_pdf2d.fits') #doctest: +SKIP
  >>> hdu.info() #doctest: +SKIP
  Filename: beast_example_phat_pdf2d.fits
  No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU       6   (2, 2)   float64
  1  Av+M_ini      1 ImageHDU         9   (50, 11, 271)   float64
  2  Av+Rv         1 ImageHDU         9   (5, 11, 271)   float64
  3  Av+Z          1 ImageHDU         9   (5, 11, 271)   float64
  4  Av+f_A        1 ImageHDU         9   (4, 11, 271)   float64
  5  Av+logA       1 ImageHDU         9   (5, 11, 271)   float64
  6  M_ini+Rv      1 ImageHDU         9   (5, 50, 271)   float64
  7  M_ini+Z       1 ImageHDU         9   (5, 50, 271)   float64
  ...
  >>> hdu['Rv+f_A'].data[0,:,:]  # 2D PDF for star 0 #doctest: +SKIP
  array([[0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
         [6.86784697e-01, 2.94159452e-01, 1.88093274e-02, 0.00000000e+00],
         [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.46523236e-04],
         [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
         [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00]])
  >>> hdu['Rv+f_A'].data[-2,:,:]  # corresponding Rv bin values #doctest: +SKIP
  array([[2., 2., 2., 2.],
         [3., 3., 3., 3.],
         [4., 4., 4., 4.],
         [5., 5., 5., 5.],
         [6., 6., 6., 6.]])
  >>> hdu['Rv+f_A'].data[-1,:,:]  # corresponding f_A bin values #doctest: +SKIP
  array([[0.25, 0.5 , 0.75, 1.  ],
         [0.25, 0.5 , 0.75, 1.  ],
         [0.25, 0.5 , 0.75, 1.  ],
         [0.25, 0.5 , 0.75, 1.  ],
         [0.25, 0.5 , 0.75, 1.  ]])


Log Likelihood file
===================

(to be added)
