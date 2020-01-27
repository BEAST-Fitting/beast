###########
Output file
###########

Below are details regarding the output statistics files produced by the BEAST


Columns in BEAST stats files
============================

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

* M_ini: initial stellar mass (in solar masses)
* logA: log10 of the stellar age (in years)
* Z: stellar metallicity
* M_act: current stellar mass (in solar masses)
* logL: log10 of the stellar luminosity
* logT: log10 of the stellar effective temperature
* logg: log10 of the stellar surface gravity
* mbol: bolometric magnitude
* radius: stellar radius

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
