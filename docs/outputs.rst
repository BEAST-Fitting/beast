###########
Output file
###########

Below are details regarding the output statistics files produced by the BEAST


Columns in BEAST stats files
============================

ID Name
  * IAU suggested naming scheme used (example: PHAT J113759.63+421022.03)

Data Parameters
---------------

RA
  * right ascention from photometry catalog

DEC
  * declination from photometry catalog

field 
  * field in the brick

inside_brick 
  * inside the brick boundaries
inside_chipgap
  * in ACS chip gap

Photometry (units are *flux* [not mag])
  (units are normalized Vega fluxes, e.g., flux/flux_vega)
  
  HST_WFC3_F275W
  HST_WFC3_F336W
  HST_ACS_WFC_F475W
  HST_ACS_WFC_F814W
  HST_WFC3_F110W
  HST_WFC3_F160W

BEAST goodness-of-fit metrics
-----------------------------
Pmax 
  * maxium probability of nD PDF

chi2min 
  * minimum value of chisqr


BEAST Fitting parameters
------------------------

The results come in three flavors

X_Best 
  * best fit value ["traditional" values]

X_Exp 
  * expectation value (average weighted by 1D PDF) [best when not using uncertainties]

X_p50 
  * 50% value from 1D PDF
X_p16 
  * 16% value from 1D PDF (minus 1 sigma)
X_p84 
  * 84% value from 1D PDF (plus 1 sigma)

[best when using uncertainites]
[use p50 - (p50-p16) + (p84 - p50) when quoting results with uncertainties]

Dust Parameters 
---------------

(first 3 primary, others derived)

Av = A(V) 
  * visual extinction in magnitudes
Rv  
  *  R(V) = A(V)/E(B-V) = ratio of total to selective extinction
f_A  
  *  fraction in extinction curve from A component (MW)

Rv_A  
  *  R(V)_A = R(V) of A component of BEAST R(V)-f_A model of extinction curves

Stellar Parameters
------------------

(first 3 primary, others derived)

M_ini  
  * initial stellar mass in solar masses
logA  
  * log10 of the stellar age (in years)
Z  
  * stellar metallicity 

M_act  
  * current stellar mass in solar masses
logL  
  * log10 of the stellar luminosity
logT  
  * log10 of the stellar effective temperature
logg  
  * log10 of the stellar surface gravity
mbol  
  * bolometric magnitude(????)
radius  
  * stellar radius

BEAST Model predicted values
----------------------------

[same flavors as the BEAST fit parameters]

logHST_WFC3_F275W_nd  
  * log10 of the unextinguished WFC3 F275W flux
logHST_WFC3_F275W_wd  
  * log10 of the extinguished WFC3 F275W flux
logHST_WFC3_F336W_nd  
  * log10 of the unextinguished WFC3 F336W flux
logHST_WFC3_F336W_wd  
  * log10 of the extinguished WFC3 F336W flux
logHST_ACS_WFC_F475W_nd  
  * log10 of the unextinguished ACS F475W flux
logHST_ACS_WFC_F475W_wd  
  * log10 of the extinguished ACS F475W flux
logHST_ACS_WFC_F814W_nd  
  * log10 of the unextinguished ACS F814W flux
logHST_ACS_WFC_F814W_wd  
  * log10 of the extinguished ACS F814W flux
logHST_WFC3_F110W_nd  
  * log10 of the unextinguished WFC3 F110W flux
logHST_WFC3_F110W_wd  
  * log10 of the extinguished WFC3 F110W flux
logHST_WFC3_F160W_nd  
  * log10 of the unextinguished WFC3 F160W flux
logHST_WFC3_F160W_wd  
  * log10 of the extinguished WFC3 F160W flux

logGALEX_FUV_nd  
  * log10 of the unextinguished GALEX FUV flux
logGALEX_FUV_wd  
  * log10 of the extinguished GALEX FUV flux
logGALEX_NUV_nd  
  * log10 of the unextinguished GALEX FUV flux
logGALEX_NUV_wd  
  * log10 of the extinguished GALEX FUV flux

logF_UV_6_13e_nd  
  * log10 of the unextinguished flux between 6 and 13 eV
logF_UV_6_13e_wd  
  * log10 of the extinguished flux between 6 and 13 eV

logF_QION_nd  
  * log10 of the unextinguished ionizing flux (***do not use for PHAT results - incorrect***)
logF_QION_wd  
  * log10 of the extinguished ionizing flux (***do not use for PHAT results - incorrect***)

Extras
------

specgrid_indx  
  * index of model in the spectral grid
Pmax_indx  
  * index in BEAST grid of Pmax
chi2min_indx  
  * index in BEAST grid of chi2min
