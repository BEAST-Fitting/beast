#################
BEAST Grid Inputs
#################

Below are details regarding grid choices for stellar evolution models,
spectral libraries, and dust extinction laws.

Stellar Evolution Models
========================

Stellar evolution models provide the BEAST model grid with stellar parameters
(L, Teff, log g) as a function of age, metallicity, and stellar mass.  The
BEAST is currently implemented to obtain this information from stellar
isochrone sets.  These isochrones are obtained at run-time via webforms that
yield interpolated data products and stored as CSV files.

Choices:

* Padova/PARSEC
   * The Padova family of stellar evolution models spans multiple generations
     of  work (Bertelli+94, Girardi+00, Girardi+02; Marigo+08, Girardi+10),
     including the most recent set of PARSEC and COLIBRI models (Bressan+12,
     Marigo+17).
   * PARSEC models accessed via
     `CMD webform <http://stev.oapd.inaf.it/cgi-bin/cmd>`_ by ``ezpadova``
     module.
   * Options = ``modeltype`` (default='parsec12s_r14' for PARSEC+COLIBRI):
     isochrone model set - alternate choices=['parsec12s' for PARSEC 1.2S,
     '2010' for Marigo08+Girardi10]; ``filterPMS`` (default=False): remove
     pre-main sequence stars (M < 12 Msun).
   * Details = age: 0 to 13.5 Gyr; metallicity: [M/H] from -2.2 to +0.5
     (0.0001≤Z≤0.02; 122 pts), also super-solar available (Z = 0.03,0.04,0.06);
     solar abundance: Z_sun=0.0152 based on Caffau+11.

* MIST
   * The MIST models (Choi+16) are computed using MESA (Paxton+11,13,15).
   * MIST Isochrone Interpolator accessed via
     `webform <http://waps.cfa.harvard.edu/MIST/interp_isos.html>`_ by
     ``ezmist`` module.
   * Options = ``rotation`` (default='vvcrit0.0'): select rotating or
     non-rotating models - alternate choice=['vvcrit0.4' for initial rotation
     of 0.4 x critical rotation].
   * Details = age: log(Age) from 5.0 to 10.3 (default sample = 0.05 dex);
     metallicity: [Fe/H] from -4 to +0.5 (default sample = 0.25 dex), all
     solar scaled in Version 1.2; solar abundance: Z_sun=0.0142 based on 
     Asplund+09 (Z_proto=0.0142, Z_photo=0.0134).

Spectral Models
===============

Stellar atmosphere models provide the BEAST model grid with stellar spectra
(surface flux as function of wavelength) as a function of Teff, log g, and
metallicity.  The BEAST currently uses theoretical spectral libraries stored
as static libraries.  Model files store spectra in flux format with units of
erg/s/cm2/AA.

Choices:

* `Kurucz`_
   * ATLAS9 model atmospheres by Castelli & Kurucz (2004; CK04), the industry
     standard for LTE stellar atmospheres that spans a wide range of L and
     Teff parameter space and a broad range of wavelenghts (XXX) at low
     spectral resolution (20 Ang/pix; 2x sampling of R~100 at 4000 Ang).
     Z grid spans -2.5 to +0.5 at 0.5 dex sampling. Models are computed using
     solar abundances from Grevesse & Sauval (1998; Z=0.0169).

* `Tlusty`_
   * Non-LTE hot star (Teff > 15,000 K) atmosphere models (OSTAR and BSTAR) by
     Lanz & Hubeny (2003, 2007), using spectra computed for
     `Cloudy <http://nova.astro.umd.edu/Tlusty2002/tlusty-frames-cloudy.html>`_
     by Peter van Hoof at R~900 (sampled at R~1800). Z grid spans 5 pts from
     1/10-2x Solar, plus additional 4 pts from 1/1000-1/30x Solar for OSTAR
     grid.  Models are computed using solar abundances from Grevesse & Sauval
     (1998; Z=0.0169).

* `BTSettl`_
   * PHOENIX model atmospheres by Allard et al., specializing in cool star
     atmospheres (Teff < 6000 K). Intrinsically high-res, resampled to
     2 Ang/pix grid (medres) and to match CK04 grid (lores). Models are
     computed using solar abundances from Asplund+09 (Z=0.0134).
   * Adjustable Parameter = ``medres`` (default=True): 2 Ang/pix resolution,
     or False for CK04 matched wavelength grid.

* `BOSZ`_
   * Future Addition -- ATLAS9 model atmospheres computed by Bohlin+17
     providing enhancement to CK04 in terms of spectral resolution, wavelength
     coverage, grid density.

* `Munari`_
   * ATLAS9 model atmospheres available at higher spectral resolution than
     CK04, but over limited wavelength range.

Recommendations: Tlusty+Kurucz as default, providing non-LTE models at high
temperatures and standard ATLAS9 models elsewhere.  BTSettl library provides
improvement at low temperatures, while BOSZ (coming soon) or Munari provide
higher spectral resolution if required.

Dust Extinction
=================

Dust extinction is applied to model spectra before bandpass convolution,
providing fully self-consistent treatment of reddening.

Choices:

* Gordon+16
   * ``extinction.Gordon16_RvFALaw()``: Mixture model of MW (Type A;
     Fitzpatrick 99) and SMC (Type B; Gordon+03) extinction curves.
   * Adjustable parameters include: A_V, R_V, and f_A.

* Fitzpatrick 99
   * ``extinction.Fitzpatrick99()``: Default Milky Way extinction curve model.
   * Adjustable parameters: R_V

* Gordon+03
   * ``extinction.Gordon03_SMCBar()``: Empirically-derived SMC Bar dust
     extinction curve.
   * Adjustable A_V, R_V fixed at 2.74

* Cardelli, Clayton, and Mathis 89
   * ``extinction.Cardelli89()``: Well-known Milky Way extinction curve model,
     but advise use of ``Fitzpatrick99()`` model instead.
   * Adjustable parameters: R_V

* `dust_extinction`_ Package Extinction Curves
   * ``extinction.Generalized_DustExt()``: Wrapper for any extinction curve
     class available via dust_extinction python package.
   * Select curve class via string parameter: `curve`
   * Example call, for Fitzpatrick 04: ``extinction.Generalized_DustExt('F04')``
   * R_V-dependent models available: Cardelli+89=CCM89, O'Donnell94=O94,
     Fitzpatrick99=F99, Fitzpatrick04=F04, MaizApellaniz+14=M14
   * Average model available: Gordon+03's SMC Bar Avg = G03_SMCBar; Gordon+03's
     LMC Avg = G03_LMCAvg; Gordon+03's LMC2 Supershell Avg = G03_LMC2;
     Gordon, Cartledge, & Clayton 09's MW Avg = GCC09_MWAvg

 .. _BTSettl: https://phoenix.ens-lyon.fr/Grids/BT-Settl/
 .. _TLusty: http://nova.astro.umd.edu/Tlusty2002/database/
 .. _Munari: http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=J%2FA%2BA%2F442%2F1127
 .. _Kurucz: http://www.stsci.edu/hst/observatory/crds/castelli_kurucz_atlas.html
 .. _BOSZ: https://archive.stsci.edu/prepds/bosz/
 .. _dust_extinction: https://dust-extinction.readthedocs.io/
