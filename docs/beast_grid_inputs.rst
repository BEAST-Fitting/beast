#################
BEAST Grid Inputs
#################

Below are details regarding grid choices for stellar evolution models,
spectral libraries, and dust extinction laws.

Stellar Evolution Models
========================

Stellar evolution models provide the BEAST model grid with stellar parameters
(L, T\ :sub:`eff`, log g) as a function of age, metallicity, and stellar mass.  The
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
(surface flux as function of wavelength) as a function of T\ :sub:`eff`, log g, and
metallicity.  The BEAST currently uses theoretical spectral libraries stored
as static libraries.  Model files store spectra in flux format with units of
erg/s/cm2/AA.

Choices:

* `Kurucz`_
    ATLAS9 model atmospheres by `Castelli & Kurucz (2004; CK04) <https://ui.adsabs.harvard.edu/abs/2004A%26A...419..725C/abstract>`_, the industry
    standard for LTE stellar atmospheres that spans a wide range of L and
    T\ :sub:`eff` parameter space and a broad range of wavelengths (0.10-10.0 micron) at low
    spectral resolution (20 Ang/pix; 2x sampling of R~100 at 4000 Ang).
    Z grid spans -2.5 to +0.5 at 0.5 dex sampling. Models are computed using
    solar abundances from `Grevesse & Sauval (1998; Z=0.0169) <https://ui.adsabs.harvard.edu/abs/1998SSRv...85..161G/abstract>`_.

* `Tlusty`_
    Non-LTE hot star (T\ :sub:`eff` > 15,000 K) atmosphere models (OSTAR and BSTAR) by
    `Lanz & Hubeny (2003 <https://ui.adsabs.harvard.edu/abs/2003ApJS..146..417L/abstract>`_, `2007 <https://ui.adsabs.harvard.edu/abs/2007ApJS..169...83L/abstract>`_), using spectra computed for
    Cloudy
    by Peter van Hoof at R~900 (sampled at R~1800). Z grid spans 5 pts from
    1/10-2x Solar, plus additional 4 pts from 1/1000-1/30x Solar for OSTAR
    grid.  Models are computed using solar abundances from `Grevesse & Sauval (1998; Z=0.0169) <https://ui.adsabs.harvard.edu/abs/1998SSRv...85..161G/abstract>`_.

* BTSettl
   * PHOENIX model atmospheres by `Allard et al. 2016 <https://ui.adsabs.harvard.edu/abs/2016sf2a.conf..223A/abstract>`_, specializing in cool star
     atmospheres (T\ :sub:`eff` < 6000 K). Intrinsically high-res, resampled to
     2 Ang/pix grid (medres) and to match CK04 grid (lores). Models are
     computed using solar abundances from `Asplund et al. 2009 <https://ui.adsabs.harvard.edu/abs/2009ARA%26A..47..481A/abstract>`_ (Z=0.0134).
   * Adjustable Parameter = ``medres`` (default=True): 2 Ang/pix resolution,
     or False for CK04 matched wavelength grid.

* `BOSZ`_
    Future Addition -- ATLAS9 model atmospheres computed by `Bohlin et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017AJ....153..234B/abstract>`_
    providing enhancement to CK04 in terms of spectral resolution, wavelength
    coverage, grid density.

* `Munari`_
    ATLAS9 model atmospheres available at higher spectral resolution than
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

* Generalized_RvFaLaw (recommended)
   * allows for any choice of the extinction curve model for the A and B components
   * recommended
      * A: F19 - based on spectroscopy in UV and optical and does a bit better than F99 in the optical
      * B: G03_SMCBar - best average for the SMC "bumpless" extinction curve
      * In beast_settings file:
        extLaw = extinction.Generalized_RvFALaw(ALaw=extinction.Generalized_DustExt(curve='F19'), BLaw=extinction.Generalized_DustExt(curve='G03_SMCBar'))

* `dust_extinction`_ Package Extinction Curves (recommended)
   * ``extinction.Generalized_DustExt()``: Wrapper for any extinction curve
     class available via dust_extinction python package.
   * Select curve class via string parameter: `curve`
   * Example call, for Fitzpatrick+19: ``extinction.Generalized_DustExt('F19')``
   * For the most possible models, see `dust_extinction docs <https://dust-extinction.readthedocs.io/en/stable/>`_

* `beast.physicsmodel.dust.extinction_extension` (recommended of sub 912 A extinction needed)
   * all the observation based models stop at or before 912 A as hydrogen dominates extinction below this wavelength
   * for work that interested in sub 912 A information (e.g., ionizing photon measurements), these models extend the
     dust extinction to shorter than 912 A wavelengths by smoothly merging dust grain models with an observed extinction model
   * `F19_D03_extension`: Extension of Fitzpatrick+19 Milky Way Rv dependent model with Draine03 grain models.
   * `G03_SMCBar_WD01_extension`: Extension of Gordon+03 SMCBar average with the Weingarter & Draine01 SMCBar grain model.

* Gordon+16 (original mixture dust extinction model)
   * ``extinction.Gordon16_RvFALaw()``: Mixture model of MW (Type A;
     Fitzpatrick 99) and SMC (Type B; Gordon+03) extinction curves.
   * Adjustable parameters include: A_V, R_V, and f_A.

* Fitzpatrick 99 (may be deprecated, superceded by dust_extinction package)
   * ``extinction.Fitzpatrick99()``: Default Milky Way extinction curve model.
   * Adjustable parameters: R_V

* Gordon+03 (may be deprecated, superceded by dust_extinction package)
   * ``extinction.Gordon03_SMCBar()``: Empirically-derived SMC Bar dust
     extinction curve.
   * Adjustable A_V, R_V fixed at 2.74

* Cardelli, Clayton, and Mathis 89 (may be deprecated, superceded by dust_extinction package)
   * ``extinction.Cardelli89()``: Well-known Milky Way extinction curve model,
     but advise use of ``dust_extinction F19`` model instead.
   * Adjustable parameters: R_V

 .. _TLusty: http://tlusty.oca.eu/
 .. _Munari: https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/A%2bA/442/1127
 .. _Kurucz: http://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas
 .. _BOSZ: https://archive.stsci.edu/prepds/bosz/
 .. _dust_extinction: https://dust-extinction.readthedocs.io/
