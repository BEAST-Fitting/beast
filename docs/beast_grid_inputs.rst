.. _beast_grid_inputs:

#################
BEAST Grid Inputs
#################

Below are details regarding grid choices for stellar evolution models,
spectral libraries, and dust extinction laws.

Stellar Evolution Models
========================

Stellar evolution models provide the BEAST model grid with stellar parameters
(L, T\ :sub:`eff`, log g) as a function of age, metallicity, and stellar mass.
The BEAST supports using evolutionary tracks directly or via isochrones.

Evolutionary Tracks
-------------------

* MIST
   * The MIST models (Choi+16) are computed using MESA (Paxton+11,13,15). [TBR]
   * Stored as static FITS files obtained via :ref:`BEAST library files <library-files>`.

Isochrones
----------

These isochrones are obtained at run-time via web queries to the CMD site that
yield interpolated data products and stored as CSV files.

Choices:

* Padova/PARSEC
   * The Padova family of stellar evolution models spans multiple generations
     of  work (Bertelli+94, Girardi+00, Girardi+02; Marigo+08, Girardi+10),
     including the most recent set of PARSEC and COLIBRI models (Bressan+12,
     Marigo+17).
   * PARSEC models accessed via
     `CMD webform <https://stev.oapd.inaf.it/cgi-bin/cmd>`_ by ``ezpadova``
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
     `webform <https://mist.science/interp_isos.html>`_ by
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
as static libraries obtained via :ref:`BEAST library files <library-files>`.
Model files store spectra in flux format with units of erg/s/cm2/AA.

Multiple libraries can be combined and the order defines the preference at the 
same T\ :sub:`eff`, log g, and metallicity.  For example: 
``osl = stellib.Tlusty2025() + stellib.BOSZ2024()``
means that Tlusty2025 models are preferred over BOSZ2024 models for the same stellar parameters.

For choices see :doc:`Stellar Atmosphere Models <atmosphere_models>`.

Recommendations: Tlusty2025+BOSZ2024+Aringer2016 as default.  This provides the widest
coverage with the best models used in different temperature ranges.

Dust Extinction
=================

Dust extinction is applied to model spectra before bandpass convolution,
providing fully self-consistent treatment of reddening.

Choices:

* Generalized_RvFaLaw (recommended)
   * allows for any choice of the extinction curve model for the A and B components
   * recommended
      * A: G23 - based on spectroscopy from the far-UV to the mid-IR
      * B: G24_SMCAvg - best average for the SMC "bumpless" extinction curve
      * In beast_settings file:
    
        * ALaw = extinction.Generalized_DustExt(curve='G23')
        * BLaw = extinction.Generalized_DustExt(curve='G24_SMCAvg')
        * extLaw = extinction.Generalized_RvFALaw(ALaw=ALaw, BLaw=Blaw)

* `dust_extinction`_ Package Extinction Curves (recommended)
   * ``extinction.Generalized_DustExt()``: Wrapper for any extinction curve
     class available via dust_extinction python package.
   * Select curve class via string parameter: `curve`
   * Example call, for Fitzpatrick+19: ``extinction.Generalized_DustExt('G23')``
   * For the most possible models, see `dust_extinction docs <https://dust-extinction.readthedocs.io/en/stable/>`_

* `beast.physicsmodel.dust.extinction_extension` (recommended of sub 912 A extinction needed)
   * all the observation based models stop at or before 912 A as hydrogen dominates extinction below this wavelength
   * For work that interested in sub 912 A information (e.g., ionizing photon measurements), these models extend the
     dust extinction to 100 A by smoothly merging dust grain models with an observed extinction model
     starting at 1675 A.  Thus, the extinction between around 1150 A and 1675 A is a combination of the observed extinction curve and
     the dust grain model extinction curve.  Most observational based extinction curves do not go all the way to 912 A due to lack
     of measurements.
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

 .. _dust_extinction: https://dust-extinction.readthedocs.io/
