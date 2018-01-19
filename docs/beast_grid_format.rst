############
File Formats
############

physicsmodel grid file
======================

Three datasets are present:

   * grid: parameters of the seds (see below)
      - table N parameters x M models
   * lamb: wavelengths of bands
      - vector X bands
   * seds: fluxes in the requested bands [ergs/cm^2/s/A]
      - table X bands x M models

Grid Parameters
---------------

stellar parameters
^^^^^^^^^^^^^^^^^^

Direct Grid Parameters:
   * M_ini: initial mass [M_sun]
   * logA: stellar age in [log10(years)]
   * Z: metallicity [what units/convenction?]

Ancillary Parameters:
   * logL: integrated luminosity of stars [log(???units???)]
   * logT: stellar atmosphere T_eff [log(K)]
   * logg: stellar atmosphere log(g) [log(cm^2/s)???]
   * radius: stellar radius [R_sun????]
   * M_act: actual mass at current age [M_sun]
   * mbolmag: M(bol) (??more info??) [mag]
   * osl: ????

dust extinction parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

Direct Grid Parameters:
   * A(V): extinction in V band [mag]
   * R(V): total-to-selective extinciton = A(V)/E(B-V)
   * f_A: mixture fraction between "MW" and "SMC" extinction curves

Ancillary Parameters:

weights
^^^^^^^

   * weight: combined grid and prior weights used directly in summation
     marginalization
   * grid_weight: weighting in summation marginalization for flat priors
   * prior_weight: weighting in the summation marginalization for input
     priors


model fluxes
^^^^^^^^^^^^

The model fluxes are stored in log10 form with
and without dust extinction

Examples:
   * logHST_ACS_WFC_F475W_wd: flux in ACS/F475W band
     [log(ergs/cm^2/s/A)]
   * logHST_ACS_WFC_F475W_nd: intrinsic flux in ACS/F475W band
     [log(ergs/cm^2/s/A)]


traceback indices
^^^^^^^^^^^^^^^^^

These parameters are useful in mapping the SED model back to the full
grid or spectral grid.  For example, the SED model may be trimmed of
points that will never fit the data due to survey sensitivity limits.

   * fullgrid_idx: index of model in full SED grid
   * specgrid_indx: index of model in the spectral grid

misc
^^^^

   * keep: True if the model is instide the stellar atmosphere grid
     defined in T_eff and log(g) space
   * stage: ???
