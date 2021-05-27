################
Photometry Files
################

Observed Photometry file
========================

The photometry file used by the BEAST are composed of tables with one row
per source.  The table should have columns
giving the measured fluxes in all photemetric bands to be fit and the ra & dec
coordinates of each source.
The measured fluxes must be supplied in linear physical units, not magnitudes.
All bands to be fit are assumed to have a measured flux and this flux can be
negative.
The BEAST was developed for photometry where the source detection was done
simultaneously in a composite, multi-band image followed by simultaneous
PSF-fitting measurements in all bands.
As such, measured fluxes are allowed to be negative as well as positive and with
any nominal signal-to-noise.

Required columns:

- ra & dec (in deg).
  Column names: `RA` & `DEC`
- Measured fluxes for each band (usually given as rate = flux/flux_vega).
  Example column name `F475W_RATE`

.. note::
   Upper limits in one or more bands are not supported.
   Including such upper limits requires significantly more
   computationally intensive likelihood calculation.

Artificial Star Test (AST) results
==================================

The AST results file used by the BEAST is composed of tables with one row
per artificial star.
The tables should have the same columns as the photometry file *plus* columns
giving the input magnitudes in each band.  Magnitudes are used as the inputs
fluxes are all positive and this is what the program that measures the
photometry expects for AST inputs.
In addition, a `CUT_FLAG` is required to allow for the completeness to be
computed as part of the observation model.  While the completeness is not used
by the BEAST, it is critical for the MegaBEAST calculation.

Required columns:

- ra & dec (in deg).
  Column names: `RA` & `DEC`
- Measured fluxes for each band (usually given as rate = flux/flux_vega).
  Example column name `F475W_RATE`
- Input fluxes for each band (usually given as vega magnitudes).
  Example column name `F475W_VEGA`
- flag if source is cut (1 if the source is cut from the catalog).
  Column name `CUT_FLAG`

Pre-processing of files
=======================

As part of the usual workflow for BEAST fitting, pre-processing of the observed
and AST photometry files is done.
For the observed photometry file, all sources that will be fit must be removed.
Such sources may have
characteristics of known bad sources (e.g., crowding, sharpness, S/N in a band, etc.)
or not have measurements in all bands.
Not having measurements in all bands can be due to not observing that point
on the sky in all bands, saturation in some bands due to a bright source, or
a cosmic ray corrupting the measurement in one or more bands for that source.
Fitting of sources that do not have measurements in all bands is possible by
creating separate observed catalog for each set of sources with the same measurements
in a set of bands.  In other words, if there are sources with measurements in
all but one specific band (e.g., F814W), then the observed catalog with those
sources can be created and then fit with a BEAST model with all bands expect that
band (e.g., F814W).
The results of these different BEAST runs can be combined as the
BEAST fit parameters, uncertainties, and (1d and nD) likelihoods fully reflect
the different number of bands with measurements.

Any processing of the photometry file should also be done for the AST output file,
expect that the sources that do not pass should not be removed, but instead the
`CUT_FLAG` column should be set to `1` instead of `0`.
Doing the same processing on the AST file as was done to the photometry file
ensures that the selection function is the same between the observed and
AST files.  This is needed to allow the BEAST to create an observation model
that correctly includes the full photometry selection function.
