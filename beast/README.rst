BEAST Code Organization
=======================

The code is organized around the three main components of the BEAST.

- Physics Model
  - stars (stellar atmosphere and evolutionary models)
  - dust (extinction)
  - grid creation

- Observation Model
  - noise model (derived from ASTs, toothpick, trunchen)
  - photometry derivation from spectroscopy

- Fitting
  - fitting of data with the BEAST
  - trimming the grid for speed based on observations to be fitted
  - generating the 1D PDFs of fit parameters fast
  - fit metrics

In addition, there are a number of other directions containing useful
code.

- Examples
  - example use of the BEAST from grid creation to fitting

- Plotting
  - code to make nice figures of results and model characteristics

- Tests
  - test code -> we need lots more of this!!!

- Tools

To get things running:

- Run 'make libs' or 'python getlibs.py' to get the needed FITS files
  containing stellar atmospheres, filter response functions, etc.
  
    
In Development!
---------------

This code is currently in active development.

License
-------

The BEAST is licensed under a 3-clause BSD style license (see the
``licenses/LICENSE.rst`` file).
