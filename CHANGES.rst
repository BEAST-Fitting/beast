2.2 (unreleased)
================
- no changes yet

2.1 (2025-05-16)
================

- updates to support megabeast (refactoring dust priors)
- added ACS SBC, WFC3/ACS narrow, and JWST bands
- improvements source density maps (diffspike option)
- added log A(V) spacing option
- added 2D dust prior step model that also depends on distance
- added absexponential distance prior
- various updates for upstream package changes
- documentation improvements

2.0 (2022-11-11)
================

- fix bug in filter/vega files for HST bands
- improvements to AST generation
- improvements to simulated catalog generation
- improvements to splitting and merging for subgrids
- merging of stellar and dust prior code (support for megabeast dev)
- updates to batch file generation for large runs
- code restructuring to reduce complexity and duplication
- migrated to the new package template (APE17)
- moved from travis to github actions for CI
- moved from coveralls to codecov for test coverage reporting
- updates to testing structure and code
- various documentation updates and bugfixes

1.4 (2020-5-15)
===============

- adjustable stellar priors
- 2D pPDFs added
- improvements to AST generation
- AGB stellar atmosphere models
- tool scripts expanded and documented
- some scripts system callable (e.g., 'beast get_libfiles')
- examples moved to separate repository
- automated testing expanded
- various documentation updates and bugfixes

1.3 (2019-08-20)
================

- python 2 support dropped
- splitting run_beast.py up to be more pythonic
- AGB model atmospheres added
- generalized dust extinction curves (via dust_extinction package)
- can simulate observations based on physics/observation model
- observationmodel split by background or source density
- very simple splinter observationmodel added
- symlog used to support negative values in log pPDFs
- various documentation updates and bugfixes

1.2 (2018-06-22)
================

- systemic velocity of galaxy can be specified
- symlog used for flux pPDFs of with negative values
- ASTs can be based on source density or surface brightness
- model grid can be split up to allow parallel execution
   - resulting pPDFs can be reassembled
- various bugfixes

1.1 (2018-04-28)
================

- pip installable version
- library files in BEAST_PATH, ~/.beast, or in src directory beast/beast/libs
- automated regression testing
- distance added as 7th variable
- optimized generation of ASTs for toothpick model
- ASTs locations in observed images can be based on residual surface brightness image
- various bugfixes

1.0 (2017-12-29)
================

- initial release
- fully functioning code
