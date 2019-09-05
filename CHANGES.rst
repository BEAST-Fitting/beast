
1.3.2 (2019-09-5)
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
