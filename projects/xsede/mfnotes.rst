Beast issues
============

* function names in `:class:Observations` are now inconsistent
        * ex: getMags returns rates (inputs)

* what is `class:PhotCharact`? DOC!
        * independent file should be appropriate
        * too many hard coded details
        * why do we have **symmetric** dispersions? The tests we made clearly
          show skewed distributions.

* checking fits:
  a quick way to check the outputs:
  log(g/gsun) = log(M/Msun) + 4 * log(Teff) - log(L/Lsun)
