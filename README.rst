= SED FITTER v0.1 =

The SED fitter at this point is basically a *fast likelihood estimator*.


== TODO LIST ==

* find a way to compress the outputs (something like indexing sparse matrix)
* add the observation class (grab from phat cluster analysis)
* include the prior class (from the phat cluster analysis)

== GIT ==

* download the repo
> git clone git@chex.astro.washington.edu:sedfitter.git

* push changes

> git add <name of the file that changed>  # select file(s) that will be associated with the same commit message.
> git commit

== Get Stellar Libraries ==

make libs (downloads CK04 & [TBA] TLUSTY model atmospheres on padova isochrones)y

== Quick start ==

_See testunits.py for a complete example_

> from anased import *
> from extinction import Cardelli

#define some filters
# normalized names correspond to libs/filters.hd5-->root.content.fnames
# basically <FACILITY>_<INSTRUMENT>_<FILTERNAME> (in caps)
#  e.g.: 'HST_WFC3_F336W', 'GROUND_JOHNSON_U'...

> filters = 'hst_wfc3_f225w hst_wfc3_f336w hst_acs_hrc_f475w hst_acs_hrc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

# define the model grid
#  a collection of seds/spectra following the internal grid definition

> g = grid.FileSpectralGrid('libs/SEDs_basel_padovaiso.fits')
> lamb = g.lamb    #wavelenghts associated to the models

# define Av variations
> Av = numpy.arange(0., 3., 0.1)
> Av_law = Cardelli()

# get Data
# soon will be a observation class
#  > obs = observations.myObs(<...>)
#  > f, e , m = obs.getObs(<...>)
#
# Data inputs must be 3 vectors:
#   f   ndarray[float, ndim=1]  a list of fluxes 
#   e   ndarray[float, ndim=1]  a list of associated errors
#   m   ndarray[float, ndim=1]  a mask array for bad values (True = "masked") 
#
# fluxes are assumed to correspond to the grid.seds or a sub sample of it, which
# will be given to the likelihood functions

> f, e, m = getData(<...>)

# Do the computations

> r = numpy.empty( ( g.seds.shape[0], len(Av) ), dtype=float)
> for k in range(len(Av)):
  ...  r[:, k] = job(lamb[:], numpy.copy(f), numpy.copy(e), m, numpy.copy(g.seds), oAv, Av=Av[k], Rv=3.1)
 
# r contains the log-likelihood values for each model in the grid, g, per value in Av



