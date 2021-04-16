"""
Defines a generic interface to observation catalog
"""
import numpy as np
from numpy.random import default_rng

from astropy.table import Table, Column

from beast.observationmodel.vega import Vega
from beast.physicsmodel.priormodel import PriorAgeModel, PriorMassModel
from beast.physicsmodel.grid_weights_stars import compute_bin_boundaries

__all__ = ["Observations", "gen_SimObs_from_sedgrid"]


class Observations(object):
    """
    A generic class that interfaces observation catalog in a standardized way

    Attributes
    ----------
    inputFile : str
        catalog source file
    filters : list
        list of filter names (internal standards)
    filter_aliases : dict
        alias of filter names between internal and external names
    desc : str
        description of the observations
    badvalue : float
        value that tags a bad measurement that should not be used in the
        fitting
    nObs : int
        number of observations in the catalog
    """

    def __init__(
        self, inputFile, filters, obs_colnames=None, vega_fname=None, desc=None
    ):
        """
        Parameters
        ----------
        inputFile : str
            observation file
        filters : list
            interal filter names of the data
        obs_colnames : list, optional
            filter names in the observed catalog
        vega_fname : str, optional
            name of the file with the vega model spectrum
        desc : str, optional
            description of the observations
        """
        if desc is None:
            self.desc = "GENERIC: %s" % inputFile
        else:
            self.desc = desc
        self.inputFile = inputFile
        self.setFilters(filters)
        self.filter_aliases = {}
        for ik, k in enumerate(filters):
            self.filter_aliases[k] = obs_colnames[ik]
        self.readData()
        self.setVegaFluxes(filters, vega_fname=vega_fname)
        # some bad values smaller than expected
        # in physical flux units
        self.setBadValue(6e-40)

    @property
    def nObs(self):
        return len(self.data)

    def __len__(self):
        return self.nObs

    def __call__(self):
        """ Calling the object will show info """
        self.info()

    def info(self):
        """ Prints some information about the catalog """
        txt = "Data read from {s.inputFile:s}\n"
        if self.desc is not None:
            txt += "Description: {s.desc:s}\n"
        txt += "Number of records: {s.nObs:d}\n\n"
        txt += "Dataset contains:"

        print("Data read from %s " % self.inputFile)
        if self.desc is not None:
            print("Description: %s" % self.desc)
            print("Number of records: %d" % self.nObs)
            print("")
            print("Dataset contains:")

        for k in list(self.data.keys()):
            txt += "\t {0:s}\n".format(k)

        if self.filters is None:
            txt += "\n No filters given yet!"
        else:
            txt += "\n Using Filters: {s.filters}\n"

        print(txt.format(s=self))

    def __getitem__(self, *args, **kwargs):
        """ get item will generate a subsample """
        return self.data.__getitem__(*args, **kwargs)

    def keys(self):
        """ Returns dataset content names """
        return self.data.keys()

    def setDescription(self, txt):
        self.desc = txt

    def setBadValue(self, val):
        self.badvalue = val

    def getFilters(self):
        return self.filters

    def setFilters(self, filters):
        self.filters = filters

    def setVegaFluxes(self, filters, vega_fname=None):
        """
        Set vega reference fluxes for conversions

        Parameters
        ----------
        filters : list
            list of filters using the internally normalized namings
        vega_fname : str, optional
            name of the file with the vega model spectrum
        """
        # for optimization purpose: pre-compute
        with Vega(source=vega_fname) as v:
            _, vega_flux, _ = v.getFlux(filters)
        self.vega_flux = vega_flux

    def getFlux(self, num, units=False):
        """
        Flux of an observation computed from normalized vega fluxes

        Parameters
        ----------
        num : int
            index of the star in the catalog to get measurement from
        units : bool
            if set returns the fluxes with units

        Returns
        -------
        flux : ndarray[dtype=float, ndim=1]
            Measured integrated flux values throughout the filters
            in erg/s/cm^2/A
        """
        if self.vega_flux is None:
            raise ValueError("vega_flux not set, can't return fluxes")

        # case for using '_flux' result
        d = self.data[num]

        flux = (
            np.array([d[self.filter_aliases[ok]] for ok in self.filters])
            * self.vega_flux
        )

        if units is True:
            return flux * units.erg / (units.s * units.cm * units.cm * units.angstrom)
        else:
            return flux

    def getFluxerr(self, num):
        """returns the error on the flux of an observation from the number of
        counts (not used in the analysis)"""

        fluxerr = np.zeros(len(self.filters), dtype=float)

        for ek, ok in enumerate(self.filters):
            fluxerr[ek] = self.data[ok + "_err"][num]

        return fluxerr

    def getObs(self, num=0):
        """ returns the flux"""
        if self.filters is None:
            raise AttributeError("No filter set provided.")

        flux = self.getFlux(num)

        return flux

    def readData(self):
        """ read the dataset from the original source file """

        if isinstance(self.inputFile, str):
            self.data = Table.read(self.inputFile)
        else:
            self.data = self.inputFile

    def iterobs(self):
        """ yield getObs """
        for k in range(self.nObs):
            yield self.getObs(k)

    def enumobs(self):
        for k in range(self.nObs):
            yield k, self.getObs(k)


def gen_SimObs_from_sedgrid(
    sedgrid,
    sedgrid_noisemodel,
    nsim=100,
    compl_filter="max",
    complcut=None,
    magcut=None,
    ranseed=None,
    vega_fname=None,
    weight_to_use="weight",
    age_prior_model=None,
    mass_prior_model=None,
):
    """
    Generate simulated observations using the physics and observation grids.
    The priors are sampled as they give the ensemble model for the stellar
    and dust distributions (IMF, Av distribution etc.).
    The physics model gives the SEDs based on the priors.
    The observation model gives the noise, bias, and completeness all of
    which are used in simulating the observations.

    Currently written to only work for the toothpick noisemodel.

    Parameters
    ----------
    sedgrid: grid.SEDgrid instance
        model grid

    sedgrid_noisemodel: beast noisemodel instance
        noise model data

    nsim : int
        number of observations to simulate

    compl_filter : str
        Filter to use for completeness (required for toothpick model).
        Set to max to use the max value in all filters.

    complcut : float (defualt=None)
        completeness cut for only including model seds above the cut
        where the completeness cut ranges between 0 and 1.

    magcut : float (defualt=None)
        faint-end magnitude cut for only including model seds brighter
        than the given magnitude in compl_filter.

    ranseed : int
        used to set the seed to make the results reproducable,
        useful for testing

    vega_fname : string
        filename for the vega info, useful for testing

    weight_to_use : string (default='weight')
        Set to either 'weight' (prior+grid), 'prior_weight', 'grid_weight',
        or 'uniform' (this option is valid only when nsim is supplied) to
        choose the weighting for SED selection.

    age_prior_model : dict
        age prior model in the BEAST dictonary format

    mass_prior_model : dict
        mass prior model in the BEAST dictonary format

    Returns
    -------
    simtable : astropy Table
        table giving the simulated observed fluxes as well as the
        physics model parmaeters
    """
    n_models, n_filters = sedgrid.seds.shape
    flux = sedgrid.seds

    # get the vega fluxes for the filters
    _, vega_flux, _ = Vega(source=vega_fname).getFlux(sedgrid.filters)

    # cache the noisemodel values
    model_bias = sedgrid_noisemodel["bias"]
    model_unc = np.fabs(sedgrid_noisemodel["error"])
    model_compl = sedgrid_noisemodel["completeness"]

    # only use models that have non-zero completeness in all filters
    # zero completeness means the observation model is not defined for that filters/flux
    # if complcut is provided, only use models above that completeness cut
    if complcut is not None:
        finalcomplcut = complcut
    else:
        finalcomplcut = 0.0

    ast_defined = model_compl > finalcomplcut
    sum_ast_defined = np.sum(ast_defined, axis=1)
    goodobsmod = sum_ast_defined >= n_filters

    # completeness from toothpick model so n band completeness values
    # require only 1 completeness value for each model
    # max picked to best "simulate" how the photometry detection is done
    if compl_filter.lower() == "max":
        model_compl = np.max(model_compl, axis=1)
    else:
        short_filters = [
            filter.split(sep="_")[-1].upper() for filter in sedgrid.filters
        ]
        if compl_filter.upper() not in short_filters:
            raise NotImplementedError(
                "Requested completeness filter not present:"
                + compl_filter.upper()
                + "\nPossible filters:"
                + "\n".join(short_filters)
            )

        filter_k = short_filters.index(compl_filter.upper())
        print("Completeness from %s" % sedgrid.filters[filter_k])
        model_compl = model_compl[:, filter_k]

    # if magcut is provided, only use models brighter than the magnitude cut
    # in addition to the non-zero completeness criterion
    if magcut is not None:
        fluxcut_compl_filter = 10 ** (-0.4 * magcut) * vega_flux[filter_k]
        goodobsmod = (goodobsmod) & (flux[:, filter_k] >= fluxcut_compl_filter)

    # initialize the random number generator
    rangen = default_rng(ranseed)

    # if the age and mass prior models are given, use them to determine the
    # total number of stars to simulate
    model_indx = np.arange(n_models)
    if (age_prior_model is not None) and (mass_prior_model is not None):
        nsim = 0
        # logage_range = [min(sedgrid["logA"]), max(sedgrid["logA"])]
        mass_range = [min(sedgrid["M_ini"]), max(sedgrid["M_ini"])]

        # compute the total mass and average mass of a star given the mass_prior_model
        nmass = 100
        masspts = np.logspace(np.log10(mass_range[0]), np.log10(mass_range[1]), nmass)
        mass_prior = PriorMassModel(mass_prior_model)
        massprior = mass_prior(masspts)
        totmass = np.trapz(massprior, masspts)
        avemass = np.trapz(masspts * massprior, masspts) / totmass

        # compute the mass of the remaining stars at each age and
        # simulate the stars assuming everything is complete
        gridweights = sedgrid[weight_to_use]
        gridweights = gridweights / np.sum(gridweights)

        grid_ages = np.unique(sedgrid["logA"])
        age_prior = PriorAgeModel(age_prior_model)
        ageprior = age_prior(grid_ages)
        bin_boundaries = compute_bin_boundaries(grid_ages)
        bin_widths = np.diff(10 ** (bin_boundaries))
        totsim_indx = np.array([], dtype=int)
        for cage, cwidth, cprior in zip(grid_ages, bin_widths, ageprior):
            gmods = sedgrid["logA"] == cage
            cur_mass_range = [
                min(sedgrid["M_ini"][gmods]),
                max(sedgrid["M_ini"][gmods]),
            ]
            gmass = (masspts >= cur_mass_range[0]) & (masspts <= cur_mass_range[1])
            curmasspts = masspts[gmass]
            curmassprior = massprior[gmass]
            totcurmass = np.trapz(curmassprior, curmasspts)

            # compute the mass remaining at each age -> this is the mass to simulate
            simmass = cprior * cwidth * totcurmass / totmass
            nsim_curage = int(round(simmass / avemass))

            # simluate the stars at the current age
            curweights = gridweights[gmods]
            if np.sum(curweights) > 0:
                curweights /= np.sum(curweights)
                cursim_indx = rangen.choice(
                    model_indx[gmods], size=nsim_curage, p=curweights
                )

                totsim_indx = np.concatenate((totsim_indx, cursim_indx))

                nsim += nsim_curage
            # totsimcurmass = np.sum(sedgrid["M_ini"][cursim_indx])
            # print(cage, totcurmass / totmass, simmass, totsimcurmass, nsim_curage)

        totsimmass = np.sum(sedgrid["M_ini"][totsim_indx])
        print(f"number total simulated stars = {nsim}; mass = {totsimmass}")
        compl_choice = rangen.random(nsim)
        compl_indx = model_compl[totsim_indx] >= compl_choice
        sim_indx = totsim_indx[compl_indx]
        totcompsimmass = np.sum(sedgrid["M_ini"][sim_indx])
        print(f"number of simulated stars w/ completeness = {len(sim_indx)}; mass = {totcompsimmass}")

    else:  # total number of stars to simulate set by command line input

        if weight_to_use == "uniform":
            # sample to get the indices of the picked models
            sim_indx = rangen.choice(model_indx[goodobsmod], nsim)

        else:
            gridweights = sedgrid[weight_to_use][goodobsmod] * model_compl[goodobsmod]
            gridweights = gridweights / np.sum(gridweights)

            # sample to get the indexes of the picked models
            sim_indx = rangen.choice(model_indx[goodobsmod], size=nsim, p=gridweights)

        print(f"number of simulated stars = {nsim}")

    # setup the output table
    ot = Table()
    qnames = list(sedgrid.keys())
    # simulated data
    for k, filter in enumerate(sedgrid.filters):
        simflux_wbias = flux[sim_indx, k] + model_bias[sim_indx, k]

        simflux = rangen.normal(loc=simflux_wbias, scale=model_unc[sim_indx, k])

        bname = filter.split(sep="_")[-1].upper()
        fluxname = f"{bname}_FLUX"
        colname = f"{bname}_RATE"
        magname = f"{bname}_VEGA"
        ot[fluxname] = Column(simflux)
        ot[colname] = Column(ot[fluxname] / vega_flux[k])
        pindxs = ot[colname] > 0.0
        nindxs = ot[colname] <= 0.0
        ot[magname] = Column(ot[colname])
        ot[magname][pindxs] = -2.5 * np.log10(ot[colname][pindxs])
        ot[magname][nindxs] = 99.999

        # add in the physical model values in a form similar to
        # the output simulated (physics+obs models) values
        # useful if using the simulated data to interpolate ASTs
        #   (e.g. for MATCH)
        fluxname = f"{bname}_INPUT_FLUX"
        ratename = f"{bname}_INPUT_RATE"
        magname = f"{bname}_INPUT_VEGA"
        ot[fluxname] = Column(flux[sim_indx, k])
        ot[ratename] = Column(ot[fluxname] / vega_flux[k])
        pindxs = ot[ratename] > 0.0
        nindxs = ot[ratename] <= 0.0
        ot[magname] = Column(ot[ratename])
        ot[magname][pindxs] = -2.5 * np.log10(ot[ratename][pindxs])
        ot[magname][nindxs] = 99.999

    # model parmaeters
    for qname in qnames:
        ot[qname] = Column(sedgrid[qname][sim_indx])

    return ot
