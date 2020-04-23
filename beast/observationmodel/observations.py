"""
Defines a generic interface to observation catalog
"""
import numpy as np

from astropy.table import Table, Column

from beast.observationmodel.vega import Vega

__all__ = ["Observations", "gen_SimObs_from_sedgrid"]


class Observations(object):
    """
    A generic class that interfaces observation catalog in a standardized way

    Attributes
    ----------
    inputFile: str
        catalog source file

    filters: sequence(str)
        list of filter names (internal standards)

    desc: str, optional
        description of the observations

    badvalue: float, optional
        value that tags a bad measurement that should not be used in the
        fitting.

    nObs: int
        number of observations in the catalog
    """

    def __init__(self, inputFile, desc=None):
        """ Generate a data interface object """
        self.inputFile = inputFile
        self.filters = None
        self.desc = desc
        self.readData()
        self.badvalue = None

    @property
    def nObs(self):
        return self.data.nrows

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

    def getMags(self, num, filters):
        raise Exception("Do not use as magnitudes")
        # return np.array([self.data[tt][num] for tt in filters])

    def getErrors(self, num, filters):
        raise Exception("Do not use as magnitudes")
        # return np.array([self.data[tt + "err"][num] for tt in filters])

    def getFlux(self, num):
        """returns the flux of an observation from the number of counts"""

        flux = np.empty(len(self.filters), dtype=float)
        for ek, ok in enumerate(self.filters):
            flux[ek] = self.data[ok][num]

        return flux

    def getFluxerr(self, num):
        """returns the error on the flux of an observation from the number of
        counts (not used in the analysis)"""

        fluxerr = np.empty(len(self.filters), dtype=float)

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
        from ..external.eztables import AstroTable

        if isinstance(self.inputFile, str):
            self.data = AstroTable(self.inputFile)
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
    compl_filter="F475W",
    ranseed=None,
    vega_fname=None,
    weight_to_use='weight',
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
        filter to use for completeness (required for toothpick model)

    ranseed : int
        used to set the seed to make the results reproducable
        useful for testing

    vega_fname : string
        filename for the vega info
        usefule for testing

    weight_to_use : string (default='weight')
        Set to either 'weight' (prior+grid), 'prior_weight', or 'grid_weight' to
        choose the weighting for SED selection.

    Returns
    -------
    simtable : astropy Table
        table giving the simulated observed fluxes as well as the
        physics model parmaeters
    """
    flux = sedgrid.seds
    n_models, n_filters = flux.shape

    # hack to get things to run for now
    short_filters = [filter.split(sep="_")[-1].upper() for filter in sedgrid.filters]
    if compl_filter.upper() not in short_filters:
        raise NotImplementedError(
            "Requested completeness filter not present:"
            + compl_filter.upper()
            + "\nPossible filters:"
            + "\n".join(short_filters)
        )

    filter_k = short_filters.index(compl_filter.upper())
    print("Completeness from %s" % sedgrid.filters[filter_k])

    # cache the noisemodel values
    model_bias = sedgrid_noisemodel["bias"]
    model_unc = np.fabs(sedgrid_noisemodel["error"])
    model_compl = sedgrid_noisemodel["completeness"]

    # the combined prior and grid weights
    # using both as the grid weight needed to account for the finite size
    #   of each grid bin
    # if we change to interpolating between grid points, need to rethink this
    gridweights = sedgrid[weight_to_use] * model_compl[:, filter_k]
    # need to sum to 1
    gridweights = gridweights / np.sum(gridweights)

    # set the random seed - mainly for testing
    if not None:
        np.random.seed(ranseed)

    # sample to get the indexes of the picked models
    indx = range(n_models)
    sim_indx = np.random.choice(indx, size=nsim, p=gridweights)

    # get the vega fluxes for the filters
    _, vega_flux, _ = Vega(source=vega_fname).getFlux(sedgrid.filters)

    # setup the output table
    ot = Table()
    qnames = list(sedgrid.keys())
    # simulated data
    for k, filter in enumerate(sedgrid.filters):
        colname = "%s_RATE" % filter.split(sep="_")[-1].upper()
        simflux_wbias = flux[sim_indx, k] + model_bias[sim_indx, k]
        simflux = np.random.normal(loc=simflux_wbias, scale=model_unc[sim_indx, k])
        ot[colname] = Column(simflux / vega_flux[k])
    # model parmaeters
    for qname in qnames:
        ot[qname] = Column(sedgrid[qname][sim_indx])

    return ot
