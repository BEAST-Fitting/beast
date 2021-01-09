# use evolutionary tracks instead of isochrones as the basis of
# the stellar physicsgrid

import numpy as np
from scipy.interpolate import interp1d

from astropy.table import Table, vstack

from beast.config import __ROOT__


__all__ = ["EvolTracks", "ETParsec", "ETMist"]


class EvolTracks(object):
    """
    Stores the evolutionary tracks interpolated to input mass and
    metallicity grids with the age grid set by the evolutionary track
    calculations (may wish to revisit and allow for condensation).

    Attributes
    ----------
    logmass : log10 of the masses (solar masses)

    Z : numpy array of floats
        metallicites of tracks

    data: table
        columns include:
        M_ini - initial mass [Msun]
        M_act - actual mass [Msun]
        Z - metallicity [??]
        logL - log luminosity [Lsun]
        logg - log surface gravity [cm/s^2]
        logA - log age [years]
        logT - log surface effective temperature [K]
        stage - evolutionary stage [index]
    """

    def __init__(self):
        self.name = "<auto>"

        # define axis labels for plotting
        self.alabels = {
            "logT": "log(Teff)",
            "logg": "log(g)",
            "logL": "log(L)",
            "logA": "log(age)",
            "phase": "evol phase",
            "M_act": "log(current mass)",
            "M_ini": "log(initial mass)",
            "eep": "EEP",
        }

    def load_orig_tables(self, source):
        """
        Read the tracks from  the original files (stub)
        """
        print("not implemented")

    def plot_tracks(self, ax, xval="logT", yval="logL", linestyle="-"):
        """
        Plot the tracks with the input x, y choices

        Parameters
        ----------
        ax : matplotlib axis
            where to plot

        xval : str, optional
            what data for x

        xval : str, optional
            what data for y

        linestyle : string
            matplotlib linestyle
        """
        if xval not in self.data.keys():
            raise ValueError("xval choice not in data table")
        if yval not in self.data.keys():
            raise ValueError("yval choice not in data table")

        # get uniq M_ini values
        uvals, indices = np.unique(self.data["M_ini"], return_inverse=True)
        for k, cval in enumerate(uvals):
            cindxs = np.where(k == indices)
            ax.plot(
                self.data[xval][cindxs], self.data[yval][cindxs], linestyle=linestyle
            )

        ax.set_xlabel(self.alabels[xval])
        ax.set_ylabel(self.alabels[yval])

        if xval == "logT":
            xmin, xmax = ax.get_xlim()
            if xmin < xmax:
                ax.set_xlim(xmax, xmin)

    def grid_metrics(self, check_keys=["logL", "logT", "logg"]):
        """
        Compute metrics of the grid
        Primarily to determine how well parameter space is covered

        Parameters
        ----------
        check_keys : string array
            keys in grid to generage metrics for

        Returns
        -------
        metrics : dictonary
            each entry has an array with [min, max, median, mean] deltas
            for that grid paramters
        """
        # loop over eep values accumulating deltas
        dvals = {}
        for cname in check_keys:
            dvals[cname] = []

        for gparam in ["M_ini"]:
            uvals, indices = np.unique(self.data[gparam], return_inverse=True)
            for k, cval in enumerate(uvals):
                (cindxs,) = np.where(k == indices)
                for cname in check_keys:
                    dvals[cname] = np.concatenate(
                        (dvals[cname], np.absolute(np.diff(self.data[cname][cindxs])))
                    )

        # compute the metrics
        metrics = {}
        for cname in check_keys:
            metrics[cname] = [
                np.min(dvals[cname]),
                np.max(dvals[cname]),
                np.median(dvals[cname]),
                np.mean(dvals[cname]),
            ]

        return metrics

    def regrid_one_met(
        self,
        logmass_range=[-1.0, 2.0],
        logmass_delta=0.05,
        condense=False,
        logL_delta=0.05,
        logT_delta=0.05,
    ):
        """
        Interpolate a set of evolutionary tracks for a single metallicity
        to a uniform grid in log(initial mass) and variable grid in stellar age.
        Use Equivalent Evolutionary Points (EEPs) values to do the mass
        interpolation.  EEPs are provide as part of the evolutionary tracks.

        Parameters
        ----------
        logmass_range : (float, float)
            range of new mass grid
            default is -1 to 2 (0.1 to 100 M_sun)

        condense : boolean
            set to condense grid using delta values
            default = False

        logmass_delta, logL_delta, logT_delta : float
            deltas for the condensed grid
            default is 0.05 for all 3
        """
        # get the unique mass values
        uvals, indices = np.unique(self.data["M_ini"], return_inverse=True)
        # print(10**uvals)

        # ensure there are evolutionary tracks spanning
        #   the min/max of the new grid --> no extrapolation
        new_min_mass = max([min(uvals), logmass_range[0]])
        new_max_mass = min([max(uvals), logmass_range[1]])

        n_new_masses = int((new_max_mass - new_min_mass) / logmass_delta) + 1
        new_mass_vals = np.linspace(new_min_mass, new_max_mass, num=n_new_masses)
        # print(n_new_masses, len(uvals))
        # print(10**new_mass_vals)

        # setup the new grid
        new_grid = {}
        for cname in self.data.keys():
            new_grid[cname] = np.array([])

        # loop over eep values and interopolate to new mass grid
        # along constant eep tracks
        uvals, indices = np.unique(self.data["eep"], return_inverse=True)
        for k, cval in enumerate(uvals):
            (cindxs,) = np.where(k == indices)
            cur_masses = self.data["M_ini"][cindxs]

            # only interpolate for masses defined for the current eep
            (new_gindxs,) = np.where(
                np.logical_and(
                    min(cur_masses) <= new_mass_vals, new_mass_vals <= max(cur_masses)
                )
            )

            for cname in self.data.keys():
                if cname == "eep":
                    vals = np.full((len(new_gindxs)), cval)
                elif cname == "M_ini":
                    vals = new_mass_vals[new_gindxs]
                else:
                    f = interp1d(cur_masses, self.data[cname][cindxs])
                    vals = f(new_mass_vals[new_gindxs])

                new_grid[cname] = np.concatenate((new_grid[cname], vals))

        # update the grid
        self.data = new_grid

        # setup the condensed grid
        if condense:
            new_grid = {}
            for cname in self.data.keys():
                new_grid[cname] = np.array([])

            # loop over each mass track and condense
            uvals, indices = np.unique(self.data["M_ini"], return_inverse=True)
            for k, cval in enumerate(uvals):
                (cindxs,) = np.where(k == indices)
                delta_logL = np.absolute(np.diff(self.data["logL"][cindxs]))
                delta_logT = np.absolute(np.diff(self.data["logT"][cindxs]))
                nindxs = [0]
                cdelt_logL = 0.0
                cdelt_logT = 0.0
                for i in range(len(delta_logL)):
                    cdelt_logL += delta_logL[i]
                    cdelt_logT += delta_logT[i]
                    if (cdelt_logL > logL_delta) or (cdelt_logT > logT_delta):
                        nindxs.append(i)
                        cdelt_logL = delta_logL[i]
                        cdelt_logT = delta_logT[i]

                if not max(nindxs) == len(delta_logL) - 1:
                    nindxs.append(len(delta_logL) - 1)

                for cname in self.data.keys():
                    new_grid[cname] = np.concatenate(
                        (new_grid[cname], self.data[cname][cindxs][nindxs])
                    )

            # update the grid
            self.data = new_grid


class ETMist(EvolTracks):
    """
    MIST evolutionary tracks
    """

    def __init__(self):
        super().__init__()
        self.source = "MIST"

        self.orig_FeH = [0.0, 0.25, 0.5]
        self.orig_files = [
            f"{__ROOT__}MIST/MIST_FeH{cstr:.2f}.fits" for cstr in self.orig_FeH
        ]

    def get_orig_tables(self):
        """
        Read the tracks from the original files

        Returns
        -------
        orig_tracks : astropy Table
            Table with evolutionary track info as columns for all metallicities
        """
        if isinstance(self.orig_files, list):
            files = self.orig_files
        else:
            files = [self.orig_files]

        itables = [Table.read(cfile) for cfile in files]
        if len(itables) > 1:
            orig_tracks = vstack(itables)
        else:
            orig_tracks = itables[0]

        return orig_tracks

    def get_evoltracks(self, masses, metals=None, FeHs=None):
        """
        Get the evolutionary tracks for the specified ages, initial masses,
        and metallicities.

        Parameters
        ----------
        masses : list
            Initial masses for grid
        metal, FeH : list
            At least one needs to be set for the grid metallicities

        Returns
        -------
        type
            Description of returned object.
        """
        orig_tracks = self.get_orig_tracks()

        # first interpolate for mass spacing

        # then interpolate for metallicity spacing


class ETParsec(EvolTracks):
    """
    EvolTracks specific to PARSEC calculations
    """

    source = "PARSEC"

    def load_orig_tables(self, files):
        """
        Read the tracks from  the original files

        files : str
            file or files with the evolutionary track calculations
            often each file is for a single initial mass
        """
        if not type(files) is list:
            files = [files]

        mass_act = np.array([])
        logL = np.array([])
        logT = np.array([])
        for cfile in files:
            a = Table.read(cfile, format="ascii")
            mass_act = np.concatenate((mass_act, a["MASS"].data))
            logL = np.concatenate((logL, a["LOG_L"].data))
            logT = np.concatenate((logT, a["LOG_TE"].data))

        self.data = {}
        self.data["M_act"] = np.array(mass_act)
        self.data["logL"] = np.array(logL)
        self.data["logT"] = np.array(logT)
