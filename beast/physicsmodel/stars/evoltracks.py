# use evolutionary tracks instead of isochrones as the basis of
# the stellar physicsgrid

import numpy as np
from scipy.interpolate import interp1d

from astropy.table import Table, QTable

from beast.config import __ROOT__, solar_metalicity


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

    def plot_tracks(self, ax, xval="logT", yval="logL", trackval=0, linestyle="-", color="k"):
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

        trackval : int, optional
            which set of tracks to plot

        linestyle : string
            matplotlib linestyle

        color : string
            matplotlib color
        """
        if xval not in self.data[trackval].keys():
            raise ValueError("xval choice not in data table")
        if yval not in self.data[trackval].keys():
            raise ValueError("yval choice not in data table")

        # get uniq M_ini values
        uvals, indices = np.unique(self.data[trackval]["M_ini"], return_inverse=True)
        for k, cval in enumerate(uvals):
            cindxs = np.where(k == indices)
            # ax.plot(
            #     self.data[xval][cindxs], self.data[yval][cindxs], linestyle=linestyle, color=color,
            # )
            ax.plot(
                self.data[trackval][xval][cindxs],
                self.data[trackval][yval][cindxs],
                "o",
                color=color,
                markersize=2,
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
            keys in grid to generate metrics for

        Returns
        -------
        metrics : dictionary
            each entry has an array with [min, max, median, mean] deltas
            for that grid parameters
        """
        # loop over eep values accumulating deltas
        dvals = {}
        for cname in check_keys:
            dvals[cname] = []

        metrics = []
        for cdata in self.data:
            uvals, indices = np.unique(cdata["M_ini"], return_inverse=True)
            for k, cval in enumerate(uvals):
                (cindxs,) = np.where(k == indices)
                for cname in check_keys:
                    dvals[cname] = np.concatenate(
                        (dvals[cname], np.absolute(np.diff(cdata[cname][cindxs])))
                    )

            # compute the metrics
            tmetrics = {}
            for cname in check_keys:
                tmetrics[cname] = np.array([
                    np.min(dvals[cname]),
                    np.max(dvals[cname]),
                    np.median(dvals[cname]),
                    np.mean(dvals[cname]),
                ])
            metrics.append(tmetrics)

        return metrics

    @staticmethod
    def regrid_one_met(
        one_track,
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
        uvals, indices = np.unique(one_track["M_ini"], return_inverse=True)
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
        for cname in one_track.keys():
            new_grid[cname] = np.array([])

        # loop over eep values and interopolate to new mass grid
        # along constant eep tracks
        uvals, indices = np.unique(one_track["eep"], return_inverse=True)
        for k, cval in enumerate(uvals):
            (cindxs,) = np.where(k == indices)
            cur_masses = one_track["M_ini"][cindxs]

            # only interpolate for masses defined for the current eep
            (new_gindxs,) = np.where(
                np.logical_and(
                    min(cur_masses) <= new_mass_vals, new_mass_vals <= max(cur_masses)
                )
            )

            for cname in one_track.keys():
                if cname == "eep":
                    vals = np.full((len(new_gindxs)), cval)
                elif cname == "M_ini":
                    vals = new_mass_vals[new_gindxs]
                else:
                    f = interp1d(cur_masses, one_track[cname][cindxs])
                    vals = f(new_mass_vals[new_gindxs])

                new_grid[cname] = np.concatenate((new_grid[cname], vals))

        # update the grid
        one_track = new_grid

        # setup the condensed grid
        if condense:
            new_grid = {}
            for cname in one_track.keys():
                new_grid[cname] = np.array([])

            # loop over each mass track and condense
            uvals, indices = np.unique(one_track["M_ini"], return_inverse=True)
            for k, cval in enumerate(uvals):
                (cindxs,) = np.where(k == indices)
                delta_logL = np.absolute(np.diff(one_track["logL"][cindxs]))
                delta_logT = np.absolute(np.diff(one_track["logT"][cindxs]))
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

                for cname in one_track.keys():
                    new_grid[cname] = np.concatenate(
                        (new_grid[cname], one_track[cname][cindxs][nindxs])
                    )

            # update the grid
            one_track = new_grid

        return one_track


class ETMist(EvolTracks):
    """
    MIST evolutionary tracks
    """

    def __init__(self):
        super().__init__()
        self.name = "MIST EvolTracks"
        self.source = "MIST"

        # fmt: off
        self.orig_FeH = np.array([-4.00, -3.50, -3.00, -2.50, -2.00, -1.75,
                                  -1.50, -1.25, -1.00, -0.75, -0.25, 0.0,
                                  0.25, 0.5])
        # fmt: on
        self.orig_files = [
            f"{__ROOT__}MIST/MIST_FeH{cstr:.2f}_vvcrit0.4.fits"
            for cstr in self.orig_FeH
        ]

        self.logmass_range = np.log10(np.array([0.1, 300.0]))
        self.z_range = (
            10 ** (np.array([min(self.orig_FeH), max(self.orig_FeH)]))
            * solar_metalicity
        )

    def load_orig_tables(self, FeH=None, filename=None):
        """
        Read the tracks from the original files

        Returns
        -------
        orig_tracks : astropy Table
            Table with evolutionary track info as columns for all metallicities
        """
        if filename is None:
            if isinstance(self.orig_files, list):
                files = self.orig_files
            else:
                files = [self.orig_files]
        else:
            files = filename

        itables = [QTable.read(cfile) for cfile in files]

        return itables

    def get_evoltracks(
        self, mass_info, metals, condense=False, logT_delta=0.05, logL_delta=0.05
    ):
        """
        Get the evolutionary tracks for the specified ages, initial masses,
        and metallicities.

        Parameters
        ----------
        mass_info : list
            info for mass range [logminmass, logmaxmass, deltalogmass]
        metals : list
            metallicities
        logL_delta, logT_delta : float
            deltas for the condensed grid
            default is 0.05 for both

        Returns
        -------
        type
            Description of returned object.
        """
        # determine the min/max metallilcities needed from original grid
        # FeH = np.log10(np.array(metals) / solar_metalicity)
        # 0.5 is make sure the orig FeH includes one metallicity beyond the values
        # gvals = (self.orig_FeH >= min(FeH) - 0.5) & (self.orig_FeH <= max(FeH) + 0.5)
        # print(self.orig_FeH)
        # print(self.orig_FeH[gvals])
        # exit()

        # get the as computed evolutionary tracks
        orig_data = self.load_orig_tables()
        print(orig_data[0])

        # interpolate for requested mass spacing
        print(mass_info)
        for k, ctrack in enumerate(orig_data):
            orig_data[k] = self.regrid_one_met(
                ctrack,
                logmass_range=mass_info[0:2],
                logmass_delta=mass_info[2],
                condense=condense,
                logT_delta=logT_delta,
                logL_delta=logL_delta,
            )
        print(len(orig_data[0]["M_act"]))
        exit()

        # then interpolate for metallicity spacing

        # condense grid if requested


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
