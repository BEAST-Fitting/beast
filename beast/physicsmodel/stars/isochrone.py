"""
Isochrone class

Intent to implement a generic module to manage isochrone mining from various
sources.
"""
import numpy as np
from numpy import interp
from numpy import log10
from scipy import interpolate
from astropy import units
import tables
from astropy.table import Table
from numpy.lib import recfunctions

# from beast.external.eztables import Table
# from beast.external.eztables.table import recfunctions
from beast.config import __ROOT__
from beast.physicsmodel.stars.ezpadova import parsec
from beast.physicsmodel.stars.ezmist import mist

__all__ = ["Isochrone", "padova2010", "pegase", "ezIsoch", "PadovaWeb", "MISTWeb"]


class Isochrone(object):
    def __init__(self, name="", *args, **kwargs):
        self.name = name

        # define axis labels for plotting
        self.alabels = {
            "logT": "log(Teff)",
            "logg": "log(g)",
            "logL": "log(L)",
            "logA": "log(age)",
            "stage": "evol phase",
            "M_act": "log(current mass)",
            "M_ini": "log(initial mass)",
        }

    def metalToFeH(self, metal):
        """
        Convert Z to [Fe/H] values

        For example:
           Zsun = 0.02 will give [Fe/H]sun = -4.33

        For reference:
           Z = [ 0.0004, 0.004, 0.008, 0.02, 0.05 ]
           [Fe/H] = [ -1.7  , -0.7 , -0.4 , 0   , 0.4  ]
        """
        return np.log10(metal / 0.02)

    def FeHtometal(self, feh):
        """
        Convert Z to [Fe/H] values
        """
        return 10 ** feh * 0.02

    def _get_isochrone(self, *args, **kwargs):
        """ Retrieve isochrone from the original source
            internal use to adapt any library
        """
        pass

    def _get_continuous_isochrone(self, *args, **kwargs):
        """ Return a resampled isochrone accounting for variations
            useful for continuous sampling
        """
        # define the maximum allowable difference between points
        dm = kwargs.pop("dm", 0.01)
        dt = kwargs.pop("dt", 0.01)
        dl = kwargs.pop("dl", 0.01)

        iso = self._get_isochrone(*args, **kwargs)
        logT, logg, logL, logM = (iso["logT"], iso["logg"], iso["logL"], iso["logM"])

        # compute vector of discrete derivaties for each quantity
        # and the final number of points
        npts = (np.abs(np.divide(np.diff(logM), dm))).astype(int)
        npts += (np.abs(np.divide(np.diff(logT), dt))).astype(int)
        npts += (np.abs(np.divide(np.diff(logL), dl))).astype(int)
        idx = np.hstack([[0], np.cumsum(npts + 1)])
        # set up vectors for storage
        ntot = (npts + 1).sum()
        newm = np.zeros(ntot, dtype=float)
        newdm = np.zeros(ntot, dtype=float)
        newt = np.zeros(ntot, dtype=float)
        newg = np.zeros(ntot, dtype=float)
        newl = np.zeros(ntot, dtype=float)

        for i in range(len(npts)):
            a, b = idx[i], idx[i] + npts[i] + 1
            if npts[i] > 0:
                # construct new 1d grids in each dimension, being careful
                #   about endpoints
                # append them to storage vectors
                newm[a:b] = np.linspace(
                    logM[i], logM[i + 1], npts[i] + 1, endpoint=False
                )
                newt[a:b] = np.linspace(
                    logT[i], logT[i + 1], npts[i] + 1, endpoint=False
                )
                newg[a:b] = np.linspace(
                    logg[i], logg[i + 1], npts[i] + 1, endpoint=False
                )
                newl[a:b] = np.linspace(
                    logL[i], logL[i + 1], npts[i] + 1, endpoint=False
                )
                newdm[a:b] = (
                    np.ones(npts[i] + 1) * (logM[i + 1] - logM[i]) / (npts[i] + 1)
                )
            else:
                # if the maximumum allowable difference is small,
                # then just store the good point
                newm[a] = logM[i]
                newt[a] = logT[i]
                newg[a] = logg[i]
                newl[a] = logL[i]
                newdm[a] = logM[i + 1] - logM[i]
        # tack on the last point on the grid, as the loop is one element short
        newm[-1] = logM[-1]
        newt[-1] = logT[-1]
        newg[-1] = logg[-1]
        newl[-1] = logL[-1]
        newdm[-1] = logM[-1] - logM[-2]

        table = Table(dict(logM=newm, logT=newt, logg=newg, logL=newl, dlogm=newdm))

        for k in list(iso.header.keys()):
            table.header[k] = iso.header[k]

        table.header["NAME"] = "Resampled " + table.header["NAME"]

        table.header["dlogT"] = dt
        table.header["dlogM"] = dm
        table.header["dlogg"] = dl

        return table

    def plot(self, ax, xval="logT", yval="logL", linestyle="-"):
        """
        Plot the isochrones with the input x, y choices

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

        # get uniq logA values
        uvals, indices = np.unique(self.data["logA"], return_inverse=True)
        for k, cval in enumerate(uvals):
            cindxs = np.where(k == indices)
            ax.plot(
                self.data[xval][cindxs], self.data[yval][cindxs], linestyle=linestyle
            )

        if xval in ["M_ini", "M_act"]:
            ax.set_xscale("log")
        if yval in ["M_ini", "M_act"]:
            ax.set_yscale("log")

        ax.set_xlabel(self.alabels[xval])
        ax.set_ylabel(self.alabels[yval])

        if xval == "logT":
            xmin, xmax = ax.get_xlim()
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

        for gparam in ["logA"]:
            uvals, indices = np.unique(self.data[gparam], return_inverse=True)
            for k, cval in enumerate(uvals):
                cindxs, = np.where(k == indices)
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


class padova2010(Isochrone):
    def __init__(self):
        super().__init__()
        self.name = "Padova 2010 (Marigo 2008 + Girardi 2010)"
        self.source = __ROOT__ + "/padova2010.iso.fits"
        self._load_table_(self.source)
        self.ages = 10 ** np.unique(self.data["logA"])
        self.Z = np.unique(self.data["Z"])

    def _load_table_(self, source):
        t = Table(self.source)
        data = {}
        for k in list(t.keys()):
            data[k] = t[k]
        # Alias columns
        data["logM"] = log10(np.asarray(data["M_ini"]))
        data["logg"] = np.asarray(data["logG"])
        data["logT"] = np.asarray(data["logTe"])
        data["logL"] = np.asarray(data["logL/Lo"])
        data["logA"] = np.asarray(data["log(age/yr)"])
        # clean columns
        data.pop("log(age/yr)")
        data.pop("M_ini")
        data.pop("logG")
        data.pop("logTe")
        data.pop("logL/Lo")

        self.data = Table(data, name="Isochrone from %s" % self.name)

    def _get_isochrone(self, age, metal=None, FeH=None, masses=None, *args, **kwargs):
        """ Retrieve isochrone from the original source
            internal use to adapt any library
        """
        # make sure unit is in years and then only give the value (no units)
        _age = int(units.Quantity(age, units.year).value)

        # if hasUnit(age):
        #    _age = int(age.to('yr').magnitude)
        # else:
        #    _age = int(age * inputUnit.to('yr').magnitude)

        assert (metal is not None) | (FeH is not None), "Need a chemical par. value."

        if (metal is not None) & (FeH is not None):
            print("Warning: both Z & [Fe/H] provided, ignoring [Fe/H].")

        if metal is None:
            metal = self.FeHtometal(FeH)

        assert metal in self.Z, "Metal %f not find in %s" % (metal, self.Z)

        data = {}
        t = self.data.selectWhere("*", "(Z == _z)", condvars={"_z": metal})
        if _age in self.ages:
            # no interpolation, isochrone already in the file
            t = t.selectWhere("*", "(logA == _age)", condvars={"_age": log10(_age)})
            for kn in list(t.keys()):
                data[kn] = np.asarray(t[kn])
        else:
            # interpolate between isochrones
            d = (self.ages - float(_age)) ** 2
            a1, a2 = self.ages[np.argsort(d)[:2]]
            # print "Warning: Interpolation between %d and %d Myr" % (a1, a2)
            r = np.log10(_age / a1) / np.log10(a2 / a1)

            t1 = t.selectWhere("*", "logA == _age", condvars={"_age": log10(a1)})
            t2 = t.selectWhere("*", "logA == _age", condvars={"_age": log10(a2)})

            stop = min(t1.nrows, t2.nrows)

            for kn in list(t1.keys()):
                y2 = t2[kn][:stop]
                y1 = t1[kn][:stop]
                data[kn] = y2 * r + y1 * (1.0 - r)
                del y1, y2

        # mass selection
        if masses is not None:
            # masses are expected in logM for interpolation
            if masses.max() > 2.3:
                _m = np.log10(masses)
            else:
                _m = masses
            data_logM = data["logM"][:]
            for kn in data:
                data[kn] = interp(_m, data_logM, data[kn])

        del t
        table = Table(data, name="Isochrone from %s" % self.name)
        table.header["metal"] = metal
        table.header["time"] = _age
        return table


class pegase(Isochrone):
    def __init__(self):
        super().__init__()
        self.name = "Pegase.2 (Fioc+1997)"
        self.source = __ROOT__ + "/pegase.iso.hd5"
        self.data = tables.openFile(self.source)
        self.ages = np.sort(
            np.asarray([k.attrs.time for k in self.data.root.Z02]) * 1e6
        )
        self.Z = np.asarray(
            [
                float("0." + k[1:])
                for k in self.data.root._g_listGroup(self.data.getNode("/"))[0]
            ]
        )

    def __getstate__(self):
        self.data.close()
        self.data = None
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__ = d
        self.data = tables.openFile(self.source)

    def __del__(self):
        if self.data is not None:
            self.data.close()

    def _get_isochrone(self, age, metal=None, FeH=None, masses=None, *args, **kwargs):
        """ Retrieve isochrone from the original source
            internal use to adapt any library
        """
        # make sure unit is in years and then only give the value (no units)
        _age = int(units.Quantity(age, units.year).value)

        #        if hasUnit(age):
        #            _age = int(age.to('Myr').magnitude)
        #        else:
        #            _age = int(age * inputUnit.to('Myr').magnitude)

        assert (metal is not None) | (FeH is not None), "Need a chemical par. value."

        if (metal is not None) & (FeH is not None):
            print("Warning: both Z & [Fe/H] provided, ignoring [Fe/H].")

        if metal is None:
            metal = self.FeHtometal(FeH)

        assert metal in self.Z, "Metal %f not find in %s" % (metal, self.Z)
        # node = self.data.getNode('/Z' + str(metal)[2:])

        data = {}
        if age in self.ages:
            # no interpolation, isochrone already in the file
            t = self.data.getNode("/Z" + str(metal)[2:] + "/a" + str(_age))
            for kn in t.colnames:
                data[kn] = t.col(kn)
        else:
            # interpolate between isochrones
            d = (self.ages - float(age)) ** 2
            a1, a2 = np.sort(self.ages[np.argsort(d)[:2]] * 1e-6)
            # print "Warning: Interpolation between %d and %d Myr" % (a1, a2)
            r = np.log10(_age / a1) / np.log10(a2 / a1)

            t1 = self.data.getNode("/Z" + str(metal)[2:] + "/a" + str(int(a1)))
            t2 = self.data.getNode("/Z" + str(metal)[2:] + "/a" + str(int(a2)))

            stop = min(t1.nrows, t2.nrows)

            for kn in t1.colnames:
                y2 = t2.col(kn)[:stop]
                y1 = t1.col(kn)[:stop]
                data[kn] = y2 * r + y1 * (1.0 - r)
                del y1, y2

        # mass selection
        if masses is not None:
            # masses are expected in logM for interpolation
            if masses.max() > 2.3:
                _m = np.log10(masses)
            else:
                _m = masses
            data_logM = data["logM"][:]
            for kn in data:
                data[kn] = interp(_m, data_logM, data[kn])

        table = Table(data, name="Isochrone from %s" % self.name)
        table.header["metal"] = metal
        table.header["time"] = _age * 1e6
        return table


class ezIsoch(Isochrone):
    """ Trying to make something that is easy to manipulate
    This class is basically a proxy to a table (whatever format works best)
    and tries to keep things coherent.
    """

    def __init__(self, source, interp=False):
        super().__init__()
        self.name = "<auto>"
        self.source = source
        self._load_table_(self.source)
        # round because of precision noise
        self.logages = np.unique(np.round(self.data["logA"], 6))
        self.ages = np.round(10 ** self.logages)
        self.Z = np.unique(np.round(self.data["Z"], 6))
        self.interpolation(interp)

#    def selectWhere(self, *args, **kwargs):
#        return self.data.selectWhere(*args, **kwargs)

    def interpolation(self, b=None):
        if b is not None:
            if hasattr(self, "interp"):
                print("Do not use interpolation yet, at your own risks!!")
            self.interp = bool(b)
        else:
            return self.interp

    def _load_table_(self, source):
        tdata = Table.read(self.source, format="csv", comment="#")
        self.data = tdata[np.isfinite(tdata["logA"])]

    def __getitem__(self, key):
        return self.data[key]

    def _get_t_isochrone(self, age, metal=None, FeH=None, masses=None, *args, **kwargs):
        """ Retrieve isochrone from the original source
            internal use to adapt any library
        """
        # make sure unit is in years and then only give the value (no units)
        _age = int(units.Quantity(age, units.year).value)

        #        if hasUnit(age):
        #            _age = int(age.to('yr').magnitude)
        #        else:
        #            _age = int(age * inputUnit.to('yr').magnitude)

        _logA = np.log10(_age)

        assert (metal is not None) | (FeH is not None), "Need a chemical par. value."

        if (metal is not None) & (FeH is not None):
            print("Warning: both Z & [Fe/H] provided, ignoring [Fe/H].")

        if metal is None:
            metal = self.FeHtometal(FeH)

        if self.interpolation():
            # Do the actual nd interpolation

            # Maybe already exists?
            if (metal in self.Z) & (_age in self.ages):
                t = self.selectWhere(
                    "*",
                    "(round(Z, 6) == {0}) & (round(logA, 6) == {1})".format(
                        metal, _logA
                    ),
                )
                if t.nrows > 0:
                    return t
            # apparently not
            # find 2 closest metal values
            ca1 = self.ages <= _age
            ca2 = self.ages > _age
            cz1 = self.Z <= metal
            cz2 = self.Z > metal
            if metal in self.Z:
                # perfect match in metal, need to find ages
                if _age in self.ages:
                    return self.selectWhere(
                        "*",
                        "(round(Z, 6) == {0}) & (round(logA, 6) == {1})".format(
                            metal, _logA
                        ),
                    )
                elif (True in ca1) & (True in ca2):
                    # bracket on _age: closest values
                    a1, a2 = (
                        np.log10(max(self.ages[ca1])),
                        np.log10(min(self.ages[ca2])),
                    )
                    iso = self.selectWhere(
                        "*",
                        "(Z == 0.02) & ( (abs(logA - {0}) < 1e-4) | (abs(logA - {1}) < 1e-4 )  )".format(
                            a1, a2
                        ),
                    )
                    if masses is None:
                        _logM = np.unique(iso["logM"])
                    else:
                        _logM = masses

                    # define interpolator
                    points = np.array([self[k] for k in "logA logM Z".split()]).T
                    values = np.array([self[k] for k in list(self.data.keys())]).T
                    _ifunc = interpolate.LinearNDInterpolator(points, values)

                    pts = np.array([(_logA, logMk, metal) for logMk in _logM])
                    r = _ifunc(pts)
                    return Table(r)
                else:
                    raise Exception("Age not covered by the isochrones")
            elif (True in cz1) & (True in cz2):
                # need to find closest Z
                pass
            return
        else:
            # find the closest match
            _Z = self.Z[((metal - self.Z) ** 2).argmin()]
            # _logA = np.log10(self.ages[((_age - self.ages) ** 2).argmin()])
            _logA = self.logages[((np.log10(_age) - self.logages) ** 2).argmin()]
            tab = self.data.selectWhere(
                "*", "(round(Z, 6) == {0}) & (round(logA,6) == {1})".format(_Z, _logA)
            )
            # mass selection
            if masses is not None:
                # masses are expected in logM for interpolation
                # if masses.max() > 2.3:
                #    _m = np.log10(masses)
                # else:
                _m = masses
                data_logM = tab["logM"][:]
                # refuse extrapolation!
                # ind = np.where(_m <= max(data_logM))
                data = {}
                for kn in list(tab.keys()):
                    data[kn] = interp(_m, data_logM, tab[kn], left=np.nan, right=np.nan)
                return Table(data)


class PadovaWeb(Isochrone):
    def __init__(
        self,
        Zref=None,
        modeltype="parsec12s_r14",
        filterPMS=False,
        filterBad=False,
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.name = "Padova CMD isochrones"
        if Zref is None:
            if modeltype.startswith("parsec"):
                Zref = 0.0152
            else:
                Zref = 0.019
        self.Zref = Zref
        self.modeltype = modeltype
        self.filterPMS = filterPMS
        self.filterBad = filterBad

    def _get_isochrone(self, age, metal=None, FeH=None, *args, **kwargs):
        """ Retrieve isochrone from the original source
            internal use to adapt any library
        """
        # make sure unit is in years and then only give the value (no units)
        _age = int(units.Quantity(age, units.year).value)

        #        if hasUnit(age):
        #            _age = int(age.to('yr').magnitude)
        #        else:
        #            _age = int(age * inputUnit.to('yr').magnitude)

        assert (metal is not None) | (FeH is not None), "Need a chemical par. value."

        if (metal is not None) & (FeH is not None):
            print("Warning: both Z & [Fe/H] provided, ignoring [Fe/H].")

        if metal is None:
            metal = self.FeHtometal(FeH)

        iso_table = parsec.get_one_isochrone(
            _age, metal, ret_table=True, model=self.modeltype
        )
        iso_table = self._clean_cols(iso_table)
        iso_table = self._filter_iso_points(
            iso_table, filterPMS=self.filterPMS, filterBad=self.filterBad
        )

        return iso_table

    def _clean_cols(self, iso_table):
        """clean column names, remove unnecessary columns"""
        # Rename Columns
        if self.modeltype == "parsec12s_r14":
            # PARSEC+COLIBRI Column Names
            iso_table.add_column("logA", np.log10(iso_table["Age"][:]))
            iso_table.add_column("logT", iso_table["logTe"][:])
            iso_table.add_column("M_ini", iso_table["Mini"][:])
            iso_table.add_column("M_act", iso_table["Mass"][:])
            iso_table.add_column("stage", iso_table["label"][:])
            iso_table.remove_columns(["Age", "logTe", "Mini", "Mass", "label"])
            # Remove age-specific Z, rename Zini as Z
            iso_table.remove_columns(["Z"])
            iso_table.add_column("Z", iso_table["Zini"][:])
            iso_table.remove_columns(["Zini"])

        else:
            # Padova (Girardi10, Marigo08, etc), Old PARSEC Column Names
            iso_table.add_column("logA", iso_table["logageyr"][:])
            iso_table.add_column("logL", iso_table["logLLo"][:])
            iso_table.add_column("logT", iso_table["logTe"][:])
            iso_table.add_column("logg", iso_table["logG"][:])
            iso_table.remove_columns(["logageyr", "logLLo", "logTe", "logG"])

        # Remove phot columns and unnecessary properties
        filternames = "U UX B BX V R I J H K L M".split()
        theorycols = [
            "C/O",
            "M_hec",
            "int_IMF",
            "period",
            "pmode",
            "CO",
            "C_O",
            "period0",
            "period1",
            "McoreTP",
            "tau1m",
        ]
        # removing mass loss outputs
        theorycols += ["logMdot", "Mloss"]
        abundcols = "X Y Xc Xn Xo Cexcess".split()
        drop = theorycols + abundcols + filternames + [s + "mag" for s in filternames]
        # make sure columns exist
        iso_table.remove_columns([x for x in drop if x in iso_table])

        # polish the header
        iso_table.setUnit("logA", "yr")
        iso_table.setComment("logA", "Age")
        iso_table.setUnit("logT", "K")
        iso_table.setComment("logT", "Effective temperature")
        iso_table.setUnit("logL", "Lsun")
        iso_table.setComment("logL", "Luminosity")
        iso_table.setUnit("M_ini", "Msun")
        iso_table.setComment("M_ini", "Initial Mass")
        iso_table.setUnit("M_act", "Msun")
        iso_table.setComment("M_act", "Current Mass, M(t)")
        iso_table.setUnit("logg", "cm/s**2")
        iso_table.setComment("logg", "Surface gravity")
        iso_table.setComment("stage", "Evolutionary Stage")
        iso_table.setComment("Z", "Metallicity")
        # iso_table.setUnit('logMdot', 'Msun/yr')
        # iso_table.setComment('logMdot', 'Mass loss')

        return iso_table

    def _filter_iso_points(self, iso_table, filterPMS=False, filterBad=False):
        """ Filter bad points and PMS points
            Bad points known to affect pre-PARSEC isochronesself.
            Selection is an empirical definition.
        """
        # Filter pre-ms stars
        if filterPMS:
            cond = "~((M_ini < 12.) & (stage == 0))"
            iso_table = iso_table.selectWhere("*", cond)

        # Filter bad points for pre-PARSEC, PadovaCMDVersion < 2.7 isochrones
        if filterBad:
            if not self.modeltype.startswith("parsec"):
                cond = "~((logL > 3.) & (M_act < 1.) & (log10(M_ini / M_act) > 0.1))"
                iso_table = iso_table.selectWhere("*", cond)
            else:
                print("No bad point filtering for PARSEC models.")

        return iso_table

    def _get_t_isochrones(self, logtmin, logtmax, dlogt, Z=0.0152):
        """ Generate a proper table directly from the PADOVA website

        Parameters
        ----------
        logtmin: float
            log-age min (age in yr)

        logtmax: float
            log-age max (age in yr)

        dlogt: float
            log-age step to request

        Z: float or sequence
            single value of list of values of metalicity Z

        returns
        -------
        tab: eztable.Table
            the table of isochrones
        """
        if not hasattr(Z, "__iter__"):
            iso_table = parsec.get_t_isochrones(
                max(6.0, logtmin), min(10.13, logtmax), dlogt, Z, model=self.modeltype
            )
            iso_table.header["NAME"] = "PadovaCMD Isochrones: " + self.modeltype
            if "Z" not in iso_table:
                iso_table.add_column("Z", np.ones(iso_table.nrows) * Z)

            # rename cols, remove phot and other unnecessary cols
            iso_table = self._clean_cols(iso_table)

            # filter iso data: pre-ms and bad points
            iso_table = self._filter_iso_points(
                iso_table, filterPMS=self.filterPMS, filterBad=self.filterBad
            )

        else:
            iso_table = self._get_t_isochrones(logtmin, logtmax, dlogt, Z[0])
            iso_table.header["NAME"] = "PadovaCMD Isochrones: " + self.modeltype

            if len(Z) > 1:
                more = [
                    self._get_t_isochrones(logtmin, logtmax, dlogt, Zk).data
                    for Zk in Z[1:]
                ]
                iso_table.data = recfunctions.stack_arrays(
                    [iso_table.data] + more, usemask=False, asrecarray=True
                )

        return iso_table


class MISTWeb(Isochrone):
    def __init__(self, Zref=0.0142, rotation="vvcrit0.0", *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "MESA/MIST isochrones"
        self.Zref = Zref
        self.rotation = rotation

    def _get_isochrone(self, age, metal=None, FeH=None, *args, **kwargs):
        """ Retrieve isochrone from the original source
            internal use to adapt any library
        """
        # make sure unit is in years and then only give the value (no units)
        _age = int(units.Quantity(age, units.year).value)

        #        if hasUnit(age):
        #            _age = int(age.to('yr').magnitude)
        #        else:
        #            _age = int(age * inputUnit.to('yr').magnitude)

        assert (metal is not None) | (FeH is not None), "Need a chemical par. value."

        if (metal is not None) & (FeH is not None):
            print("Warning: both Z & [Fe/H] provided, ignoring [Fe/H].")

        if FeH is None:
            FeH = self.metaltoFeH(metal)

        iso_table = mist.get_one_isochrone(
            _age, FeH, v_div_vcrit=self.rotation, age_scale="log10", ret_table=True
        )
        iso_table = self._clean_cols(iso_table)

        return iso_table

    def _clean_cols(self, iso_table):
        """clean column names, remove unnecessary columns"""
        # Rename Columns
        iso_table.add_column("logA", iso_table["log10_isochrone_age_yr"][:])
        iso_table.add_column("logT", iso_table["log_Teff"][:])
        iso_table.add_column("logL", iso_table["log_L"][:])
        iso_table.add_column("M_ini", iso_table["initial_mass"][:])
        iso_table.add_column("M_act", iso_table["star_mass"][:])
        iso_table.add_column("logg", iso_table["log_g"][:])
        iso_table.add_column("stage", iso_table["phase"][:])
        iso_table.remove_columns(
            [
                "log10_isochrone_age_yr",
                "log_Teff",
                "log_L",
                "log_g",
                "initial_mass",
                "star_mass",
                "phase",
            ]
        )

        # Remove phot columns and unnecessary properties
        extracol1 = "star_mdot he_core_mass c_core_mass log_LH log_LHe log_R".split()
        extracol2 = "log_center_T log_center_Rho center_gamma center_h1 center_he4 center_c12".split()
        extracol3 = "surface_h1 surface_he3 surface_he4 surface_c12 surface_o16".split()
        drop = extracol1 + extracol2 + extracol3
        # make sure columns exist
        iso_table.remove_columns([x for x in drop if x in iso_table])

        # polish the header
        iso_table.setUnit("logA", "yr")
        iso_table.setComment("logA", "Age")
        iso_table.setUnit("logT", "K")
        iso_table.setComment("logT", "Effective temperature")
        iso_table.setUnit("logL", "Lsun")
        iso_table.setComment("logL", "Luminosity")
        iso_table.setUnit("M_ini", "Msun")
        iso_table.setComment("M_ini", "Initial Mass")
        iso_table.setUnit("M_act", "Msun")
        iso_table.setComment("M_act", "Current Mass, M(t)")
        iso_table.setUnit("logg", "cm/s**2")
        iso_table.setComment("logg", "Surface gravity")
        iso_table.setComment("stage", "Evolutionary Stage")
        iso_table.setComment("Z", "Metallicity")
        # iso_table.setUnit('logMdot', 'Msun/yr')
        # iso_table.setComment('logMdot', 'Mass loss')

        return iso_table

    def _get_t_isochrones(self, logtmin, logtmax, dlogt, Z=0.0142):
        """ Generate a proper table directly from the PADOVA website

        Parameters
        ----------
        logtmin: float
            log-age min (age in yr)

        logtmax: float
            log-age max (age in yr)

        dlogt: float
            log-age step to request

        Z: float or sequence
            single value of list of values of metalicity Z

        returns
        -------
        tab: eztable.Table
            the table of isochrones
        """

        if not hasattr(Z, "__iter__"):
            iso_table = mist.get_t_isochrones(
                max(5.0, logtmin),
                min(10.13, logtmax),
                dlogt,
                v_div_vcrit=self.rotation,
                FeH_value=np.log10(Z / self.Zref),
            )
            iso_table.header["NAME"] = "MESA/MIST Isochrones"
            if "Z" not in iso_table:
                iso_table.add_column("Z", np.ones(iso_table.nrows) * Z)

            # rename cols, remove phot and other unnecessary cols
            iso_table = self._clean_cols(iso_table)

        else:
            iso_table = self._get_t_isochrones(logtmin, logtmax, dlogt, Z[0])
            iso_table.header["NAME"] = "MESA/MIST Isochrones"

            if len(Z) > 1:
                more = [
                    self._get_t_isochrones(logtmin, logtmax, dlogt, Zk).data
                    for Zk in Z[1:]
                ]
                iso_table.data = recfunctions.stack_arrays(
                    [iso_table.data] + more, usemask=False, asrecarray=True
                )

        return iso_table
