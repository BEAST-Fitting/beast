"""
Extinction Curves
"""
import numpy as np
from scipy import interpolate

from astropy import units

import dust_extinction.parameter_averages as dustext_par
import dust_extinction.averages as dustext_avg
from dust_extinction.helpers import _test_valid_x_range

from beast.config import __ROOT__
import beast.physicsmodel.dust.extinction_extension as dustext_extend

__all__ = [
    "ExtinctionLaw",
    "Cardelli89",
    "Fitzpatrick99",
    "Gordon03_SMCBar",
    "Gordon16_RvFALaw",
    "Generalized_RvFALaw",
    "Generalized_DustExt",
]

libdir = __ROOT__


class ExtinctionLaw(object):
    """
    Extinction Law Template Class

    Attributes
    ----------
    name : string
       Name identifying the extinction law
    """

    def __init__(self):
        self.name = "None"

    def function(self, lamb, *arg, **kwargs):
        """
        function to provide extinction curve given wavelength(s)
        """
        raise NotImplementedError

    def __call__(self, *args, **kwargs):
        return self.function(*args, **kwargs)

    def isvalid(self, *args, **kwargs):
        """
        Check if the current arguments are in the validity domain of the law
        Must be redefined if any restriction applies to the law
        """
        return True


class Cardelli89(ExtinctionLaw):
    """
    Cardelli89 Milky Way R(V) dependent Extinction Law

    From Cardelli, Clayton, and Mathis (1989, ApJ, 345, 245).

    Fitzpatrick99 should be used instead except for historical purposes
    as this newer law is based on 10x more observations and a better
    treatment of the optical/NIR photometry based portion of the curves.
    """

    def __init__(self):
        super().__init__()
        self.name = "Cardelli89"
        self.x_range = [0.3, 10.0]  # inverse microns

    def function(self, lamb, Av=1.0, Rv=3.1, Alambda=True, **kwargs):
        """
        Cardelli89 extinction Law

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which to evaluate the law.

        Av: float
            desired A(V) (default: 1.0)

        Rv: float
            desired R(V) (default: 3.1)

        Alambda: bool
            if set returns +2.5*1./log(10.)*tau, tau otherwise

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation as a function of wavelength
            depending on Alambda option +2.5*1./log(10.)*tau,  or tau
        """
        # ensure the units are in angstrom
        _lamb = units.Quantity(lamb, units.angstrom).value

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([lamb])
        else:
            _lamb = lamb[:]

        # convert to wavenumbers
        x = 1.0e4 / _lamb

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, self.x_range, self.name)

        # init variables
        a = np.zeros(np.size(x))
        b = np.zeros(np.size(x))
        # Infrared (Eq 2a,2b)
        ind = np.where((x >= 0.3) & (x < 1.1))
        a[ind] = 0.574 * x[ind] ** 1.61
        b[ind] = -0.527 * x[ind] ** 1.61
        # Optical & Near IR
        # Eq 3a, 3b
        ind = np.where((x >= 1.1) & (x < 3.3))
        y = x[ind] - 1.82
        a[ind] = (
            1.0
            + 0.17699 * y
            - 0.50447 * y ** 2
            - 0.02427 * y ** 3
            + 0.72085 * y ** 4
            + 0.01979 * y ** 5
            - 0.77530 * y ** 6
            + 0.32999 * y ** 7
        )
        b[ind] = (
            1.41338 * y
            + 2.28305 * y ** 2
            + 1.07233 * y ** 3
            - 5.38434 * y ** 4
            - 0.62251 * y ** 5
            + 5.30260 * y ** 6
            - 2.09002 * y ** 7
        )
        # UV
        # Eq 4a, 4b
        ind = np.where((x >= 3.3) & (x <= 8.0))
        a[ind] = 1.752 - 0.316 * x[ind] - 0.104 / ((x[ind] - 4.67) ** 2 + 0.341)
        b[ind] = -3.090 + 1.825 * x[ind] + 1.206 / ((x[ind] - 4.62) ** 2 + 0.263)

        ind = np.where((x >= 5.9) & (x <= 8.0))
        y = x[ind] - 5.9
        Fa = -0.04473 * y ** 2 - 0.009779 * y ** 3
        Fb = 0.21300 * y ** 2 + 0.120700 * y ** 3
        a[ind] = a[ind] + Fa
        b[ind] = b[ind] + Fb
        # Far UV
        # Eq 5a, 5b
        ind = np.where((x > 8.0) & (x <= 10.0))
        # Fa = Fb = 0
        y = x[ind] - 8.0
        a[ind] = -1.073 - 0.628 * y + 0.137 * y ** 2 - 0.070 * y ** 3
        b[ind] = 13.670 + 4.257 * y - 0.420 * y ** 2 + 0.374 * y ** 3

        # Case of -values x out of range [0.289,10.0]
        ind = np.where((x > 10.0) | (x < 0.3))
        a[ind] = 0.0
        b[ind] = 0.0

        # Return Extinction vector
        # Eq 1
        if Alambda:
            return (a + b / Rv) * Av
        else:
            # return( 1./(2.5 * 1. / np.log(10.)) * ( a + b / Rv ) * Av)
            return 0.4 * np.log(10.0) * (a + b / Rv) * Av


class Fitzpatrick99(ExtinctionLaw):
    """
    Fitzpatrick99 Milky Way R(V) dependent Extinction Law

    From Fitzpatrick (1999, PASP, 111, 63).

    R(V) dependent extinction curve that explicitly deals with optical/NIR
    extinction being measured from broad/medium band photometry.
    Code based on fm_unred.pro from the IDL astronomy library

    Extended to below 912 A using Draine et al. dust grain models.
    """

    def __init__(self):
        super().__init__()
        self.name = "Fitzpatrick99"
        self.x_range = [0.3, 10.0]

    def function(self, lamb, Av=1, Rv=3.1, Alambda=True, **kwargs):
        """
        Fitzpatrick99 extinction Law

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which to evaluate the law.

        Av: float
            desired A(V) (default 1.0)

        Rv: float
            desired R(V) (default 3.1)

        Alambda: bool
            if set returns +2.5*1./log(10.)*tau, tau otherwise

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation as a function of wavelength
            depending on Alambda option +2.5*1./log(10.)*tau,  or tau
        """
        # ensure the units are in angstrom
        _lamb = units.Quantity(lamb, units.angstrom).value

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([lamb])
        else:
            _lamb = lamb[:]

        # convert to wavenumbers
        x = 1.0e4 / _lamb

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, self.x_range, self.name)

        # initialize values
        c2 = -0.824 + 4.717 / Rv
        c1 = 2.030 - 3.007 * c2
        c3 = 3.23
        c4 = 0.41
        x0 = 4.596
        gamma = 0.99

        k = np.zeros(np.size(x))

        # compute the UV portion of A(lambda)/E(B-V)
        xcutuv = 10000.0 / 2700.0
        xspluv = 10000.0 / np.array([2700.0, 2600.0])
        ind = np.where(x >= xcutuv)

        if np.size(ind) > 0:
            k[ind] = (
                c1
                + (c2 * x[ind])
                + c3
                * ((x[ind]) ** 2)
                / (((x[ind]) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((x[ind]) ** 2))
            )
            yspluv = (
                c1
                + (c2 * xspluv)
                + c3
                * ((xspluv) ** 2)
                / (((xspluv) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((xspluv) ** 2))
            )

            # FUV portion
            fuvind = np.where(x >= 5.9)
            k[fuvind] += c4 * (
                0.5392 * ((x[fuvind] - 5.9) ** 2)
                + 0.05644 * ((x[fuvind] - 5.9) ** 3)
            )

            k[ind] += Rv
            yspluv += Rv

        # Optical/NIR portion

        ind = np.where(x < xcutuv)
        if np.size(ind) > 0:
            xsplopir = np.zeros(7)
            xsplopir[0] = 0.0
            xsplopir[1:7] = 10000.0 / np.array(
                [26500.0, 12200.0, 6000.0, 5470.0, 4670.0, 4110.0]
            )

            ysplopir = np.zeros(7)
            ysplopir[0:3] = np.array([0.0, 0.26469, 0.82925]) * Rv / 3.1

            ysplopir[3:7] = np.array(
                [
                    np.poly1d([2.13572e-04, 1.00270, -4.22809e-01])(Rv),
                    np.poly1d([-7.35778e-05, 1.00216, -5.13540e-02])(Rv),
                    np.poly1d([-3.32598e-05, 1.00184, 7.00127e-01])(Rv),
                    np.poly1d(
                        [1.19456, 1.01707, -5.46959e-03, 7.97809e-04, -4.45636e-05][
                            ::-1
                        ]
                    )(Rv),
                ]
            )

            tck = interpolate.splrep(
                np.hstack([xsplopir, xspluv]), np.hstack([ysplopir, yspluv]), k=3
            )
            k[ind] = interpolate.splev(x[ind], tck)

        # convert from A(lambda)/E(B-V) to A(lambda)/A(V)
        k /= Rv

        # setup the output
        if Alambda:
            return k * Av
        else:
            return k * Av * (np.log(10.0) * 0.4)


class Gordon03_SMCBar(ExtinctionLaw):
    """
    Gordon03 SMCBar extinction curve

    Average of 4 SMCBar extinction curves from
    Gordon et al. 2003 (ApJ, 594, 279).

    Note that Rv has no impact on this law: according to Gordon et al (2003),
    the average value of Rv is fixed to 2.74 +/- 0.13
    """

    def __init__(self):
        super().__init__()
        self.name = "Gordon03_SMCBar"
        self.Rv = 2.74
        self.x_range = [0.3, 10.0]

    def function(
        self, lamb, Av=1, Rv=2.74, Alambda=True, **kwargs
    ):
        """
        Gordon03_SMCBar extinction law

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which to evaluate the law.

        Av: float
            desired A(V) (default 1.0)

        Rv: float
            desired R(V) (default 2.74)

        Alambda: bool
            if set returns +2.5*1./log(10.)*tau, tau otherwise

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation as a function of wavelength
            depending on Alambda option +2.5*1./log(10.)*tau,  or tau
        """
        # ensure the units are in angstrom
        _lamb = units.Quantity(lamb, units.angstrom).value

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([lamb])
        else:
            _lamb = lamb[:]

        # convert to wavenumbers
        x = 1.0e4 / _lamb

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, self.x_range, self.name)

        # set Rv explicitly to the fixed value
        Rv = self.Rv

        c1 = -4.959 / Rv
        c2 = 2.264 / Rv
        c3 = 0.389 / Rv
        c4 = 0.461 / Rv
        x0 = 4.6
        gamma = 1.0

        k = np.zeros(np.size(x))

        # UV part
        xcutuv = 10000.0 / 2700.0
        xspluv = 10000.0 / np.array([2700.0, 2600.0])

        ind = np.where(x >= xcutuv)
        if np.size(ind) > 0:
            k[ind] = (
                1.0
                + c1
                + (c2 * x[ind])
                + c3
                * ((x[ind]) ** 2)
                / (((x[ind]) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((x[ind]) ** 2))
            )
            yspluv = (
                1.0
                + c1
                + (c2 * xspluv)
                + c3
                * ((xspluv) ** 2)
                / (((xspluv) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((xspluv) ** 2))
            )

        # FUV portion
        ind = np.where(x >= 5.9)
        if np.size(ind) > 0:
            k[ind] += c4 * (
                0.5392 * ((x[ind] - 5.9) ** 2) + 0.05644 * ((x[ind] - 5.9) ** 3)
            )

        # Opt/NIR part
        ind = np.where(x < xcutuv)
        if np.size(ind) > 0:
            xsplopir = np.zeros(9)
            xsplopir[0] = 0.0
            xsplopir[1:10] = 1.0 / np.array(
                [2.198, 1.65, 1.25, 0.81, 0.65, 0.55, 0.44, 0.37]
            )

            # Values directly from Gordon et al. (2003)
            # ysplopir =  np.array([0.0,0.016,0.169,0.131,0.567,0.801,
            #                      1.00,1.374,1.672])
            # K & J values adjusted to provide a smooth,
            #      non-negative cubic spline interpolation
            ysplopir = np.array(
                [0.0, 0.11, 0.169, 0.25, 0.567, 0.801, 1.00, 1.374, 1.672]
            )

            tck = interpolate.splrep(
                np.hstack([xsplopir, xspluv]), np.hstack([ysplopir, yspluv]), k=3
            )
            k[ind] = interpolate.splev(x[ind], tck)

        if Alambda:
            return k * Av
        else:
            return k * Av * (np.log(10.0) * 0.4)


class Gordon16_RvFALaw(ExtinctionLaw):
    """
    Gordon16 RvFA extinction law

    Mixture of Milky Way R(V) dependent extinction law and
    Gordon et al. (2003) SMCBar average extinction curve.

    This extinction curve model encompasses the average behavior of
    measured extinction curves in the Milky Way, LMC, and SMC.

    Implemented as a mixture between Fitzpatrick99 and Gordon03_SMCBar
    classes
    """

    def __init__(self):
        super().__init__()
        self.ALaw = Fitzpatrick99()
        self.BLaw = Gordon03_SMCBar()
        self.name = "Gordon16_RvFALaw"

        self.x_range = [
            np.max([self.ALaw.x_range[0], self.BLaw.x_range[0]]),
            np.min([self.ALaw.x_range[1], self.BLaw.x_range[1]]),
        ]

    def function(self, lamb, Av=1, Rv=3.1, Alambda=True, f_A=0.5, **kwargs):
        """
        Gordon16_RvFALaw

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which to evaluate the law.

        Av: float
            desired A(V) (default 1.0)

        Alambda: bool
            if set returns +2.5*1./log(10.)*tau, tau otherwise

        f_A: float
            set the mixture ratio between the two laws (default 0.5)

        Rv: float
            R(V) of mixture law (default to 3.1)

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation as a function of wavelength
            depending on Alambda option +2.5*1./log(10.)*tau,  or tau
        """
        # ensure the units are in angstrom
        _lamb = units.Quantity(lamb, units.angstrom).value

        # get the R(V) value for the A component
        Rv_A = self.get_Rv_A(Rv, f_A)

        # compute the two components
        k_A = self.ALaw.function(_lamb, Av=Av, Rv=Rv_A, Alambda=Alambda)
        k_B = self.BLaw.function(_lamb, Av=Av, Alambda=Alambda)

        # return the mixture
        return f_A * k_A + (1.0 - f_A) * k_B

    def get_Rv_A(self, Rv, f_A=0.5):
        """
        Calculate the R(V) of the A component given the R(V) of the mixture

        Rv_A is such that:

                1 / Rv = f_A / Rv_A + (1 - f_A) / Rv_B

                Rv_A = 1. / (1. / (Rv * f_A) - (1. - f_A) / (f_A * Rv_B))

                where Rv_B = 2.74 by definition (see Gordon03_SMCBar)

        Parameters
        ----------
        Rv: float
            R(V) of the mixture law

        f_A: float
            Mixture ratio between the two components [default is 0.5]

        Returns
        -------
        Rv_A : float
            R(V) of the A componet
        """

        Rv_B = self.BLaw.Rv
        if f_A > 0:
            return 1.0 / (1.0 / (Rv * f_A) - (1.0 - f_A) / (f_A * Rv_B))
        else:
            return 3.1  # not used, but needs to be non-zero

    def get_Rv(self, Rv_A, f_A=0.5):
        """
        Calculate the Rv of the mixture law given R(V) of the A component

        Rv_A is such that:

                1 / Rv = f_A / Rv_A + (1 - f_A) / Rv_B

                where Rv_B = 2.74 by definition (see Gordon03_SMCBar)

        Parameters
        ----------
        Rv_A: float
            R(V) of the A component

        f_A: float
            Mixture ratio between the two components [default is 0.5]

        Returns
        -------
        Rv : float
            R(V) of the mixture
        """

        Rv_B = self.BLaw.Rv
        return 1.0 / (f_A / Rv_A + (1 - f_A) / Rv_B)


class Generalized_RvFALaw(Gordon16_RvFALaw):
    """
    Generalized RvFA extinction law

    Mixture of R(V) dependent extinction law (`ALaw`) and
    average extinction curve (`BLaw`).

    This extinction curve model encompasses the average behavior of
    measured extinction curves in the Milky Way, LMC, and SMC.

    Implemented as ALaw=Fitzpatrick99() and BLaw=Gordon03_SMCBar()
    by default.
    """

    def __init__(self, ALaw=Fitzpatrick99(), BLaw=Gordon03_SMCBar()):
        super().__init__()
        self.ALaw = ALaw
        self.BLaw = BLaw
        self.name = "Generalized_RvFALaw:" + ALaw.name + "+" + BLaw.name

        self.x_range = [
            np.max([self.ALaw.x_range[0], self.BLaw.x_range[0]]),
            np.min([self.ALaw.x_range[1], self.BLaw.x_range[1]]),
        ]


class Generalized_DustExt(ExtinctionLaw):
    """
    Generalized extinction curve class to import classes from
    dust_extinction package.

    Accepts class name as string input (`curve`) for all `average` and
    Rv-dependent `parameter_averages` extinction curve classes.
    """

    def __init__(self, curve="F04"):
        super().__init__()
        self.name = "dustextpkg_" + curve
        if curve in dustext_par.__all__:
            self.extcurve_class = getattr(dustext_par, curve)
            self.hasRvParam = True
        elif curve in dustext_avg.__all__:
            self.extcurve_class = getattr(dustext_avg, curve)
            self.hasRvParam = False
            self.Rv = self.extcurve_class.Rv
        elif curve in dustext_extend.__all__:
            self.extcurve_class = getattr(dustext_extend, curve)
            if hasattr(self.extcurve_class, "Rv_range"):
                self.hasRvParam = True
            else:
                self.hasRvParam = False
                self.Rv = self.extcurve_class.Rv
        else:
            raise ValueError(
                curve
                + " class not found. \n"
                + "Valid dust_extinction package classes: "
                + " ".join(dustext_par.__all__ + dustext_avg.__all__)
            )

        self.x_range = self.extcurve_class.x_range

    def function(self, lamb, Av=1, Rv=3.1, Alambda=True, **kwargs):
        """
        Generalized Extinction Law

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which to evaluate the law.

        Av: float
            desired A(V) (default 1.0)

        Rv: float
            desired R(V) (default 3.1)
            ignored if self.hasRvParam=False; defaults to self.Rv

        Alambda: bool
            if set returns +2.5*1./log(10.)*tau, tau otherwise

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation as a function of wavelength
            depending on Alambda option +2.5*1./log(10.)*tau,  or tau
        """
        # ensure the units are in angstrom
        _lamb = units.Quantity(lamb, units.angstrom).value

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([lamb])
        else:
            _lamb = lamb[:]

        if self.hasRvParam:
            extcurve_obj = self.extcurve_class(Rv=Rv)
        else:
            extcurve_obj = self.extcurve_class()

        if Alambda:
            return extcurve_obj(_lamb * units.angstrom) * Av
        else:
            return extcurve_obj(_lamb * units.angstrom) * Av * (np.log(10.0) * 0.4)
