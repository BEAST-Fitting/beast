"""
Extinction Curves
"""
import numpy as np

from astropy import units

import dust_extinction.parameter_averages as dustext_par
import dust_extinction.averages as dustext_avg

from beast.config import __ROOT__
import beast.physicsmodel.dust.extinction_extension as dustext_extend

__all__ = [
    "ExtinctionLaw",
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


class Gordon16_RvFALaw(ExtinctionLaw):
    """
    Gordon16 RvFA extinction law from dust_extinction package.

    This extinction curve model encompasses the average behavior of
    measured extinction curves in the Milky Way, LMC, and SMC.

    Wrapper around dust_extinction.parameter_averages.G16.
    """

    def __init__(self):
        super().__init__()
        self.name = "Gordon16_RvFALaw"
        self.x_range = [0.3, 10.0]
        self.Rv = 3.1

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
        _lamb = units.Quantity(lamb, units.angstrom).value

        RvA = self.get_Rv_A(Rv, f_A)

        g16_model = dustext_par.G16(RvA=RvA, fA=f_A)

        if Alambda:
            result = g16_model(_lamb * units.angstrom) * Av
        else:
            result = g16_model(_lamb * units.angstrom) * Av * (np.log(10.0) * 0.4)

        return result

    def get_Rv_A(self, Rv, f_A=0.5):
        """
        Calculate the R(V) of the A component given the R(V) of the mixture

        Rv_A is such that:

                1 / Rv = f_A / Rv_A + (1 - f_A) / Rv_B

                Rv_A = 1. / (1. / (Rv * f_A) - (1. - f_A) / (f_A * Rv_B))

                where Rv_B = 2.74 by definition (see G03_SMCBar)

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

        Rv_B = 2.74
        if f_A > 0:
            return 1.0 / (1.0 / (Rv * f_A) - (1.0 - f_A) / (f_A * Rv_B))
        else:
            return 3.1

    def get_Rv(self, Rv_A, f_A=0.5):
        """
        Calculate the Rv of the mixture law given R(V) of the A component

        Rv_A is such that:

                1 / Rv = f_A / Rv_A + (1 - f_A) / Rv_B

                where Rv_B = 2.74 by definition (see G03_SMCBar)

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

        Rv_B = 2.74
        return 1.0 / (f_A / Rv_A + (1 - f_A) / Rv_B)


class Generalized_RvFALaw(Gordon16_RvFALaw):
    """
    Generalized RvFA extinction law

    Mixture of R(V) dependent extinction law (`ALaw`) and
    average extinction curve (`BLaw`).

    This extinction curve model encompasses the average behavior of
    measured extinction curves in the Milky Way, LMC, and SMC.

    Implemented as ALaw=Generalized_DustExt('F19') and BLaw=Generalized_DustExt('G03_SMCBar')
    by default.
    """

    def __init__(self, ALaw=None, BLaw=None):
        super().__init__()
        if ALaw is None:
            ALaw = Generalized_DustExt("G23")
        if BLaw is None:
            BLaw = Generalized_DustExt("G03_SMCBar")
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
        _lamb = units.Quantity(lamb, units.angstrom).value

        if isinstance(_lamb, float) or isinstance(_lamb, np.float64):
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
