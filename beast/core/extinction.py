"""
EXTINCTION MANAGEMENT
Here are defined extinction functions that relies
on given laws.
"""
import numpy as np
from scipy import interpolate, interp

from ..observationmodel import phot
from ..tools.helpers import val_in_unit
from ..external.ezunits import unit
from ..config import __ROOT__

__version__ = '1.0'
__all__ = ['Calzetti', 'Cardelli', 'ExtinctionLaw', 'Fitzpatrick99', 'Gordon03_SMCBar']

libdir = __ROOT__ + '/libs/'

class ExtinctionLaw(object):
    """ Template class """
    def __init__(self):
        self.name = 'None'

    def function(self, lamb, *arg, **kwargs):
        """ expected to contain a function of lambda that return the
        extinction values
        """
        raise NotImplementedError

    def inFilter(self, names, filterLib=None, *args, **kwargs):
        """
        returns the extinction value for a given filter band or filter color
        colors (e.g. U, U-B ...)

        Parameters
        ----------
        names: str or list(str) or list(filters)
            filter names or filter instances to evaluate. a name can be a color such as 'U-B'

        filterLib: filepath
            path to the filter library hd5 file (default is the internal library)

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation value or array of values
        """
        if (type(names) == str):
            if '-' in names:
                _filters = phot.load_filters(names.replace(' ', '').split('-'), interp=True, filterLib=filterLib)
                return float(self.function(_filters[0].cl, *args, **kwargs) - self.function(_filters[1].cl,*args, **kwargs))
            else:
                _filters = phot.load_filters([names], interp=True, filterLib=filterLib)
                return float(self.function(_filters[0].cl, *args, **kwargs))

        elif hasattr(names, '__iter__'):
            if (type(names[0]) == str):
                # assumes all entires are str
                lst = np.unique(' '.join(names).replace('-', ' ').split())
                _filters = phot.load_filters(lst, interp=True, filterLib=filterLib)
                #in case of python 2.6xx, explicit loop
                d = {}
                for lk, fk in zip(lst, _filters):
                    d[lk] = self.function(fk.cl, *args, **kwargs)
                return np.asarray([ float(eval(lk, d)) for lk in names ])
            else:
                # assumes list of Filter instances
                return np.asarray([ float(self.function(fk.cl, *args, **kwargs)) for fk in names ])

    def __call__(self, *args, **kwargs):
        return self.function(*args, **kwargs)

    def isvalid(self, *args, **kwargs):
        """ Check if the current arguments are in the validity domain of the law
        Must be redefined if any restriction applies to the law
        """
        return True


class Calzetti(ExtinctionLaw):
    """
    Calzetti et al.  (2000, ApJ 533, 682) developed a recipe for dereddening the
    spectra of galaxies where massive stars dominate the radiation output, valid
    between 0.12 to 2.2 microns.
    Extrapolation down to 0.0912 microns

    Note that the supplied color excess should be that derived for the
    stellar  continuum, :math:`EBV(stars)`, which is related to the reddening
    derived from the gas, :math:`EBV(gas)`, via the Balmer decrement by
    :math:`EBV(stars) = 0.44 \\times EBV(gas)`

    :math:`R_V` - Ratio of total to selective extinction, default is 4.05.
    Calzetti et al. (2000) estimate :math:`R_V = 4.05 \pm 0.80` from optical-IR
    observations of 4 starbursts.
    """
    def __init__(self):
        self.name = 'Calzetti'

    def function(self, lamb, Av=1, Rv=4.05, Alambda=True, **kwargs):
        """
        Returns Alambda or tau for a Calzetti law Lamb is input in Angstroms

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which evaluate the law.

        Av: float
            desired A(V) (default 1.0)

        Rv: float
            desired R(V) (default 4.05)

        Alambda: bool
            if set returns +2.5 * 1. / log(10.) * tau, tau otherwise

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation as a function of wavelength
            depending on Alambda option +2.5*1./log(10.)*tau,  or tau
        """
        # handle units
        _lamb = val_in_unit('lamb', lamb, 'angstrom').magnitude

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([lamb])
        else:
            _lamb = lamb[:]

        x = 1.e4 / _lamb  # wavenumber in um^-1
        k = np.zeros(np.size(x))

        ind = np.where( (_lamb >= 0.630 ) & (_lamb <= 2.2) )
        k[ind] = 2.659 * (-1.857 + 1.040 * x[ind]) + Rv

        ind = np.where((_lamb >= 0.0912 ) & (_lamb < 0.630) )
        k[ind] = 2.659 * (-2.156 + 1.509 * x[ind] - 0.198 * x[ind] ** 2 + 0.011 * x[ind] ** 3 ) + Rv

        if Alambda:
            return k
        else:
            return 10 ** (0.4 * k)


class Cardelli(ExtinctionLaw):
    """ Cardelli, Clayton, and Mathis (1989, ApJ, 345, 245)"""
    def __init__(self):
        self.name = 'Cardelli'

    def function(self, lamb, Av=1., Rv=3.1, Alambda=True, **kwargs):
        """
        Cardelli extinction Law

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which evaluate the law.

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
        _lamb = val_in_unit('lamb', lamb, 'angstrom').magnitude

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([lamb])
        else:
            _lamb = lamb[:]

        # init variables
        x = 1.e4 / _lamb  # wavenumber in um^-1
        a = np.zeros(np.size(x))
        b = np.zeros(np.size(x))
        # Infrared (Eq 2a,2b)
        ind = np.where((x >= 0.3) & (x < 1.1))
        a[ind] =  0.574 * x[ind] ** 1.61
        b[ind] = -0.527 * x[ind] ** 1.61
        # Optical & Near IR
        # Eq 3a, 3b
        ind = np.where((x >= 1.1) & (x <= 3.3))
        y = x[ind] - 1.82
        a[ind] = 1. + 0.17699 * y - 0.50447 * y ** 2 - 0.02427 * y ** 3 + 0.72085 * y ** 4 + 0.01979 * y ** 5 - 0.77530 * y ** 6 + 0.32999 * y ** 7
        b[ind] =      1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3 - 5.38434 * y ** 4 - 0.62251 * y ** 5 + 5.30260 * y ** 6 - 2.09002 * y ** 7
        # UV
        # Eq 4a, 4b
        ind = np.where((x >= 3.3) & (x <= 8.0))
        a[ind] =  1.752 - 0.316 * x[ind] - 0.104 / ((x[ind] - 4.67) ** 2 + 0.341)
        b[ind] = -3.090 + 1.825 * x[ind] + 1.206 / ((x[ind] - 4.62) ** 2 + 0.263)

        ind = np.where((x >= 5.9) & (x <= 8.0))
        Fa     = -0.04473 * (x[ind] - 5.9) ** 2 - 0.009779 * (x[ind] - 5.9) ** 3
        Fb     =  0.21300 * (x[ind] - 5.9) ** 2 + 0.120700 * (x[ind] - 5.9) ** 3
        a[ind] = a[ind] + Fa
        b[ind] = b[ind] + Fb
        # Far UV
        # Eq 5a, 5b
        ind = np.where((x >= 8.0) & (x <= 10.0))
        # Fa = Fb = 0
        a[ind] = -1.073 - 0.628 * (x[ind] - 8.) + 0.137 * ((x[ind] - 8.) ** 2) - 0.070 * (x[ind] - 8.) ** 3
        b[ind] = 13.670 + 4.257 * (x[ind] - 8.) + 0.420 * ((x[ind] - 8.) ** 2) + 0.374 * (x[ind] - 8.) ** 3

        # Case of -values x out of range [0.3,10.0]
        ind = np.where((x > 10.0) | (x < 0.3))
        a[ind] = 0.0
        b[ind] = 0.0

        # Return Extinction vector
        # Eq 1
        if (Alambda):
            return( ( a + b / Rv ) * Av)
        else:
            # return( 1./(2.5 * 1. / np.log(10.)) * ( a + b / Rv ) * Av)
            return( 0.4 * np.log(10.) * ( a + b / Rv ) * Av)


class Fitzpatrick99(ExtinctionLaw):
    """
    Fitzpatrick (1999, PASP, 111, 63)
    R(V) dependent extinction curve that explicitly deals with optical/NIR
    extinction being measured from broad/medium band photometry.
    Based on fm_unred.pro from the IDL astronomy library
    """
    def __init__(self):
        self.name = 'Fitzpatrick99'

    def function(self, lamb, Av=1, Rv=3.1, Alambda=True, draine_extend=False, **kwargs):
        """
        Fitzpatrick99 extinction Law
        Lamb is input in Anstroms

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which evaluate the law.

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
        _lamb = val_in_unit('lamb', lamb, 'angstrom').magnitude

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([lamb])
        else:
            _lamb = lamb[:]

        c2 = -0.824 + 4.717 / Rv
        c1 = 2.030 - 3.007 * c2
        c3 = 3.23
        c4 = 0.41
        x0 = 4.596
        gamma = 0.99

        x = 1.e4 / _lamb
        k = np.zeros(np.size(x))

        # compute the UV portion of A(lambda)/E(B-V)
        xcutuv = 10000.0 / 2700.
        xspluv = 10000.0 / np.array([2700., 2600.])
        ind = np.where(x >= xcutuv)

        if np.size(ind) > 0:
            k[ind] = c1 + (c2 * x[ind]) + c3 * ((x[ind]) ** 2) / ( ((x[ind]) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((x[ind]) ** 2 ))
            yspluv = c1 + (c2 * xspluv) + c3 * ((xspluv) ** 2) / ( ((xspluv) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((xspluv) ** 2 ))

            # FUV portion
            if not draine_extend:
                fuvind = np.where(x >= 5.9)
                k[fuvind] += c4 * (0.5392 * ((x[fuvind] - 5.9) ** 2) + 0.05644 * ((x[fuvind] - 5.9) ** 3))

            k[ind] += Rv
            yspluv += Rv

        # Optical/NIR portion

        ind = np.where(x < xcutuv)
        if np.size(ind) > 0:
            xsplopir = np.zeros(7)
            xsplopir[0] = 0.0
            xsplopir[1: 7] = 10000.0 / np.array([26500.0, 12200.0, 6000.0, 5470.0, 4670.0, 4110.0])

            ysplopir = np.zeros(7)
            ysplopir[0: 3] = np.array([0.0, 0.26469, 0.82925]) * Rv / 3.1

            ysplopir[3: 7] = np.array([np.poly1d([2.13572e-04, 1.00270, -4.22809e-01])(Rv),
                                       np.poly1d([-7.35778e-05, 1.00216, -5.13540e-02])(Rv),
                                       np.poly1d([-3.32598e-05, 1.00184,  7.00127e-01])(Rv),
                                       np.poly1d([ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04,
                                       -4.45636e-05][::-1])(Rv)])

            tck = interpolate.splrep(np.hstack([xsplopir, xspluv]), np.hstack([ysplopir, yspluv]), k=3)
            k[ind] = interpolate.splev(x[ind], tck)

        # convert from A(lambda)/E(B-V) to A(lambda)/A(V)
        k /= Rv

        # FUV portion from Draine curves
        if draine_extend:
            fuvind = np.where(x >= 5.9)
            tmprvs = np.arange(2.,6.1,0.1)
            diffRv = Rv - tmprvs
            if min(abs(diffRv)) < 1e-8:
                l_draine, k_draine = np.loadtxt(libdir+'MW_Rv%s_ext.txt' % ("{0:.1f}".format(Rv)), usecols=(0,1), unpack=True)
            else: 
                add, = np.where(diffRv < 0.)
                Rv1 = tmprvs[add[0]-1]
                Rv2 = tmprvs[add[0]]
                l_draine, k_draine1 = np.loadtxt(libdir+'MW_Rv%s_ext.txt' % ("{0:.1f}".format(Rv1)), usecols=(0,1), unpack=True)
                l_draine, k_draine2 = np.loadtxt(libdir+'MW_Rv%s_ext.txt' % ("{0:.1f}".format(Rv2)), usecols=(0,1), unpack=True)
                frac = diffRv[add[0]-1]/(Rv2-Rv1) 
                k_draine = (1. - frac)*k_draine1 + frac*k_draine2
            
            dind = np.where((1./l_draine) >= 5.9)
            k[fuvind] = interp(x[fuvind],1./l_draine[dind][::-1],k_draine[dind][::-1])

        # setup the output
        if (Alambda):
            return(k * Av)
        else:
            return(k * Av * (np.log(10.) * 0.4))


class Gordon03_SMCBar(ExtinctionLaw):
    """ Gordon et al. 2003 (ApJ, 594:279-293)
    Note that Rv has no impact on this law: according to Gordon et al (2003),
    the average value of Rv is fixed to 2.74 +/- 0.13
    """
    def __init__(self):
        self.name = 'Gordon03_SMCBar'
        self.Rv = 2.74

    def function(self, lamb, Av=1, Rv=2.74, Alambda=True, draine_extend=False,  **kwargs):
        """
        Lamb is input in Anstroms
        Note that Rv is not given as a variable in the paper of reference

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which evaluate the law.

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
        _lamb = val_in_unit('lamb', lamb, 'angstrom').magnitude

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([lamb])
        else:
            _lamb = lamb[:]

        c1 = -4.959 / Rv
        c2 = 2.264 / Rv
        c3 = 0.389 / Rv
        c4 = 0.461 / Rv
        x0 = 4.6
        gamma = 1.0

        x = 1.e4 / _lamb
        k = np.zeros(np.size(x))

        # UV part
        xcutuv = 10000.0 / 2700.
        xspluv = 10000.0 / np.array([2700., 2600.])

        ind = np.where(x >= xcutuv)
        if np.size(ind) > 0:
            k[ind] = 1.0 + c1 + (c2 * x[ind]) + c3 * ((x[ind]) ** 2) / ( ((x[ind]) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((x[ind]) ** 2 ))
            yspluv = 1.0 + c1 + (c2 * xspluv) + c3 * ((xspluv) ** 2) / ( ((xspluv) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((xspluv) ** 2 ))
            
            # FUV portion  
            if draine_extend:
                ind = np.where(x >= 5.9)
                l_draine, k_draine = np.loadtxt(libdir+'SMC_Rv2.74_norm.txt',usecols=(0,1),unpack=True)
                dind = np.where((1./l_draine) >= 5.9)
                k[ind] = interp(x[ind],1./l_draine[dind][::-1],k_draine[dind][::-1])
            else:
                k[ind] += c4 * (0.5392 * ((x[ind] - 5.9) ** 2) + 0.05644 * ((x[ind] - 5.9) ** 3))

        # Opt/NIR part
        ind = np.where(x < xcutuv)
        if np.size(ind) > 0:
            xsplopir = np.zeros(9)
            xsplopir[0] = 0.0
            xsplopir[1: 10] = 1.0 / np.array([2.198, 1.65, 1.25, 0.81, 0.65, 0.55, 0.44, 0.37])

            # Values directly from Gordon et al. (2003)
            # ysplopir =  np.array([0.0,0.016,0.169,0.131,0.567,0.801,1.00,1.374,1.672])
            # K & J values adjusted to provide a smooth, non-negative cubic spline interpolation
            ysplopir = np.array([0.0, 0.11, 0.169, 0.25, 0.567, 0.801, 1.00, 1.374, 1.672])

            tck = interpolate.splrep(np.hstack([xsplopir, xspluv]), np.hstack([ysplopir, yspluv]), k=3)
            k[ind] = interpolate.splev(x[ind], tck)

        if (Alambda):
            return(k * Av)
        else:
            return(k * Av * (np.log(10.) * 0.4 ))


class RvFbumpLaw(ExtinctionLaw):
    """ Mixture of extinction laws allowing to vary the bump amplitude in the
    extinction law
    Default is a Mixture between Fitzpatrick99 and Gordon03_SMCBar

   returns f_bump * RvLaw(*args, **kwargs) + (1. - f_bump) * NoBumpLaw(*args, **kwargs)
    """
    def __init__(self, RvLaw=None, NoBumpLaw=None, name=None):
        """ Constructor

        Parameters
        ----------
        RvLaw: ExtinctionLaw
            Component which models attenuation related to the bump
            (default Fitzpatrick99)

        NoBumpLaw: ExtinctionLaw
            Bumpless component of the models
            (default Gordon03_SMCBar)
        """
        self.RvLaw = RvLaw or Fitzpatrick99()
        self.NoBumpLaw = NoBumpLaw or Gordon03_SMCBar()
        self.name = name or 'RvFbumpLaw'

    def function(self, lamb, Av=1, Rv_A=None, Alambda=True, f_A=0.5, Rv_B=None,
                 Rv=None, **kwargs):
        """
        Lamb as to be in Angstroms!!!

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which evaluate the law.

        Av: float
            desired A(V) (default 1.0)

        Alambda: bool
            if set returns +2.5*1./log(10.)*tau, tau otherwise

        f_A: float
            set the mixture ratio between the two laws (default 0.5)

        Rv_A: float
            extinction param. on the RvLaw (Bump)

        Rv_B: float
            extinction param. on the bumpless component

        Rv: float
            effective R(V) according to the mixture

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation as a function of wavelength
            depending on Alambda option +2.5*1./log(10.)*tau,  or tau

        .. math::
            f_A * RvLaw(*args, **kwargs) + (1. - f_A) * NoBumpLaw(*args, **kwargs)
        """
        _lamb = val_in_unit('lamb', lamb, 'angstrom').magnitude

        if Rv_A is None:
            Rv_A = getattr(self.RvLaw, 'Rv', None)

        if Rv_B is None:
            Rv_B = getattr(self.NoBumpLaw, 'Rv', None)
        #print(Rv, Rv_A, Rv_B, f_A)

        if sum([Rv_A is None, Rv_B is None, Rv is None]) >= 2:
            raise ValueError('Must provide at least 2 Rv values')

        if Rv_A is None:
            Rv_A = self.get_Rv_A(Rv, f_A, Rv_B)
        if Rv_B is None:
            Rv_B = self.get_Rv_B(Rv, Rv_A, f_A)

        #print(Rv, Rv_A, Rv_B, f_A)

        return f_A * self.RvLaw.function(_lamb, Av=Av, Rv=Rv_A, Alambda=Alambda) + (1. - f_A) * self.NoBumpLaw.function(_lamb, Av=Av, Alambda=Alambda, Rv=Rv_B)

    def isvalid(self, Av=None, Rv=None, f_A=0.5, Rv_A=None, Rv_B=None):
        """ Test the validity of an extinction vector (Av, Rv, Rv_A, Rv_B, fbump)

        .. math::
            Law = f_A * RvLaw (lamb, Av=Av, Rv=Rv_A) + (1. - f_A) * NoBumpLaw(lamb, Av=Av, Rv=Rv_B)

        The validity impose :math:`R_V` ranges and  to be a fraction (i.e.,
        between 0 and 1)

        At least 2 out of the 3 :math:`R_V` values must be provided, the 3rd
        will be computed if missing.

        Parameters
        ----------
        Av: float
            Av value (any value is allowed, even <0)

        Rv, Rv_A, Rv_B: float, float, float
            effective Rv, RvLaw component and bumpless component Rv values, respectively.
            At least 2 must be provided.

        f_A: float
            Mixture ratio between the two components

        Returns
        -------
        r: bool
            True, if the values a coherent with the definition.
        """

        if Rv_B is None and hasattr(self.NoBumpLaw, 'Rv'):
            Rv_B = self.NoBumpLaw.Rv

        if Rv_A is None and hasattr(self.RvLaw, 'Rv'):
            Rv_A = self.RvLaw.Rv

        # if we do not have at least 2 of the 3 Rvs defined then it's invalid
        if sum([Rv_A is None, Rv_B is None, Rv is None]) >= 2:
            return False

        if Rv_A is None:
            Rv_A = self.get_Rv_A(Rv, f_A, Rv_B=Rv_B)
        if Rv is None:
            Rv = self.get_Rv(Rv_A, f_A, Rv_B=Rv_B)
        if Rv_B is None:
            Rv_B = self.get_Rv_B(Rv, Rv_A, f_A)

        # f_A is a fraction and any Rv is limited to [2.0, 6.0]
        return (0. <= f_A <= 1.) & (2.0 <= Rv_B <= 6.0) & (2.0 <= Rv_A <= 6.0) & (2.0 <= Rv <= 6.0)

    def get_Rv_A(self, Rv, f_A=0.5, Rv_B=None):
        """ Returns the equivalent Rv to use in the bump component
            Law = f_A * RvLaw (lamb, Av=Av, Rv=Rv_A) + (1. - f_A) * NoBumpLaw(lamb, Av=Av, Rv=Rv_B)

            and Rv_A is such that:

            ..math::

                1 / Rv = f_A / Rv_A + (1 - f_A) / Rv_B

                Rv_A = 1. / (1. / (Rv * f_A) - (1. - f_A) / (f_A * Rv_B))

            not that Gordon03_SMCBar has a fixed Rv=2.74
        """

        if Rv_B is None and hasattr(self.NoBumpLaw, 'Rv'):
            Rv_B = self.NoBumpLaw.Rv

        return 1. / (1. / (Rv * f_A) - (1. - f_A) / (f_A * Rv_B))

    def get_Rv(self, Rv_A=None, f_A=0.5, Rv_B=None):
        """ Returns the equivalent effective Rv according to the mixture

        ..math::
            Law = f_A * RvLaw (lamb, Av=Av, Rv=Rv_A) + (1. - f_A) * NoBumpLaw(lamb, Av=Av, Rv=Rv_B)

        and Rv is such that:

        ..math::

            1 / Rv = f_A / Rv_A + (1 - f_A) / Rv_B

            Rv_A = 1. / (1. / (Rv * f_A) - (1. - f_A) / (f_A * Rv_B))

        note that Gordon03_SMCBar has a fixed Rv=2.74
        """
        if Rv_B is None and hasattr(self.NoBumpLaw, 'Rv'):
            Rv_B = self.NoBumpLaw.Rv

        if Rv_A is None and hasattr(self.RvLaw, 'Rv'):
            Rv_A = self.RvLaw.Rv

        return 1. / (f_A / Rv_A + (1 - f_A) / Rv_B)

    def get_Rv_B(self, Rv, Rv_A=None, f_A=0.5):
        """ Returns the equivalent Rv to use in the bumpless component

        .. math::
            Law = f_A * RvLaw (lamb, Av=Av, Rv=Rv_A) + (1. - f_A) * NoBumpLaw(lamb, Av=Av, Rv=Rv_B)

        and Rv_B is such that

        .. math::
            1 / Rv = f_A / Rv_A + (1 - f_A) / Rv_B

            Rv_A = 1. / (1. / (Rv * f_A) - (1. - f_A) / (f_A * Rv_B))

        note that Gordon03_SMCBar has a fixed Rv=2.74
        """
        if Rv_A is None and hasattr(self.RvLaw, 'Rv'):
            Rv_A = self.RvLaw.Rv

        return (1. - f_A) / (1. / Rv - f_A / Rv_A)


def testunit():
    # check that things look correct
    # -> make some plots
    import pylab

    x = np.arange(0.1, 10, 0.1)   # in um^-1
    lamb = 1.e4 / x * unit['angstrom']

    # ccm  = Cardelli()
    f99  = Fitzpatrick99()
    gsmc = Gordon03_SMCBar()

    fig = pylab.figure()
    plt = fig.add_subplot(1, 1, 1)

    Rv_vals = np.arange(2, 6, dtype=float)
    for Rv in Rv_vals:
        # yccm = ccm.function(lamb, Rv=Rv)
        yf99 = f99.function(lamb, Rv=Rv)

        # plt.plot(x,yccm,label='CCM, Rv=%0.1f' % (Rv) )
        plt.plot(x, yf99, label='F99, Rv=%0.1f' % (Rv) )

    ygsmc = gsmc.function(lamb)
    plt.plot(x, ygsmc, label='G. SMC')

    mixlaw = RvFbumpLaw()
    ymix = mixlaw(lamb, Rv=3.1, f_bump=0.75)
    plt.plot(x, ymix, label='Mixture f(bump)=0.75')

    ymix = mixlaw(lamb, Rv=3.1, f_bump=0.5)
    plt.plot(x, ymix, label='Mixture f(bump)=0.5')

    ymix = mixlaw(lamb, Rv=3.1, f_bump=0.25)
    plt.plot(x, ymix, label='Mixture f(bump=0.25')

    plt.set_ylabel('A($\lambda$)/A(V)')
    plt.set_xlabel('1/x [$\mu$m$^{-1}$]')

    plt.legend(loc=0, frameon=False)

    pylab.show()

if __name__ == "__main__":
    testunit()
