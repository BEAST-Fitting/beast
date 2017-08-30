"""
Attenuation Curves
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from scipy import interpolate, interp

from astropy import units

from ...observationmodel import phot

from .extinction import ExtinctionLaw

__all__ = ['Calzetti00']

class Calzetti00(ExtinctionLaw):
    """
    Calzetti et al.  (2000, ApJ 533, 682) developed a recipe for dereddening the
    spectra of galaxies where massive stars dominate the radiation output, valid
    between 0.12 to 2.2 microns.
    Extrapolation down to 0.0912 microns

    Note that the supplied color excess should be that derived for the
    stellar  continuum, `EBV(stars)`, which is related to the reddening
    derived from the gas, `EBV(gas)`, via the Balmer decrement by
    `EBV(stars) = 0.44 \\times EBV(gas)`

    `R_V` - Ratio of total to selective extinction, default is 4.05.
    Calzetti et al. (2000) estimate `R_V = 4.05 +/- 0.80` from optical-IR
    observations of 4 starbursts.
    """
    def __init__(self):
        self.name = 'Calzetti'

    def function(self, lamb, Av=1, Rv=4.05, Alambda=True, **kwargs):
        """
        Returns A(lambda) or tau for the Calzetti Law

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
        # ensure the units are in angstrom
        _lamb = units.Quantity(lamb, units.angstrom).value
        #_lamb = val_in_unit('lamb', lamb, 'angstrom').magnitude

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([lamb])
        else:
            _lamb = lamb[:]

        x = 1.e4 / _lamb  # wavenumber in um^-1
        k = np.zeros(np.size(x))

        ind = np.where( (_lamb >= 0.630 ) & (_lamb <= 2.2) )
        k[ind] = 2.659 * (-1.857 + 1.040 * x[ind]) + Rv

        ind = np.where((_lamb >= 0.0912 ) & (_lamb < 0.630) )
        k[ind] = 2.659 * (-2.156 + 1.509 * x[ind] - 0.198 * x[ind] ** 2 +
                          0.011 * x[ind] ** 3 ) + Rv

        if Alambda:
            return k
        else:
            return 10 ** (0.4 * k)


