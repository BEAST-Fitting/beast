import copy
import warnings
import numpy as np
import astropy.units as u

from dust_extinction.helpers import _test_valid_x_range
from dust_extinction.parameter_averages import F19
from dust_extinction.averages import G03_SMCBar
from dust_extinction.grain_models import D03, WD01

__all__ = ["F19_D03_extension", "G03_SMCBar_WD01_extension"]


class SpectralUnitsWarning(UserWarning):
    pass


def _get_x_in_wavenumbers(in_x):
    """
    Convert input x to wavenumber given x has units.
    Otherwise, assume x is in waveneumbers and issue a warning to this effect.
    Parameters
    ----------
    in_x : astropy.quantity or simple floats
        x values
    Returns
    -------
    x : floats
        input x values in wavenumbers w/o units
    """
    # handles the case where x is a scaler
    in_x = np.atleast_1d(in_x)

    # check if in_x is an astropy quantity, if not issue a warning
    if not isinstance(in_x, u.Quantity):
        warnings.warn(
            "x has no units, assuming x units are inverse microns",
            SpectralUnitsWarning
        )

    # convert to wavenumbers (1/micron) if x input in units
    # otherwise, assume x in appropriate wavenumber units
    with u.add_enabled_equivalencies(u.spectral()):
        x_quant = u.Quantity(in_x, 1.0 / u.micron, dtype=np.float64)

    # strip the quantity to avoid needing to add units to all the
    #    polynomical coefficients
    return x_quant.value


class F19_D03_extension(F19):
    r"""
    dust_extinction.parameter_averages.F19 model extended to shorter
    wavelengths using the dust_extinction.grain_models.D03 models.

    Parameters
    ----------
    None

    Raises
    ------
    InputParameterError
       Input Rv values outside of defined range

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        import astropy.units as u

        from beast.physicsmodel.dust.extinction_extension import F19_D03_extension
        from dust_extinction.grain_models import D03

        fig, ax = plt.subplots()

        # temp model to get the correct x range
        text_model = F19_D03_extension()

        # generate the curves and plot them
        x = np.arange(text_model.x_range[0], text_model.x_range[1], 0.1) / u.micron

        Rvs = [2.0, 3.0, 4.0, 5.0, 6.0]
        for cur_Rv in Rvs:
            ext_model = F19_D03_extension(Rv=cur_Rv)
            ax.plot(1.0 / x, ext_model(x), label="F19_D03_ext R(V) = " + str(cur_Rv))

        pmods = ["MWRV31", "MWRV40", "MWRV55"]
        for cmod in pmods:
            dmod = D03(modelname=cmod)
            ax.plot(1.0 / x, dmod(x), label=f"D03 {cmod}", linestyle="dashed", color="black")

        ax.set_xlabel(r"$\lambda$ [$\mu m$]")
        ax.set_ylabel(r"$A(x)/A(V)$")

        ax.set_xscale("log")

        ax.legend(loc="best")
        plt.show()
    """

    # update the wavelength range (in micron^-1)
    x_range = [0.3, 1.0 / 0.01]

    def evaluate(self, in_x, Rv):
        """
        F19_D03_extension function

        Parameters
        ----------
        in_x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in wavenumbers [1/micron]

           internally wavenumbers are used

        Returns
        -------
        axav: np array (float)
            A(x)/A(V) extinction curve [mag]

        Raises
        ------
        ValueError
           Input x values outside of defined range
        """
        # convert to wavenumbers (1/micron) if x input in units
        # otherwise, assume x in appropriate wavenumber units
        x = _get_x_in_wavenumbers(in_x)

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, self.x_range, self.__class__.__name__)

        # just in case someone calls evaluate explicitly
        Rv = np.atleast_1d(Rv)

        # ensure Rv is a single element, not numpy array
        Rv = Rv[0]

        # determine the dust grain models to use for the input Rv
        if Rv < 4.0:
            d1rv = 3.1
            d2rv = 4.0
            d1mod = D03(modelname="MWRV31")
            d2mod = D03(modelname="MWRV40")
        else:
            d1rv = 4.0
            d2rv = 5.5
            d1mod = D03(modelname="MWRV40")
            d2mod = D03(modelname="MWRV55")

        # interpolate to get the model extinction for the input Rv value
        dslope = (d2mod(in_x) - d1mod(in_x)) / (d2rv - d1rv)
        dmod = d1mod(in_x) + dslope * (Rv - d1rv)

        # compute the F19 curve for the input Rv over the F19 defined wavelength range
        gvals_f19 = (x > super().x_range[0]) & (x < super().x_range[1])
        fmod = super().evaluate(in_x[gvals_f19], Rv)

        # now merge the two smoothly
        outmod = copy.copy(dmod)
        outmod[gvals_f19] = fmod

        merge_range = np.array([1.0 / 0.1675, super().x_range[1]])
        gvals_merge = (x > merge_range[0]) & (x < merge_range[1])
        # have weights be zero at the min merge and 1 at the max merge
        weights = (x[gvals_merge] - merge_range[0]) / (merge_range[1] - merge_range[0])
        outmod[gvals_merge] = (1.0 - weights) * outmod[gvals_merge] + weights * dmod[
            gvals_merge
        ]

        return outmod


class G03_SMCBar_WD01_extension(G03_SMCBar):
    r"""
    dust_extinction.averages.G03_SMCBar model extended to shorter
    wavelengths using the dust_extinction.grain_models.WD01 SMCBar model.

    Parameters
    ----------
    Rv: float
        R(V) = A(V)/E(B-V) = total-to-selective extinction

    Raises
    ------
    InputParameterError
       Input Rv values outside of defined range

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        import astropy.units as u

        from beast.physicsmodel.dust.extinction_extension import G03_SMCBar_WD01_extension
        from dust_extinction.grain_models import WD01

        fig, ax = plt.subplots()

        # define the extinction model
        ext_model = G03_SMCBar_WD01_extension()

        # generate the curves and plot them
        x = np.arange(ext_model.x_range[0], ext_model.x_range[1], 0.1) / u.micron

        ax.plot(1.0 / x, ext_model(x), label="G03 SMCBar WD01 ext")

        dmod = WD01(modelname="SMCBar")
        ax.plot(
            1.0 / x, dmod(x), label="WD01 SMCBar", linestyle="dashed", color="black"
        )

        ax.set_xlabel(r"$\lambda$ [$\mu m$]")
        ax.set_ylabel(r"$A(x)/A(V)$")

        ax.set_xscale("log")

        ax.legend(loc="best")
        plt.show()
    """

    # update the wavelength range (in micron^-1)
    x_range = [0.3, 1.0 / 0.01]

    def evaluate(self, in_x):
        """
        G03_SMCBar_WD01_extension function

        Parameters
        ----------
        in_x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in wavenumbers [1/micron]

           internally wavenumbers are used

        Returns
        -------
        axav: np array (float)
            A(x)/A(V) extinction curve [mag]

        Raises
        ------
        ValueError
           Input x values outside of defined range
        """
        # convert to wavenumbers (1/micron) if x input in units
        # otherwise, assume x in appropriate wavenumber units
        x = _get_x_in_wavenumbers(in_x)

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, self.x_range, self.__class__.__name__)

        # compute the dust grain model
        dmodel = WD01(modelname="SMCBar")
        dmod = dmodel(in_x)

        # compute the F19 curve for the input Rv over the F19 defined wavelength range
        gvals_g03 = (x > super().x_range[0]) & (x < super().x_range[1])
        fmod = super().evaluate(in_x[gvals_g03])

        # now merge the two smoothly
        outmod = copy.copy(dmod)
        outmod[gvals_g03] = fmod

        merge_range = np.array([1.0 / 0.1675, super().x_range[1]])
        gvals_merge = (x > merge_range[0]) & (x < merge_range[1])
        # have weights be zero at the min merge and 1 at the max merge
        weights = (x[gvals_merge] - merge_range[0]) / (merge_range[1] - merge_range[0])
        outmod[gvals_merge] = (1.0 - weights) * outmod[gvals_merge] + weights * dmod[
            gvals_merge
        ]

        return outmod
