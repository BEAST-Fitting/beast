import numpy as np
import astropy.units as u

from dust_extinction.parameter_averages import F19
from dust_extinction.averages import G03_SMCBar
from dust_extinction.grain_models import WD01, D03

from beast.physicsmodel.dust.extinction_extension import (
    F19_D03_extension,
    G03_SMCBar_WD01_extension,
)


def test_F19_D03_ext():
    emod = F19(Rv=3.1)
    dmod = D03(modelname="MWRV31")
    cmod = F19_D03_extension(Rv=3.1)

    x = np.arange(cmod.x_range[0], cmod.x_range[1], 0.1) / u.micron

    cmod_vals = cmod(x)
    dmod_vals = dmod(x)

    gvals_f19 = (x > emod.x_range[0] / u.micron) & (x < emod.x_range[1] / u.micron)
    emod_vals = emod(x[gvals_f19])

    # test that the combined model as the grain model values below 912 A
    # below the merge wavelengths
    gvals_ion = x > 1.0 / 0.0912 / u.micron
    np.testing.assert_allclose(cmod_vals[gvals_ion], dmod_vals[gvals_ion])

    # test that the combine dmodel has the F19 values above 1700 A
    # above the merge wavelengths
    gvals_amerge = (x < 1.0 / 0.1700 / u.micron) & (x > emod.x_range[0] / u.micron)
    gvals_amerge_f19 = x[gvals_f19] < 1.0 / 0.1700 / u.micron
    np.testing.assert_allclose(cmod_vals[gvals_amerge], emod_vals[gvals_amerge_f19])


def test_G03_SMCBar_WD01_ext():
    emod = G03_SMCBar()
    dmod = WD01(modelname="SMCBar")
    cmod = G03_SMCBar_WD01_extension()

    x = np.arange(cmod.x_range[0], cmod.x_range[1], 0.1) / u.micron

    cmod_vals = cmod(x)
    dmod_vals = dmod(x)

    gvals_g03 = (x > emod.x_range[0] / u.micron) & (x < emod.x_range[1] / u.micron)
    emod_vals = emod(x[gvals_g03])

    # test that the combined model as the grain model values below 912 A
    # below the merge wavelengths
    gvals_ion = x > 1.0 / 0.0912 / u.micron
    np.testing.assert_allclose(cmod_vals[gvals_ion], dmod_vals[gvals_ion])

    # test that the combine dmodel has the F19 values above 1700 A
    # above the merge wavelengths
    gvals_amerge = (x < 1.0 / 0.1700 / u.micron) & (x > emod.x_range[0] / u.micron)
    gvals_amerge_g03 = x[gvals_g03] < 1.0 / 0.1700 / u.micron
    np.testing.assert_allclose(cmod_vals[gvals_amerge], emod_vals[gvals_amerge_g03])
