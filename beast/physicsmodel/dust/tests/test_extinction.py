import numpy as np
import pytest

from beast.physicsmodel.dust import extinction


def test_extinction_generalRvFA_initialize():
    tlaw = extinction.Generalized_RvFALaw()
    if not isinstance(tlaw, extinction.Generalized_RvFALaw):
        raise AssertionError("Should use Generalized_RvFALaw extinction")

    lam = np.linspace(2.0e3, 1.0e4, 10)
    tlaw_vals = tlaw(lam, Av=1.0, Rv=4.0, f_A=0.8)
    orig = extinction.Gordon16_RvFALaw()
    orig_vals = orig(lam, Av=1.0, Rv=4.0, f_A=0.8)
    np.testing.assert_allclose(tlaw_vals, orig_vals)


@pytest.mark.parametrize("curve", ["F04", "G03_LMCAvg"])
def test_extinction_dustextpkg_initialize(curve):
    tlaw = extinction.Generalized_DustExt(curve)
    if not isinstance(tlaw, extinction.Generalized_DustExt):
        raise AssertionError("Should use Generalized_DustExt extinction")
