import pytest

from beast.physicsmodel.stars import stellib
from beast.physicsmodel.stars import isochrone
from beast.physicsmodel.creategrid import gen_spectral_grid_from_stellib_given_points


@pytest.mark.skip(reason="does not work yet")
def test_gen_spectral_grid_from_stellib_given_points():
    """
    Make sure it runs and returns a grid
    """
    osl = stellib.Kurucz()
    oiso = isochrone.padova2010()
    chunksize = 10000
    # as it is an interator, list does the actual loop
    list(gen_spectral_grid_from_stellib_given_points(osl, oiso.data,
                                                     chunksize=chunksize))


if __name__ == '__main__':
    test_gen_spectral_grid_from_stellib_given_points()
