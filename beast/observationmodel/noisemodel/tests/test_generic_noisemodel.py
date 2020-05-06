from astropy.tests.helper import remote_data

from beast.tests.helpers import download_rename
from beast.observationmodel.noisemodel.generic_noisemodel import get_noisemodelcat


@remote_data
def test_get_noisemodelcat():
    noise_fname = download_rename("beast_example_phat_noisemodel.grid.hd5")
    ntable = get_noisemodelcat(noise_fname)

    # check that at least the 3 basic elements are included
    expected_elements = ["error", "bias", "completeness"]
    for cexp in expected_elements:
        assert cexp in ntable.keys(), f"{cexp} values not found in noisemodel"


if __name__ == "__main__":
    test_get_noisemodelcat()
