import h5py
import numpy as np
import pytest
import tempfile

from beast.observationmodel.noisemodel.splinter import make_splinter_noise_model
from beast.physicsmodel.grid import SpectralGrid


@pytest.mark.parametrize("frac_unc", [0.01, 0.1, 0.3])
def test_splinter_noisemodel(frac_unc):

    # make super simplified model SED grid
    lamb = np.linspace(1000.0, 4000, 4)
    seds = np.logspace(-4, -3, 4)[None, :] * np.array([1, 1.5])[:, None]

    modelsedgrid = SpectralGrid(
        lamb=lamb, seds=seds, grid=[1], backend="memory"  # dummy input for now
    )

    # make splinter noisemodel
    noise_fname = tempfile.NamedTemporaryFile(suffix=f"{frac_unc}.grid.hd5").name
    # noise_fname = "/tmp/splinter_example_noisemodel_{:.2f}.grid.hd5".format(frac_unc)
    make_splinter_noise_model(noise_fname, modelsedgrid, frac_unc=frac_unc)

    # read entire noisemodel back in
    noisemodel = h5py.File(noise_fname, "r")

    # read the estimated sigma and check if close to manually computed errors
    sigma = noisemodel["error"]

    np.testing.assert_allclose(sigma, frac_unc * seds)
