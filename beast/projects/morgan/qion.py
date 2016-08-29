import numpy as np
from beast import phot, grid


def make_filter():
    """ Make a virtual filter to extact ionizing flux under 912 AA """
    # top hat function from 0 to 912AA, 0 after
    qion_lamb = np.hstack([np.linspace(0, 912, 20), np.linspace(913, 1000, 10)])
    qion_pass = np.hstack([np.ones(20), np.zeros(10)])

    Qion_filter = phot.Filter(qion_lamb, qion_pass, name='Virtual_Qion')
    return Qion_filter


fname = 'mf_ngc4214_full_new/mf_ngc4214_full_spec.grid.hd5'

g = grid.SpectralGrid(fname, backend='cache')
Qion = make_filter()
cls, seds, grid = phot.extractSEDs(g, [Qion])

h = 6.62606957 * 1e-34   # m **2 * kg / s
c = 299.792458 * 1e6     # m / s

Eion = seds / (h * c)
