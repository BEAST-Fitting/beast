"""
Base for trunchen noise model
Mainly a copy of the generic noisemodel.py
  but with the Table support provided by astropy.table

.. history::
    Started 18 Dec 2015 by Karl D. Gordon
"""

from astropy.table import Table

class NoiseModel(object):
    """ Initial class of noise models """
    def __init__(self, astfile, *args, **kwargs):
        self._aliases     = dict()
        self.data         = None

        self.astfile = astfile
        self.load_data()

    def __call__(self, sedgrid):
        raise NotImplemented

    def load_data(self):
        self.data = Table.read(self.astfile)

