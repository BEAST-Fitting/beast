from astropy.table import Table
# from beast.external.eztables.table import Table

__all__ = ["NoiseModel"]


class NoiseModel(object):
    """ Initial class of noise models """

    def __init__(self, astfile, *args, **kwargs):
        self.astfile = astfile
        self.filter_aliases = {}
        self.load_data()

    def load_data(self):
        self.data = Table.read(self.astfile)

    def fit(self, *args, **kwargs):
        raise NotImplementedError

    def __call__(self, sedgrid):
        raise NotImplementedError
