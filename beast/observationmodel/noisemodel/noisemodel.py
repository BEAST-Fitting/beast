from ...external.eztables import Table

__all__ = ['NoiseModel']

class NoiseModel(object):
    """ Initial class of noise models """
    def __init__(self, astfile, *args, **kwargs):
        self.astfile = astfile
        self.load_data()

    def load_data(self):
        self.data = Table(self.astfile)

    def fit(self, *args, **kwargs):
        raise NotImplemented

    def __call__(self, sedgrid):
        raise NotImplemented
