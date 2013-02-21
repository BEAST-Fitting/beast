import numpy as np
from sedfitter.observations import Observations


#===============================================================================
# HANDLE THE OBSERVATION DATA
#   *real* or *fake* we need to handle the inputs in a same manner
#   so that it becomes a minimal effort to change the dataset
#
#
# Class FakeObs is derived from observations.Observations and aims at being
# as general as possible
#
# Real data will only require to define something very similar (see phat.py)
#==============================================================================
class FakeObs(Observations):
    """ A quick insertion of the class that will be used eventually
        so that all the code can be used as a template for real data

        * This class reads either a file or a variable
          thanks to eztable.Table construtor flexibility.

        * Instances can be saved with calling .writeto()
          (it uses Table.write)
          TODO: save the mask and filter names
    """
    def __init__(self, fname):
        """ Generate a data interface object """
        self.filters   = None
        self.desc      = None
        self.badvalue  = None
        self.setDistanceModulus(0.)
        self.readData(fname)

    def getMags(self, num, filters):
        """ TODO update """
        return np.array([ self.data[tt][num] for tt in filters])

    def getErrors(self, num, filters):
        """ TODO update """
        return np.array([ self.data[tt + 'err'][num] for tt in filters])

    def getObs(self, num=0, err=0.05):
        #assert ( self.filters is not None), "No filter set."
        mags = self.getMags(num, self.filters)
        errs = np.ones(len(mags), dtype=float) * err
        if self.badvalue is not None:
            mask = (mags >= self.badvalue)
        else:
            mask = np.zeros(len(mags), dtype=bool)

        return mags, errs, mask

    def readData(self, fname):
        """ read the dataset from the original source file """
        from eztables import Table
        if type(fname) == str:
            self.inputFile = fname
        else:
            self.inputFile = None
        self.data = Table(fname)

    def writeto(self, *args, **kwargs):
        self.data.write(*args, **kwargs)
