# use evolutionary tracks instead of isochrones as the basis of
# the stellar physicsgrid

import numpy as np
from astropy.table import Table


__all__ = ['EvolTracks', 'ETParsec', 'ETMist']


class EvolTracks(object):
    """
    Stores the evolutionary tracks interpolated to input mass and
    metallicity grids with the age grid set by the evolutionary track
    calculations (may wish to revisit and allow for condensation).

    Attributes
    ----------
    logmass : log10 of the masses (solar masses)

    Z : numpy array of floats
        metallicites of tracks

    data: table
        columns include:
        M_ini - initial mass [Msun]
        M_act - actual mass [Msun]
        Z - metallicity [??]
        logL - log luminosity [Lsun]
        logg - log surface gravity [cm/s^2]
        logA - log age [years]
        logT - log surface effective temperature [K]
        stage - evolutionary stage [index]
    """
    def __init__(self):
        self.name = '<auto>'

        # define axis labels for plotting
        self.alabels = {'logT': 'log(Teff)',
                        'logg': 'log(g)',
                        'logL': 'log(L)',
                        'logA': 'log(age)',
                        'phase': 'evol phase',
                        'M_act': 'current mass',
                        'M_ini': 'initial mass'}

    def load_orig_tables(self, source):
        """
        Read the tracks from  the original files (stub)
        """
        print('not implemented')

    def plot_tracks(self, ax, xval='logT', yval='logL'):
        """
        Plot the tracks with the input x, y choices

        Parameters
        ----------
        ax : matplotlib axis
            where to plot

        xval : str, optional
            what data for x

        xval : str, optional
            what data for y
        """
        if xval not in self.data.keys():
            raise ValueError("xval choice not in data table")
        if yval not in self.data.keys():
            raise ValueError("yval choice not in data table")

        ax.plot(self.data[xval], self.data[yval])

        ax.set_xlabel(self.alabels[xval])
        ax.set_ylabel(self.alabels[yval])


class ETParsec(EvolTracks):
    """
    EvolTracks specific to PARSEC calculations
    """

    source = 'PARSEC'

    def load_orig_tables(self, files):
        """
        Read the tracks from  the original files

        files : str
            file or files with the evolutionary track calculations
            often each file is for a single initial mass
        """
        if not type(files) is list:
            files = [files]

        mass_act = np.array([])
        logL = np.array([])
        logT = np.array([])
        for cfile in files:
            a = Table.read(cfile, format='ascii')
            mass_act = np.concatenate((mass_act, a['MASS'].data))
            logL = np.concatenate((logL, a['LOG_L'].data))
            logT = np.concatenate((logT, a['LOG_TE'].data))

        self.data = {}
        self.data['M_act'] = np.array(mass_act)
        self.data['logL'] = np.array(logL)
        self.data['logT'] = np.array(logT)


class ETMist(EvolTracks):
    """
    MIST evolutionary tracks
    """
    source = 'MIST'

    def load_orig_tables(self, files):
        """
        Read the tracks from  the original files

        files : str
            file or files with the evolutionary track calculations
            often each file is for a single initial mass
        """
        if type(files) is not list:
            files = [files]

        mass_act = np.array([])
        logA = np.array([])
        logL = np.array([])
        logT = np.array([])
        logg = np.array([])
        phase = np.array([])
        eep = np.array([])
        for cfile in files:
            a = Table.read(cfile, format='ascii', header_start=11)
            mass_act = np.concatenate((mass_act, a['star_mass'].data))
            logA = np.concatenate((logA, np.log10(a['star_age'].data)))
            logL = np.concatenate((logL, a['log_L'].data))
            logT = np.concatenate((logT, a['log_Teff'].data))
            logg = np.concatenate((logg, a['log_g'].data))
            phase = np.concatenate((phase, a['phase'].data))
            eep = np.concatenate((eep, range(len(a))))

        self.data = {}
        self.data['M_act'] = mass_act
        self.data['M_ini'] = np.full((len(mass_act)), max(mass_act))
        self.data['logA'] = logA
        self.data['logL'] = logL
        self.data['logT'] = logT
        self.data['logg'] = logg
        self.data['phase'] = phase
        self.data['eep'] = eep

        self.alabels['eep'] = 'EEP'
