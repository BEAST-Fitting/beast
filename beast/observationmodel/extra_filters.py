import numpy as np
from phot import IntegrationFilter
import astropy.io.fits as pyfits



def make_integration_filter(startlam, endlam, dlamb, name,
                        observatory='SUDO',
                        instrument='FAKE',
                        comment=None):
    """
    Creates a constant-filter with transmission 100%
    in the energy range [startlam:endlam].
    Calculate the total number of photons under this filter.

    Parmeters
    ---------
    startlam: float, optional (default=0.)
        lower limit of integration, in AA

    endlam: float, optional (default=912.)
        upper limit of integration, in AA

    dlam: ndarray
        wavelength resolution for definition, assuming in AA

    name: string
        name of the filter

    observatory: string, (default='PSEUDO')
        name of the observatory

    instrument: string (default='PSEUDO')
        name of the instrument

    comment: string (default=None)
        comment to the filter

    Returns
    -------
    filt: Filter instance
        filter object
    """
    #create a flat 100% transmission constant filter
    #per wavelength between startlam and endlam
    lam = np.arange(startlam, endlam, dlamb)

    flam = np.ones_like(lam)
    mask_i = lam < startlam
    mask_f = lam > endlam
    flam[mask_i] = 0.
    flam[mask_f] = 0.

    #Check what are actual boundaries, based on lam grid
    lam0 = lam[np.invert(mask_i)][0]
    lam1 = lam[np.invert(mask_f)][-1]

    bandwidth = lam1 - lam0

    _name = '{obs:s}_{inst:s}_{name:s}'.format(obs=observatory, inst=instrument, name=name)
    filt = IntegrationFilter(lam, flam, name=_name)
    filt.bandwidth = bandwidth
    filt.comment = 'Multiply values by integral{lamb*dlamb}/hc'
    if comment is not None:
        filt.comment += '\n' + comment
    return filt


def make_top_hat_filter(startlam, endlam, dlamb, name,
                        observatory='SUDO',
                        instrument='FAKE',
                        comment=None):
    """
    Creates a psudo-filter with transmission 100%
    in the energy range [startlam:endlam].
    

    Parmeters
    ---------
    startlam: float, optional (default=0.)
        lower limit of integration, in AA

    endlam: float, optional (default=912.)
        upper limit of integration, in AA

    dlam: ndarray
        wavelength resolution for definition, assuming in AA

    name: string
        name of the filter

    observatory: string, (default='PSEUDO')
        name of the observatory

    instrument: string (default='PSEUDO')
        name of the instrument

    comment: string (default=None)
        comment to the filter

    Returns
    -------
    filt: Filter instance
        filter object
    """
    #create a flat 100% transmission pseudo filter
    #per wavelength between startlam and endlam
    lam = np.arange(startlam, endlam, dlamb)

    flam = np.ones_like(lam)
    mask_i = lam < startlam
    mask_f = lam > endlam
    flam[mask_i] = 0.
    flam[mask_f] = 0.

    #Check what are actual boundaries, based on lam grid
    lam0 = lam[np.invert(mask_i)][0]
    lam1 = lam[np.invert(mask_f)][-1]

    # the output from getSEDs must be multiplied by bandwidth,
    # in order to get integrated flux between startlam and endlam
    bandwidth = lam1 - lam0

    # adjust filter form to actual filter transmission per photons
    # make an effectively perfect cuts
    # integral ( lambda T F dlambda ) = integral (lambda G / lambda F dlambda)
    flam_ph = flam / lam

    _name = '{obs:s}_{inst:s}_{name:s}'.format(obs=observatory, inst=instrument, name=name)
    filt = Filter(lam, flam_ph, name=_name)
    filt.bandwidth = bandwidth
    filt.comment = 'Multiply values by bandwidth'
    if comment is not None:
        filt.comment += '\n' + comment
    return filt


def test():
    f = make_integration_filter(90., 913., 1, 'QION', observatory='SUDO', instrument='Fake')
    return f
