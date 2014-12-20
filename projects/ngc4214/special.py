import numpy as np
from beast.core.phot import Filter
import pyfits


""" Provide lambda points to reinterpolate the filters over.
While filters in filters.hd5 or vega.hd5 are loaded via the fucntion 'load_all_filters' 
in the beast/core/phot.py, the filters defined here are not loaded via the function.
So there is no reinterpolaraion process, which is done by '__load__' function called in 'load_all_filters', 
during SED extraction for these filters. Have to reinterpolate separately here over the stellar model spectrum 
wavelength points. Default is Tlusty stellar model.  
"""
tmp = pyfits.getdata('/local/tmp/spiral/uv/beast/libs/tlusty.lowres.grid.fits',0)
wavelength = tmp[len(tmp)-1]


def make_integration_filter(startlam, endlam, dlamb, name,
                        observatory='SUDO',
                        instrument='FAKE',
                        comment=None):
    """
    Creates a constant-filter with transmission 100%
    in the energy range [startlam:endlam]

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
    ifT = np.interp(wavelength, lam, flam, left=0., right=0.)

    _name = '{obs:s}_{inst:s}_{name:s}'.format(obs=observatory, inst=instrument, name=name)
    filt = Filter(wavelength, ifT, name=_name)
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
    in the energy range [startlam:endlam]

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
    ifT = np.interp(wavelength, lam, flam_ph, left=0., right=0.)

    _name = '{obs:s}_{inst:s}_{name:s}'.format(obs=observatory, inst=instrument, name=name)
    filt = Filter(wavelength, ifT, name=_name)
    filt.bandwidth = bandwidth
    filt.comment = 'Multiply values by bandwidth'
    if comment is not None:
        filt.comment += '\n' + comment
    return filt


def test():
    f = make_integration_filter(90., 913., 1, 'QION', observatory='SUDO', instrument='Fake')
    return f
