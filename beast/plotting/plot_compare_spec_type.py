import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table


def plot_compare_spec_type(
        match_file,
        one_2_one_line=True,
        savefig=False,
        plot_indiv_fit_kwargs=None
):
    """
    Plots the comparison of spectrally-typed stars with BEAST fits
    as calculated by tools/compare_spec_type.py.
    Optionally plots the 1D PDF of individual matched star.

    Note1: An easy addition would be to optionally plot the 1dpdfs of the
    matched stars with arguments supplied in plot_indiv_fit_kwargs.

    Note2: Currently compare_spec_type does not accept a range of spectral
    subtypes or luminosities (stemming from uncertainty in spectral typing),
    which would result in error bars for Teff and log(g) calculated for the
    specteally-typed star. When this is implemented, it would be good to plot
    these error bars as well.

    Parameters
    ----------
    match_file : string
        File containing the matching stats produced by compare_spec_type.py

    one_2_one_line : bool
        Optionally plot a 1-to-1 line for easy comparison

    savefig : str
        set to the file extension for the desired plot file (e.g., png, pdf, etc)

    plot_indiv_fit_kwargs : dict (default=None)
        If provided, plot the individual 1D PDF for the star(s) which show a match

    Returns
    -------

    """

    t = Table(fits.open(match_file)[1].data)

    teff_spec = t["spec_teff"][~np.isnan(t["spec_teff"])]
    teff_beast = t["beast_teff_p50"][~np.isnan(t["beast_teff_p50"])]
    teff_beast_1sig = (t["beast_teff_p84"][~np.isnan(t["beast_teff_p84"])]
                       - t["beast_teff_p16"][~np.isnan(t["beast_teff_p16"])]) / 2.

    logg_spec = t["spec_logg"][~np.isnan(t["spec_logg"])]
    logg_beast = t["beast_logg_p50"][~np.isnan(t["beast_logg_p50"])]
    logg_beast_1sig = (t["beast_logg_p84"][~np.isnan(t["beast_logg_p84"])]
                       - t["beast_logg_p16"][~np.isnan(t["beast_logg_p16"])]) / 2.

    fig, ax = plt.subplots(figsize=(10, 5))
    axt = plt.subplot(121)
    plt.errorbar(teff_beast, teff_spec, xerr=teff_beast_1sig, ms=3, fmt='o')
    plt.xlabel('$T_{eff}$ BEAST [K]', fontsize=14)
    plt.ylabel('$T_{eff}$ Lit. [K]', fontsize=14)
    txlim = axt.get_xlim()
    tylim = axt.get_ylim()

    axg = plt.subplot(122)
    plt.errorbar(logg_beast, logg_spec, xerr=logg_beast_1sig, ms=3, fmt='o')
    plt.ylabel('log(g) Lit. [cm $s^{-2}$]', fontsize=14)
    plt.xlabel('log(g) BEAST [cm $s^{-2}$]', fontsize=14)
    gxlim = axg.get_xlim()
    gylim = axg.get_ylim()

    if one_2_one_line:
        l_t = np.linspace(0, 40000, 2)
        l_g = np.linspace(0, 5, 2)
        axt.plot(l_t, l_t, ls='--', lw=1, c='gray')
        axg.plot(l_g, l_g, ls='--', lw=1, c='gray')

    axt.set_xlim(txlim)
    axt.set_ylim(tylim)
    axg.set_xlim(gxlim)
    axg.set_ylim(gylim)

    plt.subplots_adjust(wspace=0.25)

    if savefig:
        figname = match_file.replace('fits', savefig)
        fig.savefig('%s' % figname)
