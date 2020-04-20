from beast.tools.compare_spec_type import compare_spec_type
from beast.tests.helpers import download_rename, compare_tables
from astropy.tests.helper import remote_data
from astropy.table import Table

@remote_data
def test_compare_spec_type_inFOV():
    """
    Test for compare_spec_type.  The spectrally-typed stars aren't real sources,
    they're just invented for the purposes of documenting/testing the code.

    In this version, the stars are in the imaging field of view.
    """

    # download the needed files
    obs_fname = download_rename("b15_4band_det_27_A.fits")
    stats_fname = download_rename("beast_example_phat_stats.fits")

    # run compare_spec_type
    spec_type = compare_spec_type(
        obs_fname,
        stats_fname,
        [11.2335881, 11.23342557],  # RA
        [41.9001895, 41.90006316],  # Dec
        ['A', 'G'],                 # Spectral type
        [2, 7],                     # Subtype
        ['II', 'II'],               # Luminosity class
        match_radius=0.2            # Match radius (arcsec)
    )

    # expected output table
    expected_table = Table({
        'spec_ra': [11.2335881, 11.23342557],
        'spec_dec': [41.9001895, 41.90006316],
        'spec_type': ['A 2 II', 'G 7 II'],
        'spec_teff': [9000.0, 4916.666666666667],
        'spec_logg': [2.7164474106543732, 1.7184474106543735],
        'phot_cat_ind': [27, 8],
        'stats_cat_ind': [27, 8],
        'beast_teff_p50': [9046.250020338754, 4528.230977991138],
        'beast_teff_p16': [8643.670633196869, 4335.617282355577],
        'beast_teff_p84': [9536.391362054928, 4729.401710221546],
        'beast_logg_p50': [2.714286917261312, 1.7684285714285717],
        'beast_logg_p16': [2.636272525730954, 1.7014832653061227],
        'beast_logg_p84': [2.799534708811963, 1.8353738775510207],
        'teff_sigma': [-0.11488422362383206, 1.9308757510045778],
        'logg_sigma': [0.025343687546173433, -0.7465969411324851]
    })

    # compare to new table
    compare_tables(expected_table, Table(spec_type), rtol=2e-3)



@remote_data
def test_compare_spec_type_notFOV():
    """
    Test for compare_spec_type.  The spectrally-typed stars aren't real sources,
    they're just invented for the purposes of documenting/testing the code.

    In this version, the stars are NOT in the imaging field of view.
    """

    # download the needed files
    obs_fname = download_rename("b15_4band_det_27_A.fits")
    stats_fname = download_rename("beast_example_phat_stats.fits")

    # run compare_spec_type
    spec_type = compare_spec_type(
        obs_fname,
        stats_fname,
        [1.0],              # RA
        [1.0],              # Dec
        ['B'],              # Spectral type
        [4],                # Subtype
        ['V'],              # Luminosity class
        match_radius=0.2    # Match radius (arcsec)
    )

    # expected output table
    expected_table = Table({
        'spec_ra': [1.0],
        'spec_dec': [1.0],
        'spec_type': ['B 4 V'],
        'spec_teff': [None],
        'spec_logg': [None],
        'phot_cat_ind': [None],
        'stats_cat_ind': [None],
        'beast_teff_p50': [None],
        'beast_teff_p16': [None],
        'beast_teff_p84': [None],
        'beast_logg_p50': [None],
        'beast_logg_p16': [None],
        'beast_logg_p84': [None],
        'teff_sigma': [None],
        'logg_sigma': [None],
    })

    # compare to new table
    compare_tables(expected_table, Table(spec_type))
