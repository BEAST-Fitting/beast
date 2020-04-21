from beast.tools import star_type_probability
from beast.tests.helpers import download_rename, compare_tables
from astropy.tests.helper import remote_data
from astropy.table import Table

@remote_data
def test_star_type_probability_all_params():
    """
    Test for star_type_probability.py
    """

    # download the needed files
    pdf1d_fname = download_rename("beast_example_phat_pdf1d.fits")
    pdf2d_fname = download_rename("beast_example_phat_pdf2d.fits")
    star_prob_fname = download_rename("beast_example_phat_startype.fits")

    # run star_type_probability
    star_prob = star_type_probability.star_type_probability(
        pdf1d_fname,
        pdf2d_fname,
        output_filebase = None,
        ext_O_star_params = {'min_M_ini':10, 'min_Av':0.5, 'max_Av':5}
    )

    # expected output table
    expected_star_prob = Table.read(star_prob_fname)

    # compare to new table
    compare_tables(expected_star_prob, Table(star_prob))

@remote_data
def test_star_type_probability_no_Av():
    """
    Test for star_type_probability.py
    """

    # download the needed files
    pdf1d_fname = download_rename("beast_example_phat_pdf1d.fits")
    pdf2d_fname = download_rename("beast_example_phat_pdf2d_no_Av.fits")
    star_prob_fname = download_rename("beast_example_phat_startype_no_Av.fits")


    # run star_type_probability
    star_prob = star_type_probability.star_type_probability(
        pdf1d_fname,
        pdf2d_fname,
        output_filebase = None,
        ext_O_star_params = {'min_M_ini':10, 'min_Av':0.5, 'max_Av':5}
    )

    # expected output table
    expected_star_prob = Table.read(star_prob_fname)

    # compare to new table
    compare_tables(expected_star_prob, Table(star_prob))
