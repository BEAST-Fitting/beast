from astropy.tests.helper import remote_data
from beast.tools import convert_hdf5_to_fits


@remote_data
def test_convert_hd5_to_fits():

    # Pick some random .hd5 file to convert
    data_fname = download_rename("beast_example_phat_seds_trim.grid.hd5")

    # Convert the file
    convert_hd5_to_fits.st_file(data_fname)
