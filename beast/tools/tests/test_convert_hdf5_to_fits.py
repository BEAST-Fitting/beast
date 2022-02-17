import pytest
from astropy.table import Table
from astropy.io import fits
from beast.tools.convert_hdf5_to_fits import st_file

from beast.tests.helpers import download_rename, compare_tables


@pytest.mark.remote_data
def test_convert_hd5_to_fits():

    # Pick some random .hd5 file to convert
    data_fname = download_rename("M31-B09-EAST_tinychunk.phot.hdf5")
    data_fname_cache = download_rename("M31-B09-EAST_tinychunk.st.fits")

    # Convert the file
    st_file(data_fname)

    # Compare the contents of the new file to the cached version
    data = Table(fits.getdata(data_fname.replace("phot.hdf5", "st.fits")))
    data_cache = Table(fits.getdata(data_fname_cache))

    compare_tables(data_cache, data)
