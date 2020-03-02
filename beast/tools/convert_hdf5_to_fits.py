import numpy as np
from astropy.table import Table
import tables


def st_file(file_name):
    """Converts HDF5 files (specifically from Ben Williams) into FITS format.

    Parameters
    ----------
    file_name : str (default=None)
        Name of HDF5 file (will also be name of output file, with different extension)

    Returns
    -------
    None
        Saves object to new file

    """

    # Check to see if the input file has been specified correctly
    if file_name is None:
        print("Need a valid filename")
        return

    # read in the file
    print("reading in table")
    f = tables.open_file(file_name)

    # grab the relevant table
    print("grabbing data")
    tab = f.root.data.table.read()

    # names of the data set thingies
    set0 = f.root.data.table.attrs.values_block_0_kind
    set1 = f.root.data.table.attrs.values_block_1_kind
    set2 = f.root.data.table.attrs.values_block_2_kind
    print("data info:")
    print(set0, set1, set2)

    # make a dictionary to hold the data
    data_dict = {}

    # put labels in it
    for s in set0:
        data_dict[s] = []
    for s in set1:
        data_dict[s] = []
    for s in set2:
        data_dict[s] = []

    # iterate through the table to fill the dictionary
    print("filling up the dictionary with data")

    for ind in tab:

        # 0th item is the index
        # (don't need it)

        # 1st item is block 0
        for t, tag in enumerate(set0):
            data_dict[tag].append(ind[1][t])

        # 2nd item is block 1
        for t, tag in enumerate(set1):
            data_dict[tag].append(ind[2][t])

        # 3rd item is block 2
        for t, tag in enumerate(set2):
            data_dict[tag].append(ind[3][t])

    # save everything to astropy table in a fits file
    print("making astropy table")

    new_table = Table()

    for tag in set0:
        new_table[tag.upper()] = np.array(data_dict[tag])
    for tag in set1:
        new_table[tag.upper()] = np.array(data_dict[tag])
    for tag in set2:
        new_table[tag.upper()] = np.array(data_dict[tag])

    print("saving to fits file")
    new_table.write(
        file_name.replace("phot.hdf5", "st.fits"), format="fits", overwrite=True
    )


if __name__ == "__main__":
    st_file(file_name=None)
