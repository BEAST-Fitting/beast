import numpy as np
from astropy.table import Table
import tables


def st_file(file_name, autotrim=False):
    """Converts HDF5 files (specifically from Ben Williams) into FITS format.

    Parameters
    ----------
    file_name : str (default=None)
        Name of HDF5 file (will also be name of output file, with .phot.hdf5
        replaced with .st.fits)
    autotrim : bool, optional (default=False)
        Whether to automatically trim column names that are too long for FITS files,
        and will therefore make this function crash otherwise

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

    # trim too-long column names, if requested, to stop the write function crashing
    if autotrim:
        set0_trim, set1_trim, set2_trim = replace_largest_common_substring(
            set0, set1, set2
        )
        for t, tag in enumerate(set0):
            print("Renaming 0 " + tag.upper() + " to " + set0_trim[t].upper())
            new_table.rename_column(tag.upper(), set0_trim[t].upper())
        for t, tag in enumerate(set1):
            print("Renaming 1 " + tag.upper() + " to " + set1_trim[t].upper())
            new_table.rename_column(tag.upper(), set1_trim[t].upper())
        for t, tag in enumerate(set2):
            print("Renaming 2 " + tag.upper() + " to " + set2_trim[t].upper())
            new_table.rename_column(tag.upper(), set2_trim[t].upper())

    # write to file
    print("saving to fits file")
    new_table.write(file_name.replace(".hdf5", ".fits"), format="fits", overwrite=True)
    f.close()


def replace_largest_common_substring(set0, set1, set2):
    """Finds the largest contiguous substring common to column names too long
    for FITS key names

    Parameters
    ----------
    set0 : list
        set0, as produced by st_file
    set1 : list
        set1, as produced by st_file
    set2 : list
        set2, as produced by st_file

    Returns
    -------
    set0 : list
        set0, with any necessary tag trimming to allow FITS table to write
    set1 : list
        set1, with any necessary tag trimming to allow FITS table to write
    set2 : list
        set2, with any necessary tag trimming to allow FITS table to write
    """

    # Find all column names longer than 65 characters
    sets = set0 + set1 + set2
    indices = np.where(np.array([len(tag) for tag in sets]) > 65)[0]

    # If no column names are too long, returns unchanged
    if len(indices) == 0:
        return set0, set1, set2

    # Extract too-long column names
    subset = [sets[i] for i in indices]

    # Check the shortest possible shared string
    shortest_str = min(subset, key=len)
    max_len = len(shortest_str)

    # Check substrings of decreasing length; break out when longest found
    longest = False
    for length in range(max_len, 0, -1):
        for start in range(0, max_len - length + 1):
            candidate = shortest_str[start : start + length]
            if all(candidate in s for s in subset):
                longest = True
                longest_str = candidate
                break
        if longest:
            break

    # Raise exception if no fix possible
    if not longest:
        raise Exception(
            "Some column names are too long to write to FITS table, "
            + "but no shared string found to trim"
        )

    # Find maximum length of column names, and therefore how much trimming is necessary
    longest_tag = np.array([len(tag) for tag in sets]).max()
    n_char_trim = longest_tag - 65

    # If longest shared string isn't enough, raise exception
    if (longest_tag - n_char_trim) >= 68:
        raise Exception(
            "No shared string in column names is long enough to enable "
            + "enough trimming to get column names short enough to allow "
            + "FITS table to write without crashing"
        )

    # Use this to decide what replacement stirng should be, first checking if we can just replace - and _ characters
    new_str_test = longest_str.replace("-", "").replace("_", "")
    if (len(longest_str) - len(new_str_test)) >= n_char_trim:
        new_str = new_str_test
    else:
        new_str = longest_str[n_char_trim:]

    # Replace all colnames that contain our long shared string with the trimmed verison
    set0_trim, set1_trim, set2_trim = set0.copy(), set1.copy(), set2.copy()
    for i in range(len(set0)):
        set0_trim[i] = set0[i].replace(longest_str, new_str)
    for i in range(len(set1)):
        set1_trim[i] = set1[i].replace(longest_str, new_str)
    for i in range(len(set2)):
        set2_trim[i] = set2[i].replace(longest_str, new_str)

    # Return trimed column names
    return set0_trim, set1_trim, set2_trim


if __name__ == "__main__":
    st_file(file_name=None)
