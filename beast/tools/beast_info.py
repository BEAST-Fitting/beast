from os import path
import copy
import asdf

__all__ = ["add_to_beast_info_file"]


def add_to_beast_info_file(infoname, info):
    """
    Add information to the beast info file.

    Parameters
    ----------
    infoname : str
        the beast info filename

    info : dict
        dictionary to add to the beast info file
    """
    if path.exists(infoname):
        with asdf.open(infoname) as af:
            tree = copy.deepcopy(af.tree)
    else:
        tree = {}

    # update entries if they already exist, otherwise add them
    for ckey in info.keys():
        tree[ckey] = info[ckey]

    # Create the ASDF file object from our data tree
    af = asdf.AsdfFile(tree)

    # Write the data to a new file
    af.write_to(infoname)
