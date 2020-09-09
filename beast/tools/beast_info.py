from os import path
import asdf

__all__ = ["add_to_beast_info_file"]


def add_to_beast_info_file(project, info):
    """
    Add information to the beast info file.

    Parameters
    ----------
    project : str
        name of project that gives the beast info filename

    info : dict
        dictionary to add to the beast info file
    """
    infoname = f"{project}/{project}_beast_info.asdf"
    if path.exists(infoname):
        with asdf.open(infoname) as af:
            tree = af.tree
    else:
        tree = {}

    # update entries if they already exist, otherwise add them
    for ckey in info.keys():
        tree[ckey] = info[ckey]

    # Create the ASDF file object from our data tree
    af = asdf.AsdfFile(tree)

    # Write the data to a new file
    af.write_to(infoname)
