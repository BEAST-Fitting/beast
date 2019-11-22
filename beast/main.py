import inspect
import argparse

from beast.tools import get_libfiles
from beast.plotting import plot_cmd, plot_filters

# cannot get plot_filters to work as the main parameter passed (filter_names)
# needs to allow multiple strings
#
# plot_indiv_fit needs work to make it callable with defaults this way


def main():
    """
    Main script for command-line use

    Summary: gets arguments for each of the scripts listed in 'scripts'
    and parses them for the function given in the input
    """
    all_funcs = []
    scripts = [get_libfiles, plot_cmd, plot_filters]  # scripts available

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        dest="subparser_name", help="sub-command help"
    )

    for item in scripts:
        scriptname = item.__name__.split('.')[-1]
        funcs_to_subcommand = inspect.getmembers(item, inspect.isfunction)
        for cfunc in funcs_to_subcommand:
            # only add the function that has the same name as the script file
            if cfunc[0] == scriptname:
                all_funcs.append(cfunc)

    # adds arguments for all of the functions listed in scripts
    for name, func in all_funcs:
        subparser = subparsers.add_parser(name, help=func.__doc__)
        for parname, arg in inspect.signature(func).parameters.items():
            sanitized_name = parname  # .replace("_", "-")
            if arg.default == inspect.Signature.empty:
                subparser.add_argument(sanitized_name)
            else:
                subparser.add_argument("--" + sanitized_name, default=arg.default)

    # now actually parses the arguments
    args = parser.parse_args()

    # and call the function
    for name, func in all_funcs:
        if args.subparser_name == name:
            # make a dictionary containing the correct inputs to the function,
            # extracted from the parsed arguments
            funcargs = {}
            for parname in inspect.signature(func).parameters:
                funcargs[parname] = getattr(args, parname)  # .replace("_", "-"))
                # funcargs[parname] = getattr(args, parname)
            func(**funcargs)
            break  # drop out immediately, which skips the "else" below
    else:
        raise AssertionError("Invalid subparser! This should be impossible...")
