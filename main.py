import inspect
import pdb
import argparse
import sys

import beast
from beast.plotting import plot_cmd, plot_filters, plot_indiv_fit, beastplotlib
from beast.tools import get_libfiles


def main():
    """Summary gets arguments for each of the scripts listed in 'scripts' 
    and parses them for the function given in the input
    """
    all_funcs = []
    scripts = [plot_cmd, plot_filters, plot_indiv_fit]

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--global-setting", action="store_true", help="some global thingie"
    )

    subparsers = parser.add_subparsers(
        dest="subparser_name", required=True, help="sub-command help"
    )

    for item in scripts:
        funcs_to_subcommand = inspect.getmembers(item, inspect.isfunction)
        all_funcs.append(funcs_to_subcommand)
    all_funcs = [item for sublist in all_funcs for item in sublist]

    for name, func in all_funcs:
        subparser = subparsers.add_parser(name, help=func.__doc__)
        for parname, arg in inspect.signature(func).parameters.items():
            sanitized_name = parname.replace("_", "-")
            if arg.default == inspect.Signature.empty:
                subparser.add_argument(sanitized_name)
            else:
                subparser.add_argument("--" + sanitized_name, default=arg.default)

    # now actually parse the arguments
    args = parser.parse_args()

    # and call the function
    for name, func in all_funcs:
        if args.subparser_name == name:
            # make a dictionary containing the correct inputs to the function,
            # extracted from the parsed arguments
            funcargs = {}
            for parname in inspect.signature(func).parameters:
                funcargs[parname] = getattr(args, parname.replace("-", "_"))
            func(**funcargs)
            break  # drop out immediately, which skips the "else" below
    else:
        assert False, "Invalid subparser! This should be impossible..."


if __name__ == "__main__":
    main()
