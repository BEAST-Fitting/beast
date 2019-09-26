import inspect
import pdb
import argparse

import beast
from beast.plotting import plot_cmd, run_beast


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--global-setting', action='store_true',
                        help='some global thingie')

    subparsers = parser.add_subparsers(help='sub-command help')
    subparsers_plot = parser.add_subparsers(dest='plot', help='plot help')

    # create the subparsers and populate them with the correct arguments
    func_run = inspect.getmembers(run_beast, inspect.isfunction)
    func_plot = inspect.getmembers(plot_cmd, inspect.isfunction)

    def add_func(funcs_to_subcommand, subparsers):
        for name, func in funcs_to_subcommand:
            subparser = subparsers.add_parser(name, help=func.__doc__)
            for parname, arg in inspect.signature(func).parameters.items():
                # sanitized_name = parname.replace('_', '-')
                if arg.default == inspect.Signature.empty:
                    subparser.add_argument(parname)
                else:
                    subparser.add_argument('--' + parname,
                                           default=arg.default)
        return subparser

    add_func(func_run, subparsers)
    add_func(func_plot, subparsers)

    # now actually parse the arguments
    args = parser.parse_args()

    # set the global
    plot_cmd.GLOBAL_SETTING = args.global_setting

    # and call the function
    for name, func in funcs_to_subcommand:
        if args.subparser_name == name:
            # make a dictionary containing the correct inputs to the function,
            # extracted from the parsed arguments
            funcargs = {}
            for parname in inspect.signature(func).parameters:
                funcargs[parname] = getattr(args, parname)
            # pdb.set_trace()
            func(**funcargs)
            break  # drop out immediately, which skips the "else" below
    else:
        assert False, 'Invalid subparser! This should be impossible...'


if __name__ == '__main__':
    main()
