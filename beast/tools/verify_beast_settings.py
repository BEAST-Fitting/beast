#!/usr/bin/env python
"""
Script to verify wrong format of the BEAST input parameters and terminate run_beast
It must be called run_beast.py
Created by Maria Kapala on Feb 24th 2017
"""
from os.path import exists
from numpy import inf
import warnings


def verify_range(param, param_name, param_lim):
    # check if input param limits make sense

    if any(p < param_lim[0] for p in param):
        raise ValueError(param_name + " min value not physical.")

    if any(p > param_lim[1] for p in param):
        raise ValueError(param_name + " max value not physical.")


def check_grid(param, param_name, param_lim):
    # check if input param limits and grid initialisation make sense

    param_min, param_max, param_step = param[0:3]

    if param_min < param_lim[0]:
        raise ValueError(param_name + " min value not physical.")

    if param_max > param_lim[1]:
        raise ValueError(param_name + " max value not physical.")

    if param_min > param_max:
        raise ValueError(param_name + " min value greater than max")

    if (param_max - param_min) < param_step:
        if param_max - param_min == 0.0:
            warnings.warn("Note: " + param_name + " grid is single-valued.")
        else:
            raise ValueError(param_name + " step value greater than (max-min)")


def verify_one_input_format(param, param_name, param_format, param_lim):
    """
    Test for an input parameter correctness of format and limits.

    Parameters
    ----------
    param: str, float, or list(float)
        Input parameter for run_beast.py, defined in beast settings

    param_name: str
        Extact name of the param, used for printing purposes.

    param_format: str

    param_lim: list(float) or None
       pass information about any physical limitations of the param;
       [-inf, inf] when param is not constraint at all;
       None when concerns a str parameter

    """

    if "list" in param_format:
        if not isinstance(param, list):
            if param is None:
                warnings.warn(param_name + " is not defined.")
            else:
                raise TypeError(param_name + " is not in the right format - a list.")
        elif "float" in param_format:
            if len(param) > 3:
                tparam = param[0:3]
                if isinstance(param[3], str):
                    if param[3] not in ["log", "lin"]:
                        raise ValueError(f"4th element in {param_name} is not log or lin")
            else:
                tparam = param
            is_list_of_floats = all(isinstance(item, float) for item in tparam)
            if not is_list_of_floats:
                raise TypeError(
                    param_name + " is not in the right format - list of floats."
                )
            elif "grid" in param_format:
                # when param is aranged from given [min, max, step],
                # instead of a specific list of values
                check_grid(param, param_name, param_lim)
            else:
                verify_range(param, param_name, param_lim)

        if "str" in param_format:
            if not isinstance(param, str):
                raise TypeError(param_name + " is not in the right format - a string.")
            elif "file" in param_format:
                if not exists(param):
                    raise OSError(
                        param_name + " does not exist. Please provide the file path."
                    )

        if "version" in param_format:
            if not isinstance(param, float):
                raise TypeError(param_name + " is not in the right format - a float")
            elif param not in param_lim:
                raise TypeError(
                    param_name
                    + " is an invalid number, leading to version of the isochrone."
                )


def verify_input_format(settings):

    """
    Define relevant parameters, their correct names, format and limits.
    Call verify_one_input_format to test for correctness of format and
    limits.

    Parameters
    ----------
    settings: beast.tools.beast_settings.beast_settings instance
        Input parameters are initialized in beast_settings

    """

    try:
        if settings.allow_warnings:
            print("verify_input_format: using non-interrupting warnings")
        else:
            warnings.simplefilter("error", UserWarning)
    except AttributeError:
        warnings.simplefilter("error", UserWarning)

    parameters = [
        settings.z,
        settings.obsfile,
        settings.astfile,
        settings.logt,
        settings.avs,
        settings.rvs,
        settings.fAs,
    ]
    parameters_names = ["z", "obsfile", "astfile", "logt", "avs", "rvs", "fAs"]
    param_format = [
        "list_float",
        "str_file",
        "str_file",
        "list_float_grid",
        "list_float_grid",
        "list_float_grid",
        "list_float_grid",
    ]

    print(settings.oiso.name)
    if settings.oiso.name == "MESA/MIST isochrones":
        print('Working on the MIST isochrone')
        parameters_limits = [
            [0.0142E-4, 0.0142 * 10**(0.5)],
            None,
            None,
            [5, 10.3],
            [0.0, inf],
            [1.0, 7.0],
            [0.0, 1.0],
        ]
    if settings.oiso.name == "Padova CMD isochrones":
        print('Working on the PARSEC isochrone')
        parameters_limits = [
            [1E-4, 0.06],
            None,
            None,
            [-inf, 10.15],
            [0.0, inf],
            [1.0, 7.0],
            [0.0, 1.0],
        ]
    if settings.oiso.name == "MIST EvolTracks":
        print('Working on the MIST Evolutionary Tracks')
        parameters[3] = settings.logmass
        parameters_names[3] = "logmass"

        parameters_limits = [
            settings.oiso.z_range,
            None,
            None,
            settings.oiso.logmass_range,
            [0.0, inf],
            [1.0, 7.0],
            [0.0, 1.0],
        ]
    else:
        print("setup needed for ", settings.oiso.name)
        exit()

    for i, param_ in enumerate(parameters):
        verify_one_input_format(
            param_, parameters_names[i], param_format[i], parameters_limits[i],
        )


if __name__ == "__main__":  # pragma: no cover

    pass
