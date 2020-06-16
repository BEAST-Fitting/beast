#!/usr/bin/env python

# note that other needed imports will be in the settings file
from astropy import units
from beast.tools import verify_beast_settings


class beast_settings:
    """
    All of the settings for running the BEAST

    Attributes
    ----------
    beast_settings_file : string
        name of the parameter file that will be used for the input settings

    """

    def __init__(self, input_settings_file):
        """
        Parameters
        ----------
        input_settings_file : string
            Name of the input file with BEAST settings

        """

        self.settings_file = input_settings_file

        if self.settings_file is not None:
            # read in the file
            self.read_beast_settings()
            # verify parameters
            self.verify_settings()


    def read_beast_settings(self):
        """
        Read in the beast settings file and set parameters
        """

        # read everything in as strings
        with open(self.settings_file, "r") as f:
            temp_data = f.readlines()
        # remove empty lines and comments
        input_data = [
            line.strip()
            for line in temp_data
            if line.strip() != "" and line.strip()[0] != "#"
        ]
        # if parameters are defined over multiple lines, combine lines
        for i in reversed(range(len(input_data))):
            if ('import ' not in input_data[i]) and ("=" not in input_data[i]):
                input_data[i - 1] += input_data[i]
                del input_data[i]

        # list of parameters that have units
        # (and therefore need special treatment)
        params_with_units = ["distance_unit", "velocity"]

        # parse it into a dictionary
        beast_params = {}

        for i in range(len(input_data)):
            # execute imports
            if 'import ' in input_data[i]:
                exec(input_data[i])

            # extract parameter and value (as strings)
            else:
                param = input_data[i].split("=")[0].strip()
                val = input_data[i].split("=")[1].strip()
                # if it's a unit, put into Quantity
                if param in params_with_units:
                    try:
                        beast_params[param] = units.Quantity(val)
                    except:
                        beast_params[param] = eval(val)
                # if it's not a unit
                else:
                    # exec the string to get parameter values accessible
                    exec(input_data[i])
                    beast_params[param] = eval(param)

        # turn dictionary into attributes
        for key in beast_params:
            setattr(self, key, beast_params[key])


    def verify_settings(self):
        """
        Run the verification code on beast_settings
        """

        verify_beast_settings.verify_input_format(self)
