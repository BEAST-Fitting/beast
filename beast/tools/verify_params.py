#!/usr/bin/env python
"""
	Script to verify and warn about wrong format of the BEAST input parameters.
	It must be called run_beast.py 

"""

from warnings import warn
from os.path import exists


def verify_range(param, param_lim)
	# check if input param limits make sense
	param_min, param_max, param_step = param

	if param_min > param_max:
		warn(param_name + " min value greater than max")

	if (param_max-param_min) < param_step:
		warn(param_name + " step value greater than (max-min)")





def verify_input_format(param, param_name, param_lim=None):
  ''' Test for an input parameter correctness of format and limits (if provided)
  '''

  
  # check input parameter's type
  if type(param) is list:
  	is_list_of_floats = all(type(item) is float for item in my_list)
  	if not is_list_of_floats:
  		warn(param_name + " is not a right format - list of floats.")

  	elif param_lim not None:
  		# when param is aranged from given [min, max, step],
  		# instead of a specific list of values
  		verify_range(param, param_lim)

  elif type(param) is str:
  	if not exists(param):
  		warn(param_name + " does not exist. Please provide the file path.") 
  elif type(param) is float:





