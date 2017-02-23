#!/usr/bin/env python
"""
	Script to verify and warn about wrong format of the BEAST input parameters.
	It must be called run_beast.py 

"""

from warnings import warn
from os.path import exists
from numpy import inf


def verify_range(param, param_name, param_lim):
	# check if input param limits make sense
	# no constarins indicated by param_lim = [-np.inf, np.inf]
	
	if any(p < param_lim[0] for p in param):
		warn(param_name + " min value not physical.")

	if any(p > param_lim[1] for p in param):
		warn(param_name + " max value not physical.")



def check_grid(param, param_name, param_lim):
	param_min, param_max, param_step = param

	# check if input param limits make sense
	# no constarins indicated by param_lim = [-np.inf, np.inf]
	if param_min < param_lim[0]:
		warn(param_name + " min value not physical.")

	if param_max > param_lim[1]:
		warn(param_name + " max value not physical.")

	if param_min > param_max:
		warn(param_name + " min value greater than max")

	if (param_max-param_min) < param_step:
		warn(param_name + " step value greater than (max-min)")



def verify_input_format(param, param_name, param_format, param_lim=None):
	''' Test for an input parameter correctness of format and limits (if provided)
	'''
	# check input parameter's format
	if 'list' in param_format:
	  	if type(param) is not list:
	  		warn(param_name + " is not in the right format - a list.")
	  	elif 'float' in param_format:
	  		is_list_of_floats = all(type(item) is float for item in param)
	  		if not is_list_of_floats:
	  			warn(param_name + " is not in the right format - list of floats.")
	  		elif 'grid' in param_format:
	  			# when param is aranged from given [min, max, step],
				# instead of a specific list of values
				check_grid(param, param_name, param_lim)
			else:
	  			verify_range(param, param_name, param_lim)
	  			


	if 'str' in param_format:
		if type(param) is not str:
			warn(param_name + " is not in the right format - a string.") 
		elif 'file' in param_format:
			if not exists(param):
				warn(param_name + " does not exist. Please provide the file path.") 


	if 'version' in param_format:
		if type(param) is not float:
			warn(param_name + " is not in the right format - a float")
  		elif param not in param_lim:
  			warn(param_name + " is an invalid version of the isochrone.")



if __name__ == "__main__":

 	
 	bright_limits_mag = [14., 14.5, 16., 15., 16., 14., 14.5, 14., 14.]
 	sens_limits_mag = [26., 26., 27., 29., 27.5, 28., 28.5, 27., 26.]
 	z = [0.03, 0.019, 0.008, 0.004]
	obsfile = '/Users/maria/Documents/data/beast/hack_week/todo_verify_params.txt'
	astfile = 'data/fake_stars_b15_27_all.hd5'
	#logt = [6.0, 10.13, 1.0]
	logt = [6.0, 10.13, 1.]
	avs = [0.0, 10.055, 1.0]
	rvs = [2.0,6.0,1.0]
	fbumps = [0.0,1.0, 0.25]
	trackVersion = 2.7

	parameters = [bright_limits_mag, sens_limits_mag, z, obsfile, astfile, logt, avs, rvs, fbumps, trackVersion]
	parameters_names = ['bright_limits_mag', 'sens_limits_mag', 'z', 'obsfile', 'astfile', 'logt', 'avs', 'rvs', 'fbumps', 'trackVersion']
	param_format = ['list_float', 'list_float', 'list_float', 'str_file', 'str_file', 'list_float_grid', 'list_float_grid', 'list_float_grid', 'list_float_grid', 'version']
	parameters_limits = [ [-inf, inf], [-inf, inf], [0., 0.1], None, None, [-inf, 10.15], [0., inf], [1., 7.], [0., 1.], [2.3, 2.7]]
	
	for i, param_ in enumerate(parameters):
		verify_input_format(param_, parameters_names[i], param_format[i], parameters_limits[i])

    


