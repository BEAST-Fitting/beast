""" Latex export """

import os, inspect, sys
localpath = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])
import numpy as np
from .basebackend import BaseBackend
from ..core.tableheader import TableHeader
from ..table import Table

class LatexBackend(BaseBackend):
	def __init__(self):
		""" constructor """
		BaseBackend.__init__(self, tableType='tex')

	def read(self, filename, *args, **kwargs):
		raise Exception('Latex reading not implemented')


	def writeData(self, unit, data, fmt, delimiter='&'):
		""" Write data part into the opened unit """
		size = data.shape[0]
		for ik in range(size):
			unit.write( fmt % data[ik].tolist() )
			unit.write("\\\\\n")

	def write(self, tab, output='exportedData.tex', header=True,
			comments='%', verbose=False, **kwargs):
		"""
		export data to a latex Table

		inputs:
			data -- data dictionnary to export
		
		outputs:
			output -- output file (def: exportedData.dat)

		keywords:
			header    -- header object describing the file
			delimiter -- delimiter to use (def: ',')
			comment   -- comment character for header (def: '#')
			keep      -- keeps unit opened
			unit	  -- uses opened stream, if provided 
		"""
		
		if hasattr(output, 'write'):
			unit = output
		else:
			unit = open(output, 'w')

		unit.write('\\begin{table}\n \\begin{center}\n')
		if 'NAME' in tab.header:
			if tab.header['NAME'] not in ['', None, 'None']:
				unit.write('\\caption{%s}\n' % tab.header['NAME'])

		aligntxt = ''.join(['c']*tab.ncols)

		unit.write('\\begin{tabular}{%s}\n' % aligntxt)

		# write tab header
		unit.write('\\hline\\hline\n')
		colDefTxt = ''	
		_Notes = []
		for k in tab.columns:
			colDefTxt += ' & %s' % k
			if tab.columns[k].description not in ['', 'None', None]:
				_Notes.append(tab.columns[k].description)
				colDefTxt += "$^{\\rm{(%s)}}$" % len(_Notes)
		colDefTxt = colDefTxt[2:]
		unit.write('%s\\\\ \n' % colDefTxt )

		units = [ tab.columns[k].unit for k in list(tab.keys()) ]

		#Add units if present
		if ''.join(units).replace('None','') != '':
			units = ' & '.join(units).replace('None', '')
			unit.write('%s \\\\\n' % units)
		unit.write('\\hline\n')

		fmt  = ' & '.join(['%'+tab.columns[k].format for k in tab.columns])
		self.writeData(unit, tab.data, fmt, delimiter=' & ')	

		unit.write('\\hline\n')

		# add column coments
		if len(_Notes) > 0:
			unit.write('\\begin{scriptsize}\n')
			for k in range(len(_Notes)):
				unit.write("$^{\\rm(%d)}$ %s\\\\\n" % (k+1, _Notes[k]))
			unit.write('\\end{scriptsize}\n')
		unit.write('\\end{tabular}\n\\end{center}\n\\end{table}\n%end of table')

		if hasattr(output, 'write') :
			return unit
		else:
			unit.close()

		if verbose: print("Data exported into %s" % output)

