ezTables v1.0dev
================

This module provides table manipulations regardless of storage format
Complete rewrite based on my old mytables package.


* Bidirectional support (read/write): 
        **csv, ascii**
        **fits**            (requires pyfits)
        **hdf5**            (requires pytables)
        **json**            (requires json)


* only write support: 
        **latex**


* starting to add astro tools

TODO
----

* add BLOSC compression to hdf5 tables: more load on cpu but less on disk IO 


* cone search on position keys without vo tools


* support more formats
	* sqlite
	* psql
	* mysql
	* lsd
	* votable


* distributed tables: chuncks of lines distributed over n-workers
	* MPI
	* BLOSC carray, https://github.com/FrancescAlted/carray
	* multiprocessing + zmq or shared memory
	* pure MapReduce ?? (backend mp or mpi)
	* look at Durus FileStorage: https://github.com/nascheme/durus
	(persistent storage -- pickleshare + network)
	* look at ZODB (not sure about network)
	(persistent storage -- pickleshare + network)


RANDOM IDEAS
------------

* consider the idea of a DB = aggregation of tables
	* need join -- done
	* need storage/reading/writing
	* can be mixed formats!!
	* need a query parser


* develop a pytable array based for heavy computations --> better than memmap
	* pytable is able to do fast computations even with disk access
	* optimized expressions using tables.expr (BLAST)
	* faster than numpy memmap, and numexpr
	* optimized where queries: 
	'''python
	>>> for row in tbl.where('(sqrt(x**2 + y**2) <= 1) & (z < 100)'):
	... do something with row['name'], row['x']...
	'''


REVISION HISTORY
-----------------

	0.1 -- 04/03/2011, MF

		Definition of the classes related to a Table object: 
				Table, 
				TableHeader,
				TableColumn,
				TableColumnHeader, 
				TableManager,
			Default registered formats:
				asciiManager (dat, txt, ascii) 
				csvManager   (csv)
				fitsManager  (fits) -- requires pyfits
			Import/Export manages also column description/comments & units
			in all the default registered formats
	0.2 -- 07/03/2011, MF
		New formats added:
				hd5Manager (hd5, hdf5) -- requires pytables
				latexManager (tex) only writing
	0.3 -- 07/03/2011, MF
		New feature added:
			Tables are now callable to produce simple operations
			you can now generate a column based on simple scalar operations
			or also +,-,*,/ operation.
			example:  > r = table('x+2*y')  or > r = table('a/b')
			but this is not tested to follow operation priorities, besides
			it could not perform complex operations sur as >table('x/(a+b)')
	0.3.1 -- 09/03/2011, MF
		Debug: Fits export, 
				order in the units, comments fixed
				string export fixed
				multiple comment & history keywords export fixed
	0.3.2 -- 01/04/2011, MF
		New features added
			Added join & match functions to the module
	0.3.3 -- 12/04/2011, MF
		Bug correction: opening HD5 tables that does not contain UNITS...
	0.3.4 -- 09/05/2011, MF
		Bug correction: opening HD5 tables and tablename issues leading
			to errors in reading table header
	0.3.5 -- 31/08/2011, MF
		Bug correction: ascii files issues
				delimiter was not always taken into account
				now you can force a given line to be the column
				description forceHealLine keyword
	0.3.6 -- 09/09/2011, MF
		Bug correction: ascii/csv files issues when mixing with text
				replacing loadfromtxt call by a recfrom<> call which
				does the job better
				if a column name already exists it will add a suffix
				instead of overwriting the column
	0.4 -- 09/09/2011, MF
		added features:
			disp        -- pretty print (part of) the table 
			__str__	    -- returns disp with default arguments
			evalexpr    -- let you do some simple operations on the table using
				       column names as variables (incl. math symbols) and
				       external variables as well
				       TODO: need to find a way to extend math functions to
				       arrays
			where       -- traditional where function based on evalexpr
			selectWhere -- is equivalent to Table[where()] and return a table
				       with only matching conditions (can also restrict 
				       selected	fields)
			extract     -- Returns a sub-table based on a given selection of
				       lines or fields
			__call__ updated to use evalexpr
			__getitem__ returns uses getRow for a single numerical value
					    uses getCol for a string
					    uses extract for an iterable
	0.4.1 -- 15/09/2011, MF
		added features:
			stats        -- returns a table with statistics over the
					selected columns (mean, std, hpd, quantiles...)
	0.5.0-- 13/10/2011, MF
		added features:
			distributedTable-- new table class that let you use disk space
					instead of RAM and allow multiple concurrent
					access
					(Still in development)
					Currently should work for reading
					Saving changes are currently not fully
					operational, data headers are not updated in
					realtime.
			Generating tables from recarray.
			Added a idl save file reader
	0.6.0 -- 30/03/2012, MF
		added features;
			sqltables are now managed through sqlManager nd sqlite3
			you can use it to load or generate sqlite tables. 

	1.0.dev -- 01/11/2012, MF
		Major changes:
			the entire code has been changed to be based on np.recArray data storage
				-- reduces memory footprint
				-- more transparent in the usage with numpy as an array
				-- Tables are pickleable objects (so far)
				-- evalexpr is replace by python eval using the table
				   dictionnary as a globals (incl. numpy) 
		Added features:
			-- the constructor of the class is able to read files directly
			   replacing the load() function
			-- JSON import and export is now available
			-- Aliases of columns is now supported (and kept through exports)
			-- Columns are now ordered and their positions can be changed
			-- The table can be generated from ndarray, dict or Table objects (copy)
			-- the constructor can handle streams as well as files (e.g.  StringIO)
		
		Removed features:
			-- load 	replaced by the class constructor >>> Table(file)
			-- extract	replaced by __getitem__ >>> tab[fields][indexes] or tab[indexes][fields]

