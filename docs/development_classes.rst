*******
Classes
*******

The BEAST uses classes to abstract the specific implementation from the details
of storage in memory and on disk.

Grids
=====

The two main classes to use to access information in spectral and SED grids are:

- :class:`~beast.physicsmodel.grid.SEDGrid`:
  Provides the interface to grids of extinguished SEDs composed of band integrated fluxes.
  Used fairly extensively in the BEAST to interface with the core physicsmodel grid.

- :class:`~beast.physicsmodel.grid.SpectralGrid`:
  Provides the interface to grids of spectra.
  Generally used as part of generating the physics grid.
  Provides functions to apply dust extinction and compute the band integrated
  fluxes from the spectra.

The core attributes of both classes are the same.

- lamb: wavelengths of the seds.
- seds: fluxes at the wavelengths
- grid: astropy Table giving the parameters of each sed (mass, age, etc.)
- keys: column names of the grid (can be present w/o grid existing: e.g., backend="disk")
- nbytes: the memory used by the class
- see :class:`~beast.physicsmodel.grid.SEDGrid` and :class:`~beast.physicsmodel.grid.SpectralGrid`
  for more attributes

Both grids can interface with files via backends on disk with differing
levels memory usage.

- memory:
  The :class:`~beast.physicsmodel.helpers.gridbackends.MemoryBackend` reads
  contents of the file into memory.  Largest memory usage and
  fastest access to the information after reading is done.

- cache:
  The :class:`~beast.physicsmodel.helpers.gridbackends.CacheBackend` reads the
  contents into memory when that specific information is first accessed.

- disk:
  The :class:`~beast.physicsmodel.helpers.gridbackends.DiskBackend` reads the
  data from disk when it is accessed.  This supports reading only a portion
  of the data.  This allows spectral and SED grids larger than can fit into
  memory.

The current file formats that are supported by the backends are:

- fits: Supported by memory and cache backends (disk backend support TBD)

- hdf: Supported by memory, cache, and disk backends
