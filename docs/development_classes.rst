*******
Classes
*******

The BEAST uses classes to abstract the specific implementation with the interface.

Grids
=====

The two main classes to use to access information in spectral and SED grids are:

- :class:`~beast.physicsmodel.grid.SpectralGrid`:
  Provides the interface to grids of spectra.
  Generally used as part of generating the physics grid.

- :class:`~beast.physicsmodel.grid.SEDGrid`:
  Provides the interface to grids of extinguished SEDs composed of band integrated fluxes.
  Used fairly extensively in the BEAST to interface with the core physicsmodel grid.

Both grids can interface with files via "backends" on disk with differing
levels memory usage.

- "memory":
  This backend reads contents of the file into memory.  Largest memory usage and
  fastest access to the information after reading is done.

- "cache":
  This backend reads the contents into memory when that specific information is
  first accessed.

- "lazycache": TBD
  This backend is TBD.

The current file formats that are supported by the backends are:

- "fits": Supported by "memory" and "cache" backends.

- "hdf5": Supported by "memory" and "cache" backends

- "lazycache": TBD.
