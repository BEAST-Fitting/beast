#########################
Installation instructions
#########################

Requirements
============

Running the BEAST requires:

- Python >= 3.7
- Astropy >= 1.3

In turn, Astropy depends on
`other packages <http://docs.astropy.org/en/latest/install.html>`_ for
optional features. From these you will need:

- `hdf5 <http://h5py.org/>`_ to read/write `Table` objects from/to HDF5 files.

You will also need:

- `PyTables <http://www.pytables.org/>`_ to manage large amounts of data.

One easy way to obtain the above is through the AstroConda Python stack:

- First install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ which
  contains the conda package manager. Once Miniconda is installed,
  you can use the `conda` command to install any other packages and create
  environments, etc.

- Setup the AstroConda Channel

.. code-block:: console

    $ conda config --add channels http://ssb.stsci.edu/astroconda

- Install AstroConda with Python 3 (recommended)

.. code-block:: console

    $ conda create -n astroconda stsci

- Ensure the minimum Python requirements above are met. If needed, specify the minimum version required

.. code-block:: console

    $ python --version
    $ conda create -n astroconda stsci python=3.7

- Make sure that the ``PyTables`` and ``hdf5`` packages are installed

.. code-block:: console

    $ conda install -n astroconda pytables

    $ conda install -n astroconda hdf5


Installation for Users
======================

In addition to installing the code, library files also need to be installed.
See :ref:`library-files`.

.. note::
   The ``numpy`` package must be installed before installing the BEAST. This package is needed as part of the installation process, not just as a prerequisite for running the BEAST.

Using pip
---------

``beast`` can be installed using pip:

.. code-block:: console

    $ pip install beast

If you already have an older version installed, use:

.. code-block:: console

    $ pip install --upgrade beast

Or, from the master trunk on the repository, considered developmental code:

.. code-block:: console

    $ pip install git+https://github.com/BEAST-Fitting/beast.git

From source
-----------

If you are happy with your current environment, ``beast`` can also be installed from
the source code in the normal python fashion after cloning from the git repo or
downloading from Github

.. code-block:: console

     $ pip install .

If you are using conda, you may wish to create a conda environment with the
dependencies before doing the install

.. code-block:: console

     $ conda env create -n beast --file conda-environment.yml
     $ conda activate beast
     $ pip install .


.. _Dev_install:

Installation for Developers
===========================

If you plan on modifying the ``beast`` in addition to running the code, it may
be useful to create a development installation. First, create a fork of the
official ``beast`` repository and clone it:

.. code-block:: console

   $ git clone https://github.com/YourName/beast.git

Optionally, you can rename this cloned copy:

.. code-block:: console

   $ git clone https://github.com/YourName/beast.git beast-YourName

Set the value of the fork's 'upstream' to the official distribution so you
can incorporate changes made by others to your development fork. In the clone
of your fork, run the following:

.. code-block:: console

   $ git remote add upstream https://github.com/BEAST-Fitting/beast.git

In order to run a development installation, navigate to the directory in your
``beast`` repository that contains `setup.py`, and run:

.. code-block:: console

   $ pip install -e .


Installation for Production Runs
================================

Using a dedicated conda environment for production BEAST runs may be
desirable. Such an environment provides a way to ensure that
production runs are reproducible by fixing the versions of all the
software used. The instructions below assume that the `astroconda channel
<https://astroconda.readthedocs.io/>`_ is being used.

Create a conda environment, and name it differently from your default and other environments:

  .. code-block:: console

    $ conda create -n beast_production python=3.6

Activate the environment after all the packages are finished installing:

  .. code-block:: console

    $ conda activate beast_production

Install dependencies using conda (better for speed):

  .. code-block:: console

    $ conda install astropy scipy h5py matplotlib cython

Next, install the BEAST. You have three options:

Option 1: Use pip to install a stable version of the BEAST (e.g. v1.3):

  .. code-block:: console

    $ pip install beast==1.3

Option 2: Get the latest production branch, which can be ahead of the pypi version:

  .. code-block:: console

    $ pip install git+https://github.com/BEAST-Fitting/beast.git@v1.x

Option 3: If you'll be doing development, follow the instructions above (:ref:`Dev_install`). Note that you can make separate environments for development and production modes.

The BEAST production version is now ready for use. Note, you need to
activate this conda environment every time you want to use this installed
version.

.. _library-files:

BEAST Library Files
===================

For the BEAST to work properly, you need to place a set of files in a
directory. These files contain information related to filters,
stellar atmospheres, and in the future stellar evolution models.

Manual download
---------------

The required library files can be manually acquired from:

https://stsci.box.com/v/beastlibs

Note that the archive at this link contains a folder called `files`. The
*contents* of this folder are the library files required by the BEAST. These files need to be placed within (any of) the possible locations specified in :ref:`library_loc`.

Script download
---------------

Alternatively, after installing the BEAST, run the following script and the library files will be downloaded into the location specified in :ref:`library_loc`:

.. code-block:: console

     $ beast get_libfiles

.. _library_loc:

Location
--------

There are 2 possible locations for the required library files. The possible locations are
(in the order the code will search for them):

1. In a directory designated by a BEAST_LIBS environment variable.
2. In the ``.beast`` directory in the home directory of the current user (ie, ``~/.beast``);
   this is usually the easiest and will be automatically created if it doesn't exist.
