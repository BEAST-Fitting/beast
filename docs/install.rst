############
Installation
############

Requirements
============

Running the BEAST requires:

- Python >=3.5
- Astropy >=1.3

In turn, Astropy depends on
`other packages <http://docs.astropy.org/en/latest/install.html>`_ for
optional features. From these you will need:

- `hdf5 <http://h5py.org/>`_ to read/write `Table` objects from/to HDF5 files.

You will also need:

- `PyTables <http://www.pytables.org/>`_ to manage large amounts of data.

One easy way to obtain the above is through the AstroConda Python stack:

- First install `Miniconda <https://conda.io/miniconda.html>`_ which
  contains the conda package manager. Once Miniconda is installed,
  you can use the `conda` command to install any other packages and create
  environments, etc.

- Setup the AstroConda Channel

.. code-block:: console

    $ conda config --add channels http://ssb.stsci.edu/astroconda

- Install AstroConda with Python 3 (recommended)

.. code-block:: console

    $ conda create -n astroconda stsci

- Make sure that the ``PyTables`` and ``hdf5`` packages are installed

.. code-block:: console

    $ conda install -n astroconda pytables

    $ conda install -n astroconda hdf5


Installing the BEAST
====================

In addition to installing the code, library files also need to be installed.
See :ref:`library-files`.

Using pip
---------

``beast`` can also be installed using pip

.. code-block:: console

    $ pip install beast

if you already have an older version installed

.. code-block:: console

    $ pip install --upgrade beast

from the master trunk on the repository, considered developmental code

.. code-block:: console

    $ pip install git+https://github.com/BEAST-Fitting/beast.git

From source
-----------

If you are happy with your current environment, ``beast`` can be installed from
the source code in the normal python fashion after cloning from the git repo or
downloading from Github

.. code-block:: console

     $ python setup.py install

If you are using conda, you may wish to create a conda environment with the
dependencies before doing the install

.. code-block:: console

     $ conda env create -n beast --file conda-environment.yml
     $ conda activate beast
     $ python setup.py install

If you would like to modify beast, you may want to use links instead of
installing, which is best done by replacing the last line with

.. code-block:: console

     $ python setup.py develop


For developers
--------------

This option is suitable if you *plan to make code contributions* to the BEAST.
See the :ref:`beast_development` for details.

.. _library-files:

BEAST Library Files
===================

For the BEAST to work properly, you need to place a set of files in a
directory.  These files contain information related to filters,
stellar atmospheres, and in the future stellar evolution models.

.. _library_loc:

Location
--------

There are 3 possible locations for these files (in the order the code
will search for them):

1. in a directory designated by the BEAST_LIBS environment variable
2. in the '.beast' directory in the home directory of the current user
3. in the source code in 'beast/beast/libs'

Whichever of the options used, the directory needs to be manually created.

Manual download
---------------

<https://stsci.box.com/v/beastlibs>

Note that the archive at this link contains a folder called `files`. The
_contents_ of this folder are the library files required by beast. It is these
files that need to be places within (any of) the possible locations given above.

Script download
---------------

After installing the BEAST, run the following script and the library files
will be downloaded into the location specified in :ref:`library_loc`

.. code-block:: console

     $ python -m beast.tools.get_libfiles

Running Example
===============

You can find examples of BEAST runs in the ``beast/examples/`` directory.

Inside each example, there is a run_beast*.py script.

phat_small example
------------------

This example is based on a *very* small amount of old PHAT data.

If the beast has not been installed (only downloaded from github), then
In the 'phat_small' directory, place a soft link named 'beast' to where the
beast code is located.  Specifically

.. code-block:: console

    $ cd examples/phat_small

    $ ln -s beast_code_loc/beast/beast beast

If you installed Python through AstroConda, first activate the correct
AstroConda environment

.. code-block:: console

    $ source activate astroconda

Verify that the current default Python is version 3

.. code-block:: console

    $ python --version

Next, bring up the BEAST help message, which describes the available switch
options, with

.. code-block:: console

    $ ./run_beast.py -h

You should be presented with the following options::

  -h, --help              show this help message and exit
  -p, --physicsmodel      Generate the model grid
  -o, --observationmodel  Calculate the noise model
  -t, --trim              Trim the model and noise grids
  -f, --fit               Fit the observed data
  -r, --resume            Resume a run

Now launch a sample BEAST run (with flags set to run through the full
sequence of generation of physics model, observation model generation, trimming
of the grid, and fitting to the observed data) using

.. code-block:: console

  $ ./run_beast.py -potf

If the BEAST is running correctly, this command should run without errors
and should have written the output files into 'beast_example_phat/'. The result
can be plotted using

.. code-block:: console

    $ python beast/plotting/plot_indiv_fit.py beast_example_phat/beast_example_phat

The argument for this script is the prefix of the output files. The output
should look like this

.. image:: beast_example_phat_ifit_starnum_0.png
