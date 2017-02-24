Installation
============

Requirements
------------

Running the BEAST requires:

- Python 2.7
- Astropy 1.3

In turn, Astropy depends on 
`other packages <http://docs.astropy.org/en/latest/install.html>`_ for 
optional features. From these you will need:

- `hdf5 <http://h5py.org/>`_ to read/write `Table` objects from/to HDF5 files.
- other Astropy dependencies [list these here].

You will also need:

- `PyTables <http://www.pytables.org/>`_ to manage large amounts of data.

One easy way to obtain the above is through the AstroConda Python stack:

- First install `Miniconda <https://conda.io/miniconda.html>`_ which 
  contains the conda package manager. Once Miniconda is installed,
  you can use the `conda` command to install any other packages and create 
  environments, etc.

- Setup the AstroConda Channel:

``$ conda config --add channels http://ssb.stsci.edu/astroconda``

- Install AstroConda with Python 2.7:

``$ conda create -n astroconda2 python=2.7 stsci pyraf iraf``

- Make sure that the ``PyTables`` and ``hdf5`` packages are installed:

``$ sudo conda install -n astroconda2 pytables``

``$ sudo conda install -n astroconda2 hdf5``


Installing the BEAST
--------------------

For the BEAST to work properly, you need to place a 'libs' directory containing
files related to filters, photometry, stellar atmospheres, and in the future,
stellar evolution models. See **Obtaining BEAST libraries** below. 
[Place an internal crossreference here].

Option 1 
________

The following is the recommended option which will allow you to easily keep up 
with code updates. If additionally you would like to contribute to code 
enhancements, see Option 3 below.

Places BEAST on your local computer as a clone of the
`BEAST GitHub repository <https://github.com/karllark/beast>`_. To do this, go 
to the directory where you want to place the BEAST and type the following:

``$ git clone https://github.com/karllark/beast.git``

This will create a directory named 'beast' containing the BEAST.

Option 2
________

This option places BEAST in its current version on your local computer. One 
disadvantage is that you will have to manually obtain BEAST updates from the 
GitHub repository if you use this option.

On the `BEAST GitHub homepage <https://github.com/karllark/beast>`_ click on 
the 'Clone or Download' button, then select 'Download ZIP'. Unzip the 
file in the desired directory.
   
Option 3
________

This option is suitable if you *plan to make code contributions* to the BEAST.
See the `BEAST Development Documentation <http://beast.readthedocs.io/en/latest/beast_development.rst>`_
for details. [Update link]


.. _`Obtaining BEAST libraries`
Obtaining BEAST Library Files
-----------------------------

For the BEAST to work properly, you need to place a 'libs' directory containing
files related to filters, photometry, stellar atmospheres, and in the future
stellar evolution models in the 'beast/beast' directory. Follow this link to 
download the 'libs' files. [Update link]


Confirming Proper Installation and Running Sample Code
------------------------------------------------------

There is a small sample script named 'run_beast.py' located in
'beast/beast/examples/phat_small' as a quick check to confirm that BEAST 
installation is working.

In 'beast/beast/examples/phat_small', place a soft link named 'beast' 
pointing two levels up:  

``$ cd beast/beast/examples/phat_small``

``$ ln -s ../../ beast``

If you installed Python through AstroConda, first activate the correct 
AstroConda environment:

``$ source activate astroconda2``

Verify that the current default Python is version 2.7:

``$ python --version``

Now try a sample BEAST run:

``$ ./run_beast.py`` or ``$ python run_beast.py``

Optionally, you can run BEAST with one, or a combination, of these arguments

  -h, --help    show this help message and exit
  -m, --models  Generate the model grid
  -n, --noise   Calculate the noise model
  -t, --trim    Trim the model and noise grids
  -f, --fit     Fit the observed data
  -r, --resume  Resume a run

``$ ./run_beast.py -h`` or ``$ ./run_beast.py -mntf``

[Optional]: If the BEAST is running correctly there should be no errors and the 
output should be a plot which looks like this [insert plot here?]:
