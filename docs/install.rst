Installation
============

Requirements
------------

Running the BEAST requires:

- Python 2.7
- Astropy 2.0

In turn, Astopy depends on 
`other packages <http://docs.astropy.org/en/latest/install.html>`_ for 
optional features. From these you will need:

- 'hdf5' to read/write Table objects from/to HDF5 files.
- other astropy dependencies (list these here). 

You will also need:

- `pytables <http://www.pytables.org/>`_ to manage large amounts of data.

One easy way to obtain the above is through the AstroConda Python stack:

- First, install `Miniconda < https://conda.io/miniconda.html >`_ which 
  contains the conda package manager and Python. Once Miniconda is installed,
  you can use the conda command to install any other packages and create 
  environments, etc.

- Setup the AstroConda Channel
  - ``$ conda config --add channels http://ssb.stsci.edu/astroconda``

- Install AstroConda with Python 2.7
  - ``$ conda create -n astroconda2 python=2.7 stsci pyraf iraf``

- Make sure that the 'pytables' and 'hdf5' packages are installed
  - ``$ sudo conda install -n astroconda2 pytables``
  - ``$ sudo conda install -n astroconda2 hdf5``


Installing the BEAST
--------------------

You will need a `GitHub <https://github.com/>`_ account to install the BEAST and
you will also need to 
`set up Git on your computer <https://help.github.com/articles/set-up-git/>`_.
Here is a general introduction to `GitHub <https://help.github.com/>`_.

There are three installation options listed below. In addition to 
choosing one of the three options, you will also need to download a BEAST 'libs'
files. See `Obtaining BEAST libs <http://beast.readthedocs.io/en/latest/install.html#obtaining-beast-libs>`_ below.

Option 1 
--------
This option installs the package on your local machine as a clone of the
`BEAST GitHub repository <https://github.com/karllark/beast>`_
such that you can easily keep up with code updates. If you
*do not plan to make code contributions* this is the recommended installation 
method. See Option 3 below if you if you would like to contribute 
to code enhancements.


Type the following on the command line in your local directory where you want
to place the BEAST: ``$ git clone https://github.com/karllark/beast.git``
This will create a directory named 'beast' containing the package contents.
Type `ls` to confirm the BEAST has been downloaded.

Option 2
--------

The second option installs the package in its current version on your local 
computer. To obtain package updates you will have to again manually install the 
package:

Go to the `BEAST GitHub homepage <https://github.com/karllark/beast>`_ and 
click on the green 'Clone or download' button, then select 'Download ZIP'. 
Unzip the file and place the resulting 'beast-master' folder in the directory 
from which you would like to work on it.
   
Option 3
--------

The third option is suitable if you *plan to make code contributions* to the
BEAST.
   
Go to the `BEAST GitHub homepage <https://github.com/karllark/beast>`_.
Create your own copy of the BEAST repository by creating a BEAST fork by
clicking on the fork button in the upper right of the repository page. Creating 
your own copy will allow you to use and edit the code according to your needs.
Make a clone of your own BEAST repository. To do this go to
``github.com/YOURGITNAME/beast`` and click on the green `Clone or Download`
button. When the small window opens, copy the HTTPS 
(``https://github.com/YOURGITNAME/beast.git``) or SSH code 
(``git@github.com:YOURGITNAME/beast.git``). On the command line in a local 
directory where you want to place the BEAST, type `git clone` and then paste 
the copied code. This should look like
``$ git clone https://github.com/YOURGITNAME/beast.git`` or
``$ git clone git@github.com:YOURGITNAME/beast.git``. 
   
To confirm that your BEAST repository was placed in the local directory, type 
``ls`` - you should see a new folder called `BEAST` created in this directory.

For additional instructions on how you can use the BEAST and eventually make
code contributions, see the 
`BEAST Development Documentation <http://beast.readthedocs.io/en/latest/beast_development.rst>`_.


Obtaining BEAST Library Files
-----------------------------

For the BEAST to work properly, you need to place these files in the 
'beast/beast' directory. If you have the files copy them to this location. 
You should have a resulting 'beast/beast/libs' directory. Otherwise you can 
place a soft link to the files (Rubab - confirm your instructions):

``$ cd beast/beast`` (or '$ cd beast-master/beast' if you used install option 2)
``$ wget -r <www.TBD.edu/libs>``
`Obtain wget for Windows <http://gnuwin32.sourceforge.net/packages/wget.htm>`_.
You can also use:
``$ cd beast/beast``
``$ ln -s /usr/local/bin/beast-libs libs``


Confirm Proper Installation and Run Sample Code
-----------------------------------------------

There is a small sample script named 'run_beast.py' located in
'beast/beast/examples/phat_small' as a quick check to confirm the BEAST 
installation is working. 

In 'beast/beast/examples/phat_small', place a soft link named 'beast' 
pointing two levels up:  
``$ cd beast/beast/examples/phat_small``
``$ ln -s ../../ beast``

Verify that the default Python installation is version 2.7:
``$ python --version``

If you installed Python through AstroConda, first activate the correct 
AstroConda environment:
``$ source activate astroconda2``

Take a look at the basic help message of 'run_beast.py'
``$ ./run_beast.py -h``

Now try out a sample BEAST run.
``$ ./run_beast.py -potf``


If BEAST is running correctly the code should run without errors and the 
output should be a plot which looks like this (insert plot here):
