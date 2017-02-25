BEAST First Run
===============


AstroConda setup
================

Running the BEAST requires a Python 2.7 installation along with some
standard astronomy packages.

- Install Miniconda. This is a lightweight 'conda' install
  that does not acually install Python or any packages.
  Installers: < https://conda.io/miniconda.html >

- Setup the AstroConda Channel
  - $ conda config --add channels http://ssb.stsci.edu/astroconda

- Install AstroCOnda with Python 2.7
  - $ conda create -n astroconda2 python=2.7 stsci pyraf iraf

- Make sure that the 'pytables' and 'hdf5' packages are installed
  - $ sudo conda install -n astroconda2 pytables
  - $ sudo conda install -n astroconda2 hdf5


Obtaining the BEAST
===================

Fetch the BEAST distro from GitHub

- Checkout the BEAST from GitHub. This will create a directory
  named 'beast' and a sub-directory 'beast/beast'
  - $ git clone https://github.com/karllark/beast.git

- Obtain the BEAST 'libs' file and copy it into 'beast/beast',
  or place a soft link pointing to it, e.g.,
  - $ cd beast/beast
  - $ wget -r <www.TBD.edu/libs>
  - OR
  - ln -s /usr/local/bin/beast-libs libs



Hello World with the BEAST
==========================
    
- In 'beast/beast/examples/phat_small', place a soft link
  named 'beast' pointing two levels up
  - $ cd beast/beast/examples/phat_small
  - $ ln -s ../../ beast

- Activate AstroConda setup that contains Python 2.7
  - $ source activate astroconda2

- Check out the basic gelp message of 'run_beast.py'
  - $ ./run_beast.py -h

- Make BEAST say 'Hello World'!
  - $ ./run_beast.py -potf
