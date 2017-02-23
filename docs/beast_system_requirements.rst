System Requirements
===================

Running the BEAST requires a Python 2.7 installation along with some
standard astronomy packages. One easy way to obtain the pre-requisites
is through the AstroConda Python stack.

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
