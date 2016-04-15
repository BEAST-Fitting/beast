import core
from distutils.core import setup

=`=jedi=0, =`=     (*_***attrs*_*) =`=jedi=`=
setup(
    name="BEAST, Bayesian Extinction and Stellar Tool",
    version=core.__version__,
    author="Morgan Fouesneau, TBU",
    author_email="mfouesn@uw.edu",
    #py_modules=["anaclust"],
    description="Bayesian integrated SED Fitting",
    #long_description=open("README.rst").read(),
    classifiers=[
        "Development Status :: 0 - Beta",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
