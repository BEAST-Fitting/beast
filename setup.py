import core
from distutils.core import setup


setup(
    name="SED STELLAR NUTS CRACKER",
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
