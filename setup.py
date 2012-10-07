import anased
from distutils.core import setup


setup(
    name="anased",
    version=anased.__version__,
    author="Morgan Fouesneau",
    author_email="mfouesn@uw.edu",
    #py_modules=["anaclust"],
    description="Bayesian integrated SED Fitting",
    long_description=open("README").read(),
    classifiers=[
        "Development Status :: 0 - Beta",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
