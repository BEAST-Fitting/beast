EZMIST -- A python package that allows you to download MIST/MESA isochrones directly from their website
=======================================================================================================


This small package provides a direct interface to the MIST/MESA isochrone
webpage (http://waps.cfa.harvard.edu/MIST).
It compiles the URL needed to query the website and retrives the data into a
python variable.

This package has been tested on python 2.7 and python 3.

:version: 1
:author: MF

(this package is similar to EzPadova:  https://github.com/mfouesneau/ezpadova)


EXAMPLE USAGE
-------------

* Basic example of downloading a sequence of isochrones, plotting, saving
```python 
>>> r = ezmist.get_t_isochrones(6.0, 7.0, 0.05, FeH_value=0.0, theory_output='full')
>>> import pylab as plt
>>> plt.scatter(r['logT'], r['logL'], c=r['logA'], edgecolor='None')
>>> plt.show()
>>> r.write('myiso.fits')
```

Note: MIST isochrone metallicities are defined in terms of [Fe/H] (not Z) 

* getting only one isochrone
```python 
>>> r = ezmist.get_one_isochrones(1e7, 0.0, age_scale='linear')
```
