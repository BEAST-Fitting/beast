#######
Filters
#######

Filters included in the BEAST

Hubble
======

WFC3
----

Imaging, UVIS, Wide

.. plot::

   from beast.plotting.plot_filters import plot_filters

   wfc3_uvis = ["f218w","f225w", "f275w", "f336w", "f390w", "f438w",
                "f475w", "f555w", "f606w", "f625w", "f775w", "f814w"]
   filters = [f"HST_WFC3_{cfilt.upper()}" for cfilt in wfc3_uvis]
   plot_filters(filters)

Imaging, UVIS, Medium

.. plot::

   from beast.plotting.plot_filters import plot_filters

   wfc3_uvis = ["f390m", "f410m", "f467m", "f547m", "f621m", "f689m", "f763m", "f845m"]
   filters = [f"HST_WFC3_{cfilt.upper()}" for cfilt in wfc3_uvis]
   plot_filters(filters)

Imaging, UVIS, Narrow, Blue

.. plot::

   from beast.plotting.plot_filters import plot_filters

   wfc3_uvis = ["f280n", "f343n", "f373n", "f395n", "f469n", "f487n", "f502n", "f631n"]
   filters = [f"HST_WFC3_{cfilt.upper()}" for cfilt in wfc3_uvis]
   plot_filters(filters)

Imaging, UVIS, Narrow, Red 
(Note: F656N not included, error when plotting)

.. plot::

   from beast.plotting.plot_filters import plot_filters

   wfc3_uvis = ["f645n", "f657n", "f658n", "f665n", "f673n", "f680n", "f953n"]
   filters = [f"HST_WFC3_{cfilt.upper()}" for cfilt in wfc3_uvis]
   plot_filters(filters)

Imaging, IR, wide

.. plot::

   from beast.plotting.plot_filters import plot_filters

   wfc3_ir = ["f105w", "f110w", "f125w", "f140w", "f160w"]
   filters = [f"HST_WFC3_{cfilt.upper()}" for cfilt in wfc3_ir]
   plot_filters(filters)

Imaging, IR, medium

.. plot::

   from beast.plotting.plot_filters import plot_filters

   wfc3_ir = ["f098m", "f127m", "f139m", "f153m"]
   filters = [f"HST_WFC3_{cfilt.upper()}" for cfilt in wfc3_ir]
   plot_filters(filters)

ACS
---

Imaging, WFC, wide

.. plot::

   from beast.plotting.plot_filters import plot_filters

   acs_wfc = ["f435w", "f475w", "f555w", "f606w", "f625w", "f775w", "f814w"]
   filters = [f"HST_ACS_WFC_{cfilt.upper()}" for cfilt in acs_wfc]
   plot_filters(filters)


Imaging, WFC, extrawide, medium and narrow

.. plot::

   from beast.plotting.plot_filters import plot_filters

   acs_wfc = ["f850lp", "f502n", "f550m", "f658n"]
   filters = [f"HST_ACS_WFC_{cfilt.upper()}" for cfilt in acs_wfc]
   plot_filters(filters)

Imaging, SBC

.. plot::

   from beast.plotting.plot_filters import plot_filters

   acs_sbc = ["f115lp", "f125lp", "f140lp", "f150lp", "f165lp", "f122m"]
   filters = [f"HST_ACS_SBC_{cfilt.upper()}" for cfilt in acs_sbc]
   plot_filters(filters)

WFPC2
-----

Imaging, UV, Wide

.. plot::

   from beast.plotting.plot_filters import plot_filters

   wfpc2 = ["f157w", "f170w", "f185w", "f218w", "f255w", "f300w"]
   filters = [f"HST_WFPC2_{cfilt.upper()}" for cfilt in wfpc2]
   plot_filters(filters)

Imaging, Optical, Blue

.. plot::

   from beast.plotting.plot_filters import plot_filters

   wfpc2 = ["f336w", "f380w", "f439w", "f450w", "f555w", "f569w"]
   filters = [f"HST_WFPC2_{cfilt.upper()}" for cfilt in wfpc2]
   plot_filters(filters)

Imaging, Optical, Red

.. plot::

   from beast.plotting.plot_filters import plot_filters

   wfpc2 = ["f606w", "f622w", "f675w", "f702w", "f791w", "f814w"]
   filters = [f"HST_WFPC2_{cfilt.upper()}" for cfilt in wfpc2]
   plot_filters(filters)

Imaging, Medium

.. plot::

   from beast.plotting.plot_filters import plot_filters

   wfpc2 = ["f122m", "f410m", "f467m", "f547m"]
   filters = [f"HST_WFPC2_{cfilt.upper()}" for cfilt in wfpc2]
   plot_filters(filters)

Webb
====

NIRCam
------

Imaging, Very Wide Bands

.. plot::

   from beast.plotting.plot_filters import plot_filters

   jwst_nircam = ["f150w2", "f332w2"]
   filters = [f"JWST_NIRCAM_{cfilt.upper()}" for cfilt in jwst_nircam]
   plot_filters(filters)

Imaging, Wide Bands

.. plot::

   from beast.plotting.plot_filters import plot_filters

   jwst_nircam = ["f070w", "f090w", "f115w", "f150w", "f200w",
                  "f277w", "f356w", "f444w"]
   filters = [f"JWST_NIRCAM_{cfilt.upper()}" for cfilt in jwst_nircam]
   plot_filters(filters)

Imaging, Medium Bands

.. plot::

   from beast.plotting.plot_filters import plot_filters

   jwst_nircam = ["f140m", "f162m", "f182m", "f210m",
                  "f250m", "f300m", "f335m", "f360m", "f410m", "f430m", "f460m", "f480m"]
   filters = [f"JWST_NIRCAM_{cfilt.upper()}" for cfilt in jwst_nircam]
   plot_filters(filters)

Imaging, Narrow Bands

.. plot::

   from beast.plotting.plot_filters import plot_filters

   jwst_nircam = ["f164n", "f187n", "f212n",
                  "f323n", "f405n", "f466n", "f470n"]
   filters = [f"JWST_NIRCAM_{cfilt.upper()}" for cfilt in jwst_nircam]
   plot_filters(filters)

NIRISS
------

Imaging, Wide Bands

.. plot::

   from beast.plotting.plot_filters import plot_filters

   jwst_niriss = ["f090w", "f115w", "f150w", "f200w", "f277w", "f356w", "f444w"]
   filters = [f"JWST_NIRISS_{cfilt.upper()}" for cfilt in jwst_niriss]
   plot_filters(filters)

Imaging, Medium Bands

.. plot::

   from beast.plotting.plot_filters import plot_filters

   jwst_niriss = ["f140m", "f158m", "f380m", "f430m", "f480m"]
   filters = [f"JWST_NIRISS_{cfilt.upper()}" for cfilt in jwst_niriss]
   plot_filters(filters)

MIRI
----

Imaging

.. plot::

   from beast.plotting.plot_filters import plot_filters

   jwst_miri = ["f560w", "f770w", "f1000w", "f1130w", "f1280w", "f1500w", "f1800w", "f2100w", "f2550w"]
   filters = [f"JWST_MIRI_{cfilt.upper()}" for cfilt in jwst_miri]
   plot_filters(filters)

Coronagraphy

.. plot::

   from beast.plotting.plot_filters import plot_filters

   jwst_miri = ["f1065c", "f1140c", "f1550c", "f2300c"]
   filters = [f"JWST_MIRI_{cfilt.upper()}" for cfilt in jwst_miri]
   plot_filters(filters)

GALEX
=====

Imaging

.. plot::

   from beast.plotting.plot_filters import plot_filters

   plot_filters(["GALEX_FUV", "GALEX_NUV"])

