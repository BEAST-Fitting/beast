######
Priors
######

Priors on the parameters are set in the datamodel.py file using
python dictionaries to give the prior model and any related
parameters.

New prior models can be added in the specific code regions by defining
a new prior model name and creating the necessary code to support it.
The code files of interest are
beast/physicsmodel/prior_weights_[stars,dust].py.
New models added should be documented here.

Stellar
=======

Age
---

Mass
----

The mass prior is set by the choice of an Initial Mass Function (IMF).
The two mass function supported are:

Kroupa [default option] (details needed)

.. code-block:: python

  mass_prior_model = {'name': 'kroupa'}

Salpeter (details needed)

.. code-block:: python

  mass_prior_model = {'name': 'salpeter'}

Metallicity
-----------

Extinction
==========

A(V)
----

R(V)
----

f_A
---

Distance
========
