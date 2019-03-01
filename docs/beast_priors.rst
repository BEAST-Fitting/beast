.. _beast_priors:

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

The age prior is the star formation rate (SFR) and can be

1. Flat in linear age [default option]

.. code-block:: python

  age_prior_model = {'name': 'flat'}

2. Set by bins spaced in logage.  For example:

.. code-block:: python

  age_prior_model = {'name': 'bins',
                     'logages': [6.0, 7.0, 8.0, 9.0, 10.0],
                     'values': [1.0, 2.0, 1.0, 5.0, 3.0]}

3. An exponentially decreasing SFR starting at 1.0,
with a 1000 Myr time constant is:

.. code-block:: python

  age_prior_model = {'name': 'exp',
                     'A': 1.0,
                     'tau': 1000.}

Mass
----

The mass prior is set by the choice of an Initial Mass Function (IMF).
The two mass function supported are:

1. Kroupa [default option] (details needed)

.. code-block:: python

  mass_prior_model = {'name': 'kroupa'}

2. Salpeter (details needed)

.. code-block:: python

  mass_prior_model = {'name': 'salpeter'}

Metallicity
-----------

The metallicity prior can be

1. Flat [default option]

.. code-block:: python

  met_prior_model = {'name': 'flat'}

Extinction
==========

A(V)
----

The A(V) prior can be:

1. Flat [default option]

.. code-block:: python

  av_prior_model = {'name': 'flat'}

2. Lognormal with the maximum at the A(V) given by max_pos, the width
given by sigma, and the number at max given by N.

.. code-block:: python

  av_prior_model = {'name': 'lognormal',
                    'max_pos': 2.0,
                    'sigma': 1.0,
                    'N': 10.}

3. Two lognormals (see above for definition of terms)

.. code-block:: python

  av_prior_model = {'name': 'two_lognormal',
                    'max_pos1': 0.2,
                    'max_pos1': 2.0,
                    'sigma1': 0.5,
                    'sigma2': 2.0,
                    'N1': 10.,
                    'N2': 20.}

4. Exponential with decay rate 'a' and amplitude 'N'

.. code-block:: python

  av_prior_model = {'name': 'exponential',
                    'a': 2.0,
                    'N': 10.}

R(V)
----

1. Flat [default option]

.. code-block:: python

  rv_prior_model = {'name': 'flat'}

2. Lognormal with the maximum at the R(V) given by max_pos, the width
given by sigma, and the number at max given by N.

.. code-block:: python

  rv_prior_model = {'name': 'lognormal',
                    'max_pos': 2.0,
                    'sigma': 1.0,
                    'N': 10.}

3. Two lognormals (see above for definition of terms)

.. code-block:: python

  rv_prior_model = {'name': 'two_lognormal',
                    'max_pos1': 0.2,
                    'max_pos1': 2.0,
                    'sigma1': 0.5,
                    'sigma2': 2.0,
                    'N1': 10.,
                    'N2': 20.}

f_A
---

1. Flat [default option]

.. code-block:: python

  fA_prior_model = {'name': 'flat'}

2. Lognormal with the maximum at the f_A given by max_pos, the width
given by sigma, and the number at max given by N.

.. code-block:: python

  fA_prior_model = {'name': 'lognormal',
                    'max_pos': 2.0,
                    'sigma': 1.0,
                    'N': 10.}

3. Two lognormals (see above for definition of terms)

.. code-block:: python

  fA_prior_model = {'name': 'two_lognormal',
                    'max_pos1': 0.2,
                    'max_pos1': 2.0,
                    'sigma1': 0.5,
                    'sigma2': 2.0,
                    'N1': 10.,
                    'N2': 20.}

Distance
========

[TBD]
