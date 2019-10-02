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

1. Flat in linear age

.. code-block:: python

  age_prior_model = {'name': 'flat'}

or

.. code-block:: python

  age_prior_model = {'name': 'flat_linear'}

2. Flat in log age

.. code-block:: python

  age_prior_model = {'name': 'flat_log'}

3. Set by bins spaced in logage.

For example, step like priors can be specified by:

.. code-block:: python

  age_prior_model = {'name': 'bins_histo',
                     'logages': [6.0, 7.0, 8.0, 9.0, 10.0],
                     'values': [1.0, 2.0, 1.0, 5.0, 3.0]}

For example, lines connecting the bin value of the priors can be specified by:

.. code-block:: python

 age_prior_model = {'name': 'bins_interp',
                    'logages': [6.0, 7.0, 8.0, 9.0, 10.0],
                    'values': [1.0, 2.0, 1.0, 5.0, 3.0]}

4. An exponentially decreasing SFR starting at 1.0,
with a 1000 Myr time constant is:

.. code-block:: python

  age_prior_model = {'name': 'exp',
                     'A': 1.0,
                     'tau': 1000.}

Plot showing examples of the possible age prior models with the parameters given above.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from beast.physicsmodel.prior_weights_stars import compute_age_prior_weights

    fig, ax = plt.subplots()

    # logage grid from 1 Myrs to 10 Gyrs
    logages = np.linspace(6.0, 10.0)

    age_prior_models = [
        {"name": "flat"},
        {"name": "flat_log"},
        {
            "name": "bins_histo",
            "logages": [6.0, 7.0, 8.0, 9.0, 10.0],
            "values": [1.0, 2.0, 1.0, 5.0, 3.0],
        },
        {
            "name": "bins_interp",
            "logages": [6.0, 7.0, 8.0, 9.0, 10.0],
            "values": [1.0, 2.0, 1.0, 5.0, 3.0],
        },
        {"name": "exp", "A": 1.0, "tau": 100.}
    ]

    for ap_mod in age_prior_models:
        ax.plot(logages, compute_age_prior_weights(logages, ap_mod), label=ap_mod["name"])

    ax.set_ylabel("probability")
    ax.set_xlabel("log(age)")
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()


Mass
----

The mass prior is set by the choice of an Initial Mass Function (IMF).
The two mass function supported are:

1. Kroupa (details needed)

.. code-block:: python

  mass_prior_model = {'name': 'kroupa'}

2. Salpeter (details needed)

.. code-block:: python

  mass_prior_model = {'name': 'salpeter'}

Plot showing examples of the possible mass prior models with the parameters given above.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from beast.physicsmodel.prior_weights_stars import compute_mass_prior_weights

    fig, ax = plt.subplots()

    # mass grid from 0.01 to 100 solar masses (log spacing)
    masses = np.logspace(-2.0, 2.0)

    mass_prior_models = [
        {"name": "kroupa"},
        {"name": "salpeter"}
    ]

    for mp_mod in mass_prior_models:
        ax.plot(masses, compute_mass_prior_weights(masses, mp_mod), label=mp_mod["name"])

    ax.set_ylabel("probability")
    ax.set_xlabel("mass")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()

Metallicity
-----------

The metallicity prior can be

1. Flat

.. code-block:: python

  met_prior_model = {'name': 'flat'}

Plot showing examples of the possible mass prior models with the parameters given above.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from beast.physicsmodel.prior_weights_stars import compute_metallicity_prior_weights

    fig, ax = plt.subplots()

    # met grid with linear spacing
    mets = np.linspace(0.004, 0.03)

    met_prior_models = [{"name": "flat"},]

    for mp_mod in met_prior_models:
        ax.plot(mets, compute_metallicity_prior_weights(mets, mp_mod), label=mp_mod["name"])

    ax.set_ylabel("probability")
    ax.set_xlabel("metallicity")
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()


Extinction
==========

A(V)
----

The A(V) prior can be:

1. Flat

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
                    'max_pos2': 2.0,
                    'sigma1': 1.0,
                    'sigma2': 0.2,
                    'N1': 20.,
                    'N2': 50.}

4. Exponential with decay rate 'a' and amplitude 'N'

.. code-block:: python

  av_prior_model = {'name': 'exponential',
                    'a': 1.0,
                    'N': 10.}

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from beast.physicsmodel.prior_weights_dust import PriorWeightsDust

    fig, ax = plt.subplots()

    # met grid with linear spacing
    avs = np.linspace(0.0, 10.0, num=200)

    dust_prior_models = [
        {"name": "flat"},
        {"name": "lognormal", "max_pos": 2.0, "sigma": 1.0, "N": 10.0},
        {
            "name": "two_lognormal",
            "max_pos1": 0.2,
            "max_pos2": 2.0,
            "sigma1": 1.0,
            "sigma2": 0.5,
            "N1": 20.,
            "N2": 50.,
        },
        {"name": "exponential", "a": 1.0, "N": 10.0},
    ]

    for dmod in dust_prior_models:
        dmodel = PriorWeightsDust(
            avs, dmod, [1.0], {"name": "flat"}, [1.0], {"name": "flat"}
        )

        ax.plot(avs, dmodel.av_priors, label=dmod["name"])

    ax.set_ylabel("probability")
    ax.set_xlabel("A(V)")
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()

R(V)
----

1. Flat

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

1. Flat

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
