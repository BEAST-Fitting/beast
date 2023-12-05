.. _beast_priors:

######
Priors
######

Priors on the parameters are set in the beast_settings.txt file using
python dictionaries to give the prior model and any related
parameters.

All priors are normalized to have an average value of 1.  This is possible
as the priors are relevant in a relative sense, not absolute and this
avoids numerical issues with too large or small numbers.

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

  age_prior_model = {"name": "flat"}

or to set the star formation rate in M_sun/year, use

.. code-block:: python

  age_prior_model = {"name": "flat"},
                     "amp": sfr}

2. Flat in log age

.. code-block:: python

  age_prior_model = {"name": "flat_log"}

3. Set by bins spaced in logage (log10(years)).

For example, step like priors can be specified by:

.. code-block:: python

  age_prior_model = {"name": "bins_histo",
                     "x": [6.0, 7.0, 8.0, 9.0, 10.0],
                     "values": [1.0, 2.0, 1.0, 5.0, 3.0]}

Or using bin edges (where N = N_values+1) like those output by `np.histogram()`:

.. code-block:: python

  age_prior_model = {"name": "bins_histo",
                     "x": [6.0, 7.0, 8.0, 9.0, 10.0],
                     "values": [1.0, 2.0, 1.0, 5.0]}

For example, lines connecting the bin value of the priors can be specified by:

.. code-block:: python

 age_prior_model = {"name": "bins_interp",
                    "x": [6.0, 7.0, 8.0, 9.0, 10.0],
                    "values": [1.0, 2.0, 1.0, 5.0, 3.0]}

4. An exponentially decreasing SFR (in time, but here increasing with age)
with a 0.1 Gyr time constant (with `tau` parameter in Gyr):

.. code-block:: python

  age_prior_model = {"name": "exponential",
                     "tau": 0.1,
                     "amp": 1.0}

Plot showing examples of the possible age prior models with the parameters given above.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from beast.physicsmodel.priormodel import PriorAgeModel

    fig, ax = plt.subplots()

    # logage grid from 1 Myrs to 10 Gyrs
    logages = np.linspace(6.0, 10.0)

    age_prior_models = [
        {"name": "flat"},
        {"name": "flat_log"},
        {
            "name": "bins_histo",
            "x": [6.0, 7.0, 8.0, 9.0, 10.0],
            "values": [1.0, 2.0, 1.0, 5.0, 3.0],
        },
        {
            "name": "bins_interp",
            "x": [6.0, 7.0, 8.0, 9.0, 10.0],
            "values": [1.0, 2.0, 1.0, 5.0, 3.0],
        },
        {"name": "exponential", "tau": 0.1, "amp": 1.0}
    ]

    for ap_mod in age_prior_models:
        pmod = PriorAgeModel(ap_mod)
        ax.plot(logages, pmod(logages), label=ap_mod["name"])

    ax.set_ylabel("probability")
    ax.set_xlabel("log(age)")
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()


Mass
----

The mass prior is set by the choice of an Initial Mass Function (IMF).
The mass function supported are:

1. Kroupa

Functional form from Kroupa (2001, MNRAS, 322, 231) with alpha0,1,2,3 slopes
for <0.08, 0.08-0.5, 0.5-1.0, >1.0 solar masses.

No alpha0,1,2,3 values gives the defaults listed in 2nd example.

.. code-block:: python

  mass_prior_model = {"name": "kroupa"}

With explicit values for the alphas (all need to be specified).

.. code-block:: python

  mass_prior_model = {"name": "kroupa",
                      "alpha0": 0.3,
                      "alpha1": 1.3,
                      "alpha2": 2.3,
                      "alpha3": 2.3}


2. Salpeter

Functional form from Salpeter (1955, ApJ, 121, 161).

No slope value gives the default listed in 2nd example.

.. code-block:: python

  mass_prior_model = {"name": "salpeter"}

With an explicit value for the slope.

.. code-block:: python

  mass_prior_model = {"name": "salpeter",
                      "slope": 2.35}

3. Flat

There is also a flat mass prior.  This is useful for creating grids for BEAST
verification (see :doc:`Simulations <simulations>`), and should not be
used for a standard fitting run.

.. code-block:: python

  mass_prior_model = {"name": "flat"}


Plot showing examples of the possible mass prior models with the parameters given above.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from beast.physicsmodel.priormodel import PriorMassModel

    fig, ax = plt.subplots()

    # mass grid from 0.01 to 100 solar masses (log spacing)
    masses = np.logspace(-2.0, 2.0)

    mass_prior_models = [
        {"name": "kroupa"},
        {"name": "salpeter"},
        {"name": "flat"}
    ]

    for mp_mod in mass_prior_models:
        pmod = PriorMassModel(mp_mod)
        ax.plot(masses, pmod(masses), label=mp_mod["name"])

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

  met_prior_model = {"name": "flat"}

Plot showing examples of the possible metallicity prior models with the parameters given above.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from beast.physicsmodel.priormodel import PriorMetallicityModel

    fig, ax = plt.subplots()

    # met grid with linear spacing
    mets = np.linspace(0.004, 0.03)

    met_prior_models = [{"name": "flat"}]

    for mp_mod in met_prior_models:
        pmod = PriorMetallicityModel(mp_mod)
        ax.plot(mets, pmod(mets), label=mp_mod["name"])

    ax.set_ylabel("probability")
    ax.set_xlabel("metallicity")
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()

Distance
--------

The distance prior can be

1. Flat

.. code-block:: python

  distance_prior_model = {"name": "flat"}

2. Absolute(Exponential) distribution with an exponential scale height (tau) before and
after a fiducial distance (dist0) and an amplitude (amp).

.. code-block:: python

  distance_prior_model = {"name": "absexponential",
                          "dist0": 60.0*u.kpc,
                          "tau": 5.*u.kpc,
                          "amp": 1.0}

Plot showing examples of the possible distance prior models with the parameters given above.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u

    from beast.physicsmodel.priormodel import PriorDistanceModel

    fig, ax = plt.subplots()

    # met grid with linear spacing
    dists = np.arange(50., 70, 0.1) * 1e3

    met_prior_models = [
      {"name": "flat"},
      {"name": "absexponential",
        "dist0": 60.0*u.kpc,
        "tau": 5.*u.kpc,
        "amp": 1.0}]

    for mp_mod in met_prior_models:
        pmod = PriorDistanceModel(mp_mod)
        ax.plot(dists, pmod(dists), label=mp_mod["name"])

    ax.set_ylabel("probability")
    ax.set_xlabel("distance")
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

  av_prior_model = {"name": "flat"}

2. Lognormal with the maximum at the A(V) given by mean and the width
given by sigma.

.. code-block:: python

  av_prior_model = {"name": "lognormal",
                    "mean": 2.0,
                    "sigma": 1.0}

3. Two lognormals (see above for definition of terms)

.. code-block:: python

  av_prior_model = {"name": "two_lognormal",
                    "mean1": 0.2,
                    "mean2": 2.0,
                    "sigma1": 1.0,
                    "sigma2": 0.2,
                    "N1_to_N2": 1.0 / 5.0}

4. Step at a specified distance.  Distance must have units.  Models
the effect of having a dust cloud located at a certain distance.
A(V) after dist0 is amp1 + damp2.

.. code-block:: python

  av_prior_model = {"name": "step",
                    "dist0": 60 * u.kpc,
                    "amp1": 0.1,
                    "damp2": 1.0,
                    "lgsigma1": 0.05,
                    "lgsigma2": 0.05}

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from beast.physicsmodel.priormodel import PriorDustModel

    fig, ax = plt.subplots()

    # av grid with linear spacing
    avs = np.linspace(0.0, 10.0, num=200)

    dust_prior_models = [
        {"name": "flat"},
        {"name": "lognormal", "mean": 2.0, "sigma": 1.0},
        {
            "name": "two_lognormal",
            "mean1": 0.2,
            "mean2": 2.0,
            "sigma1": 1.0,
            "sigma2": 0.5,
            "N1_to_N2": 1.0 / 5.0
        },
    ]

    for dmod in dust_prior_models:
        pmod = PriorDustModel(dmod)
        ax.plot(avs, pmod(avs), label=dmod["name"])

    ax.set_ylabel("probability")
    ax.set_xlabel("A(V)")
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u

    from beast.physicsmodel.priormodel import PriorDustModel

    fig, ax = plt.subplots()

    # distance grid with linear spacing
    d1, d2 = (50.e3, 70.e3)
    dists = np.linspace(d1, d2, num=100)
    av1, av2 = (0.0, 2.0)
    avs = np.arange(av1, av2, 0.025)
    distim, avim = np.meshgrid(dists, avs)

    dustmod = {
        "name": "step",
        "dist0": 60 * u.kpc,
        "amp1": 0.1,
        "damp2": 1.0,
        "lgsigma1": 0.05,
        "lgsigma2": 0.05}

    dustprior = PriorDustModel(dustmod)
    probim = dustprior(avim, y=distim)

    ax.imshow(
        probim, origin="lower", aspect="auto", extent=[d1, d2, av1, av2], norm="log"
    )

    ax.set_ylabel("A(V)")
    ax.set_xlabel("distance [kpc]")
    ax.set_title("step")
    plt.tight_layout()
    plt.show()

R(V)
----

1. Flat

.. code-block:: python

  rv_prior_model = {"name": "flat"}

2. Lognormal with the maximum at the R(V) given by mean and the width
given by sigma.

.. code-block:: python

  rv_prior_model = {"name": "lognormal",
                    "mean": 3.1,
                    "sigma": 0.25}

3. Two lognormals (see above for definition of terms)

.. code-block:: python

  rv_prior_model = {"name": "two_lognormal",
                    "mean1": 3.1,
                    "mean1": 4.5,
                    "sigma1": 0.1,
                    "sigma2": 0.2,
                    "N1_to_N2": 2.0 / 5.0}

4. Step at a specified distance.  Distance must have units.  Models
the effect of having a dust cloud located at a certain distance.
R(V) after dist0 is amp1 + damp2.

.. code-block:: python

  rv_prior_model = {"name": "step",
                    "dist0": 60 * u.kpc,
                    "amp1": 0.1,
                    "damp2": 1.0,
                    "lgsigma1": 0.05,
                    "lgsigma2": 0.05}

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from beast.physicsmodel.priormodel import PriorDustModel

    fig, ax = plt.subplots()

    # rv grid with linear spacing
    rvs = np.linspace(2.0, 6.0, num=200)

    dust_prior_models = [
        {"name": "flat"},
        {"name": "lognormal", "mean": 3.1, "sigma": 0.25},
        {
            "name": "two_lognormal",
            "mean1": 3.1,
            "mean2": 4.5,
            "sigma1": 0.1,
            "sigma2": 0.2,
            "N1_to_N2": 2.0 / 5.0
        }
    ]

    for dmod in dust_prior_models:
        pmod = PriorDustModel(dmod)
        ax.plot(rvs, pmod(rvs), label=dmod["name"])

    ax.set_ylabel("probability")
    ax.set_xlabel("R(V)")
    ax.set_title("step")
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u

    from beast.physicsmodel.priormodel import PriorDustModel

    fig, ax = plt.subplots()

    # distance grid with linear spacing
    d1, d2 = (50.e3, 70.e3)
    dists = np.linspace(d1, d2, num=100)
    rv1, rv2 = (2.0, 6.0)
    rvs = np.arange(rv1, rv2, 0.05)
    distim, rvim = np.meshgrid(dists, rvs)

    dustmod = {
        "name": "step",
        "dist0": 60 * u.kpc,
        "amp1": 3.1,
        "damp2": 1.4,
        "lgsigma1": 0.01,
        "lgsigma2": 0.01}

    dustprior = PriorDustModel(dustmod)
    probim = dustprior(rvim, y=distim)

    ax.imshow(
        probim, origin="lower", aspect="auto", extent=[d1, d2, rv1, rv2], norm="log"
    )

    ax.set_ylabel("R(V)")
    ax.set_xlabel("distance [kpc]")
    ax.set_title("step")
    plt.tight_layout()
    plt.show()

f_A
---

1. Flat

.. code-block:: python

  fA_prior_model = {"name": "flat"}

2. Lognormal with the maximum at the f_A given by mean and the width
given by sigma.

.. code-block:: python

  fA_prior_model = {"name": "lognormal",
                    "mean": 0.8,
                    "sigma": 0.1}

3. Two lognormals (see above for definition of terms)

.. code-block:: python

  fA_prior_model = {"name": "two_lognormal",
                    "mean1": 0.1,
                    "mean1": 0.8,
                    "sigma1": 0.1,
                    "sigma2": 0.2,
                    "N1_to_N2": 2.0 / 5.0}

4. Step at a specified distance.  Distance must have units.  Models
the effect of having a dust cloud located at a certain distance.
f_A after dist0 is amp1 + damp2.

.. code-block:: python

  fA_prior_model = {"name": "step",
                    "dist0": 60 * u.kpc,
                    "amp1": 0.1,
                    "damp2": 0.8,
                    "lgsigma1": 0.1,
                    "lgsigma2": 0.01}

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from beast.physicsmodel.priormodel import PriorDustModel

    fig, ax = plt.subplots()

    # fA grid with linear spacing
    fAs = np.linspace(0.0, 1.0, num=200)

    dust_prior_models = [
        {"name": "flat"},
        {"name": "lognormal", "mean": 0.8, "sigma": 0.1},
        {
            "name": "two_lognormal",
            "mean1": 0.2,
            "mean2": 0.8,
            "sigma1": 0.1,
            "sigma2": 0.2,
            "N1_to_N2": 2.0 / 5.0
        }
    ]

    for dmod in dust_prior_models:
        pmod = PriorDustModel(dmod)
        ax.plot(fAs, pmod(fAs), label=dmod["name"])

    ax.set_ylabel("probability")
    ax.set_xlabel(r"$f_A$")
    ax.legend(loc="best")
    plt.tight_layout()
    plt.show()

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u

    from beast.physicsmodel.priormodel import PriorDustModel

    fig, ax = plt.subplots()

    # distance grid with linear spacing
    d1, d2 = (50.e3, 70.e3)
    dists = np.linspace(d1, d2, num=100)
    fA1, fA2 = (0.0, 1.0)
    fAs = np.arange(fA1, fA2, 0.01)
    distim, fAim = np.meshgrid(dists, fAs)

    dustmod = {
        "name": "step",
        "dist0": 60 * u.kpc,
        "amp1": 0.1,
        "damp2": 0.8,
        "lgsigma1": 0.1,
        "lgsigma2": 0.01}

    dustprior = PriorDustModel(dustmod)
    probim = dustprior(fAim, y=distim)

    ax.imshow(
        probim, origin="lower", aspect="auto", extent=[d1, d2, fA1, fA2], norm="log"
    )

    ax.set_ylabel(r"$f_A$")
    ax.set_xlabel("distance [kpc]")
    plt.tight_layout()
    plt.show()
