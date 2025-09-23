import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
import astropy.units as u

from beast.physicsmodel.grid_weights import compute_bin_boundaries
import beast.physicsmodel.priormodel_functions as pmfuncs


__all__ = [
    "PriorModel",
    "PriorDustModel",
    "PriorAgeModel",
    "PriorMassModel",
    "PriorMetallicityModel",
    "PriorDistanceModel",
]


class PriorModel:
    """
    Compute the priors as weights given the input grid
    """

    def __init__(self, model, allowed_models=None):
        """
        Initialize with basic information

        Parameters
        ----------
        model: dict
          Choice of model type
        """
        if (allowed_models is not None) and (model["name"] not in allowed_models):
            modname = model["name"]
            raise NotImplementedError(f"{modname} is not an allowed model")
        # save the model
        self.model = model

    def __call__(self, x, y=None):
        """
        Weights based on input model choice

        Parameters
        ----------
        x : float
            values for model evaluation
        y : float
            secondary values for 2D priors
        """
        if self.model["name"] == "flat":
            if "amp" in self.model.keys():
                amp = self.model["amp"]
            else:
                amp = 1.0
            if hasattr(x, "shape"):
                return np.full(x.shape, amp)
            else:
                return amp
        elif self.model["name"] in ["bins_histo", "bins_interp"]:
            for ckey in ["x", "values"]:
                if ckey not in self.model.keys():
                    raise KeyError(f"{ckey} not in prior model keys")
            # check if all ages within interpolation range
            mod_x = self.model["x"]
            if np.any([(cval > np.max(mod_x)) or (cval < np.min(mod_x)) for cval in x]):
                raise ValueError("requested x outside of model x range")

            if self.model["name"] == "bins_histo":
                # interpolate according to bins, assuming value is constant from i to i+1
                # and allow for bin edges input
                if len(self.model["values"]) == len(self.model["x"]) - 1:
                    self.model["values"].append(0.0)
                interfunc = interp1d(self.model["x"], self.model["values"], kind="zero")
                return interfunc(x)
            else:
                # interpolate model to grid ages
                return np.interp(
                    x, np.array(self.model["x"]), np.array(self.model["values"]),
                )
        elif self.model["name"] == "lognormal":
            for ckey in ["mean", "sigma"]:
                if ckey not in self.model.keys():
                    raise ValueError(f"{ckey} not in prior model keys")
            return pmfuncs._lognorm(x, self.model["mean"], sigma=self.model["sigma"])
        elif self.model["name"] == "two_lognormal":
            for ckey in ["mean1", "sigma1", "mean2", "sigma2"]:
                if ckey not in self.model.keys():
                    raise ValueError(f"{ckey} not in prior model keys")
            return pmfuncs._two_lognorm(
                x,
                self.model["mean1"],
                self.model["mean2"],
                sigma1=self.model["sigma1"],
                sigma2=self.model["sigma2"],
                N1=self.model["N1_to_N2"],
                N2=1.0,
            )
        elif self.model["name"] == "absexponential":
            for ckey in ["dist0", "tau", "amp"]:
                if ckey not in self.model.keys():
                    raise ValueError(f"{ckey} not in prior model keys")
            return pmfuncs._absexponential(
                x,
                dist0=self.model["dist0"].to(u.pc).value,
                tau=self.model["tau"].to(u.pc).value,
                amp=self.model["amp"],
            )
        elif self.model["name"] == "step":
            for ckey in ["dist0", "amp1", "damp2", "lgsigma1", "lgsigma2"]:
                if ckey not in self.model.keys():
                    raise ValueError(f"{ckey} not in prior model keys")
            if y is None:
                raise ValueError("y values not passed required for 2D priors")
            if len(x) != len(y):
                raise ValueError(
                    "x and y values not the same length, required for 2D priors"
                )
            return pmfuncs._step(
                x,
                y,
                dist0=self.model["dist0"].to(u.pc).value,
                amp1=self.model["amp1"],
                damp2=self.model["damp2"],
                lgsigma1=self.model["lgsigma1"],
                lgsigma2=self.model["lgsigma2"],
            )
        else:
            modname = self.model["name"]
            raise NotImplementedError(f"{modname} is not an allowed model")


class PriorDustModel(PriorModel):
    """
    Prior model for dust parameters with specific allowed models.
    """

    def __init__(self, model):
        """
        Initialize the dust prior model

        Parameters
        ----------
        model : dict
          Possible choices are flat, lognormal, two_lognormal, and exponential
        """
        super().__init__(
            model, allowed_models=["flat", "lognormal", "two_lognormal", "step"]
        )


class PriorAgeModel(PriorModel):
    """
    Prior model for age parameter with specific allowed models.
    """

    def __init__(self, model):
        """
        Initialize the stellar age prior model

        Parameters
        ----------
        model : dict
          Possible choices are flat, flat_log, bins_histo, bins_interp, and exponential
        """
        super().__init__(
            model,
            allowed_models=[
                "flat",
                "flat_log",
                "bins_histo",
                "bins_interp",
                "exponential",
            ],
        )

    def __call__(self, x):
        """
        Weights based on input model choice

        Parameters
        ----------
        x : float
            values for model evaluation
        """
        if self.model["name"] == "flat_log":
            weights = 1.0 / np.diff(10 ** compute_bin_boundaries(x))
            return weights / np.sum(weights)
        elif self.model["name"] == "exponential":
            return pmfuncs._exponential(10.0 ** x, tau=self.model["tau"] * 1e9)
        else:
            return super().__call__(x)


class PriorMetallicityModel(PriorModel):
    """
    Prior model for metallicity parameter with specific allowed models.
    """

    def __init__(self, model):
        """
        Initialize the metallicity prior model

        Parameters
        ----------
        model : dict
          Possible choices are flat
        """
        super().__init__(model, allowed_models=["flat"])


class PriorDistanceModel(PriorModel):
    """
    Prior model for distance parameter with specific allowed models.
    """

    def __init__(self, model):
        """
        Initialize the distance prior model

        Parameters
        ----------
        model : dict
          Possible choices are flat
        """
        super().__init__(model, allowed_models=["flat", "absexponential"])


class PriorMassModel(PriorModel):
    """
    Prior model for mass parameter with specific allowed models.
    """

    def __init__(self, model):
        """
        Initialize the initial mass prior model

        Parameters
        ----------
        model : dict
          Possible choices are flat, slapeter, and kroupa
        """
        super().__init__(model, allowed_models=["flat", "salpeter", "kroupa"])

    def __call__(self, x):
        """
        Weights based on input model choice

        Parameters
        ----------
        x : float
            values for model evaluation
        """
        # sort the initial mass along this isochrone
        sindxs = np.argsort(x)

        # Compute the mass bin boundaries
        mass_bounds = compute_bin_boundaries(x[sindxs], noneg=True)

        # integrate the IMF over each bin
        args = None
        if self.model["name"] == "kroupa":
            if "alpha0" in self.model.keys():  # assume other alphas also present
                args = (
                    self.model["alpha0"],
                    self.model["alpha1"],
                    self.model["alpha2"],
                    self.model["alpha3"],
                )
            imf_func = pmfuncs._imf_kroupa
        elif self.model["name"] == "salpeter":
            if "slope" in self.model.keys():
                slope = self.model["slope"]
                args = (slope,)
            imf_func = pmfuncs._imf_salpeter
        elif self.model["name"] == "flat":
            imf_func = pmfuncs._imf_flat

        # calculate the average prior in each mass bin
        mass_weights = np.zeros(len(x))
        for i, cindx in enumerate(sindxs):
            # fmt: off
            if args is not None:
                mass_weights[cindx] = (quad(imf_func, mass_bounds[i], mass_bounds[i + 1], args))[0]
            else:
                mass_weights[cindx] = (quad(imf_func, mass_bounds[i], mass_bounds[i + 1]))[0]
            # fmt: on
            mass_weights[cindx] /= mass_bounds[i + 1] - mass_bounds[i]

        # normalize to avoid numerical issues (too small or too large)
        mass_weights /= np.average(mass_weights)

        return mass_weights
