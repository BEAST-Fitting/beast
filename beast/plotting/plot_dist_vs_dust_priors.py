import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling.models import Gaussian2D
from astropy.visualization import simple_norm
import astropy.units as u

from beast.plotting.beastplotlib import initialize_parser
from beast.physicsmodel.priormodel import PriorDustModel, PriorDistanceModel

if __name__ == "__main__":  # pragma: no cover
    parser = initialize_parser()
    parser.add_argument(
        "--save_name",
        action="store",
        default="dist_vs_dust_priors_example",
        help="Save figure to file",
    )
    args = parser.parse_args()

    fig, ax = plt.subplots(2, 2, figsize=(10, 8))

    dist = np.arange(50.0, 70.0, 0.05) * 1e3

    distmod = {
        "name": "absexponential",
        "dist_0": 60.0 * u.kpc,
        "tau": 5.0 * u.kpc,
        "amp": 1.0,
    }
    # distmod = {"name": "flat"}
    distprior = PriorDistanceModel(distmod)
    ax[1, 0].plot(dist, distprior(dist), "k-")
    ax[1, 0].set_xlabel("distance [pc]")
    ax[1, 0].set_ylabel("# stars prior")

    rng = np.random.default_rng()
    npts = 100
    distsum = np.cumsum(distprior(dist))
    distsum /= distsum[-1]

    # generate a 2D Gaussian for the distribtion on sky
    skyprior = Gaussian2D(
        amplitude=1.0, x_mean=0.5, y_mean=0.5, x_stddev=0.2, y_stddev=0.2
    )
    x = np.arange(0.0, 1.0, 0.03)
    y = np.arange(0.0, 1.0, 0.03)
    alldists = None
    for i in range(len(x)):
        for j in range(len(y)):
            avmod = {
                "name": "step",
                "dist_0": 60.0 * u.kpc,
                "amp_1": 0.0,
                "amp_2": skyprior(x[i], y[j]),
            }
            avprior = PriorDustModel(avmod)
            ax[0, 0].plot(dist, avprior(dist), "k-", alpha=0.1)
            ax[0, 0].set_xlabel("distance [pc]")
            ax[0, 0].set_ylabel("A(V) prior")

            # sample from the priors
            distsamp = np.interp(rng.random(npts), distsum, dist)
            avsamp = np.interp(distsamp, dist, avprior(dist))
            if alldists is None:
                alldists = distsamp
                allavs = avsamp
            else:
                alldists = np.concatenate((alldists, distsamp))
                allavs = np.concatenate((allavs, avsamp))

    ax[0, 1].hist2d(alldists, allavs, bins=20, norm="log")
    ax[0, 1].set_xlabel("distance [pc]")
    ax[0, 1].set_ylabel("A(V) samples")

    # display the on sky distribution
    imx, imy = np.meshgrid(x, y)
    image = skyprior(imx, imy)
    # norm = simple_norm(image, 'sqrt')
    ax[1, 1].imshow(image, origin="lower")
    ax[1, 1].set_xlabel("x")
    ax[1, 1].set_ylabel("y")
    ax[1, 1].set_title("sky distribution")

    plt.tight_layout()

    if args.tex:
        plt.rc({"usetex": True})
    if args.savefig:
        fig.savefig("{}.{}".format(args.save_name, args.savefig))
    else:
        plt.show()
