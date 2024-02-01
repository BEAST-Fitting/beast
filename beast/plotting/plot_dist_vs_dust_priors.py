import numpy as np
import matplotlib.pyplot as plt

# from astropy.visualization import SqrtStretch, LogStretch, ImageNormalize
from astropy.modeling.models import Gaussian2D
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

    d1 = 50.0
    d2 = 70.0
    dist = np.arange(d1, d2, 0.5) * 1e3

    distmod = {
        "name": "absexponential",
        "dist0": 60.0 * u.kpc,
        "tau": 1.0 * u.kpc,
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

    av1 = 0.0
    av2 = 2.0
    avs = np.arange(av1, av2, 0.025)

    distim, avim = np.meshgrid(dist, avs)

    # generate a 2D Gaussian for the distribtion on sky
    skyprior = Gaussian2D(
        amplitude=1.0, x_mean=0.5, y_mean=0.5, x_stddev=0.01, y_stddev=0.01
    )
    x = np.arange(0.0, 1.0, 0.025)
    y = np.arange(0.0, 1.0, 0.025)
    alldists = None
    sumprobim = distim * 0.0

    for xi in x:
        for yi in y:
            avmod = {
                "name": "step",
                "dist0": 60.0 * u.kpc,
                "amp1": 0.1,
                "damp2": skyprior(xi, yi),
                "lgsigma1": 0.05,
                "lgsigma2": 0.05,
            }
            avprior = PriorDustModel(avmod)
            probim = avprior(avim, y=distim)
            probim /= np.sum(probim)
            sumprobim += probim * distprior(distim)

            # for visualization of result - only show amp
            av_vis = dist * 0.0 + 0.1
            av_vis[dist >= 60e3] = skyprior(xi, yi) + 0.1

            ax[0, 0].plot(dist / 1e3, av_vis, "k-", alpha=0.1)
            ax[0, 0].set_xlabel("distance [kpc]")
            ax[0, 0].set_ylabel("A(V) prior")
            ax[0, 0].set_ylim(0.0, 2.0)

            # sample from the priors
            # distsamp = np.interp(rng.random(npts), distsum, dist)
            # avsamp = np.interp(distsamp, dist, avprior(dist))
            # if alldists is None:
            #     alldists = distsamp
            #     allavs = avsamp
            # else:
            #     alldists = np.concatenate((alldists, distsamp))
            #     allavs = np.concatenate((allavs, avsamp))

    # ax[0, 1].hist2d(alldists, allavs, bins=20, norm="log")
    # norm = ImageNormalize(vmin=1e-5, vmax=1, stretch=LogStretch())
    ax[0, 1].imshow(
        sumprobim, origin="lower", aspect="auto", extent=[d1, d2, av1, av2], norm="log"
    )
    ax[0, 1].set_xlabel("distance [kpc]")
    ax[0, 1].set_ylabel("A(V)")

    # display the on sky distribution
    imx, imy = np.meshgrid(x, y)
    image = skyprior(imx, imy) + 0.1
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
