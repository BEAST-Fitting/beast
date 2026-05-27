import argparse
import numpy as np
import matplotlib.pyplot as plt

from beast.physicsmodel.stars import stellib
from helpers import get_stellib_boundaries

if __name__ == "__main__":  # pragma: no cover

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--grid",
        default="BOSZ2024",
        choices=["Tlusty2025", "BOSZ2024", "Aringer2016", "Kurucz", "Tlusty"],
        help="Grid to show",
    )
    parser.add_argument("--png", action="store_true", help="save figures to png files")
    args = parser.parse_args()

    if args.grid == "BOSZ2024":
        slib = stellib.BOSZ2024()
    elif args.grid == "Tlusty2025":
        slib = stellib.Tlusty2025()
    elif args.grid == "Aringer2016":
        slib = stellib.Aringer2026()
    elif args.grid == "Kurucz":
        slib = stellib.Kurucz()
    elif args.grid == "Tlusty":
        slib = stellib.Tlusty()

    # useful when new grids are generarted
    # results are used to define the bounding box in a model definition
    bound = get_stellib_boundaries(slib, dlogT=0.0, dlogg=0.0)
    print("Algorithmically generated bounds")
    print(bound)

    # setup the plots
    fontsize = 12
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    for cz in np.unique(slib.Z.data):
        gvals = slib.Z.data == cz
        plt.plot(slib.logT[gvals], slib.logg[gvals], "ko")
        plt.plot(bound[:, 0], bound[:, 1], "g-", label="algorithmly generated")
        slib.plot_boundary(dlogT=0.0, dlogg=0.0, label="defined in class def", alpha=0.5, color="b")
        plt.title(f"{args.grid}; z = {cz:.2e}")
        plt.xlabel("log(Teff)")
        plt.ylabel("log(g)")
        plt.legend()

        if args.png:
            plt.savefig(f"slib_coverage_{args.grid}_z_{cz:.2e}.png")
            plt.close()
        else:
            plt.show()