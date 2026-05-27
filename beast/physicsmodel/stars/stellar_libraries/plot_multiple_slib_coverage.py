import argparse
import matplotlib.pyplot as plt

from beast.physicsmodel.stars import stellib

if __name__ == "__main__":  # pragma: no cover

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--old",
        action="store_true",
        help="Show the old Tlusty/Kurucz grids",
    )
    parser.add_argument("--png", action="store_true", help="save figure to a png file")
    args = parser.parse_args()

    if args.old:
        slibs = [stellib.Tlusty(), stellib.Kurucz()]
        labels = ["Tlusty", "Kurucz"]
    else:
        slibs = [stellib.Tlusty2025(), stellib.BOSZ2024(), stellib.Aringer2016()]
        labels = ["Tlusty2025", "BOSZ2025", "Aringer2016"]
    colors = ["r", "b", "g"]

    # setup the plots
    fontsize = 12
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    plt.plot([1.0, 1.0], [1.1, 1.1])
    for slib, ccol, clab in zip(slibs, colors, labels):
        slib.plot_boundary(dlogT=0.0, dlogg=0.0, color=ccol, label=clab, alpha=0.5)

    plt.xlim(3.3, 4.8)
    plt.ylim(-1.5, 6.0)
    plt.xlabel("log(Teff)")
    plt.ylabel("log(g)")
    plt.legend()

    if args.old:
        extstr = "_old"
    else:
        extstr = ""
    if args.png:
        plt.savefig(f"slib_mulicoverage{extstr}.png")
    else:
        plt.show()
