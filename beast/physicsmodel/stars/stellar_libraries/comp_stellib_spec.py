import argparse
import matplotlib.pyplot as plt

from beast.physicsmodel.stars import stellib

if __name__ == "__main__":  # pragma: no cover

    grid_choices = ["Tlusty2025", "BOSZ2024", "Aringer2016", "Kurucz", "Tlusty"]

    parser = argparse.ArgumentParser()
    parser.add_argument("--logT", default=4.1, type=float,
                        help="comparision spectra log(Teff)")
    parser.add_argument("--logg", default=4.0, type=float,
                        help="comparision spectra log(g)")
    parser.add_argument("--met", default=0.015, type=float,
                        help="comparision spectra met")
    parser.add_argument(
        "--grid1",
        default="BOSZ2024",
        choices=grid_choices,
        help="Grid to show",
    )
    parser.add_argument(
        "--grid2",
        default="Tlusty2025",
        choices=grid_choices,
        help="Grid to show",
    )
    parser.add_argument("--png", action="store_true", help="save figures to png files")
    args = parser.parse_args()

    lib1_name = getattr(stellib, args.grid1)
    lib1 = lib1_name()
    lib2_name = getattr(stellib, args.grid2)
    lib2 = lib2_name()

    pts = {"logT" : [args.logT],
           "logg" : [args.logg],
           "Z" : [args.met],
           "logL" : [0.01]}  # dummy value

    spec1 = lib1.gen_spectral_grid_from_given_points(pts)
    spec2 = lib2.gen_spectral_grid_from_given_points(pts)

    if args.grid2 == "Aringer2016":
        div2 = 1e8
    else:
        div2 = 1.0

    # plot
    plt.plot(spec1.lamb, spec1.seds[0], label=args.grid1)
    plt.plot(spec2.lamb, spec2.seds[0] / div2, label=args.grid2)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel(r"$\lambda$ [$\AA$]")
    plt.ylabel(r"$F(\lambda)$")
    plt.legend()
    plt.title(f"log(T) = {args.logT:.2f}; log(g) = {args.logg:.2f}")
    plt.show()
