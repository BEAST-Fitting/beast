"""
Make a plot of the evolutionary tracks
"""
import matplotlib.pyplot as plt

from beast.physicsmodel.stars.evoltracks import ETParsec, ETMist
from beast.plotting.beastplotlib import initialize_parser


if __name__ == "__main__":
    parser = initialize_parser()
    parser.add_argument(
        "filename", type=str, nargs="*", help="name of file(s) with evolutionary track"
    )
    parser.add_argument(
        "--type", default="mist", choices=["mist", "parsec"], help="source of tracks"
    )
    parser.add_argument(
        "--condense",
        help="condense grid based on logM, logT, logL deltas",
        action="store_true",
    )
    args = parser.parse_args()

    fig, ax = plt.subplots(2, 2, figsize=(10, 10))

    # switch between the types of tracks
    if args.type == "parsec":
        et = ETParsec()
    else:
        et = ETMist()

    # read in the tracks
    et.data = et.load_orig_tables(filename=args.filename)

    orig_metrics = et.grid_metrics()

    # plot, plot, plot
    nls = ":"
    color = "k"
    et.plot_tracks(ax[0, 0], xval="logT", yval="logL", linestyle=nls, color=color)
    et.plot_tracks(ax[1, 1], xval="logA", yval="M_ini", linestyle=nls, color=color)
    et.plot_tracks(ax[0, 1], xval="eep", yval="M_ini", linestyle=nls, color=color)
    et.plot_tracks(ax[1, 0], xval="logT", yval="logg", linestyle=nls, color=color)

    # regrid the evolutionary tracks to uniform log(mass) and variable age
    print("size orig = ", len(et.data["M_ini"]))

    logmass_range = [-1.0, 3.0]
    condense = args.condense
    logmass_delta = 0.05
    logT_delta = 0.05
    logL_delta = 0.05
    et.regrid_one_met(
        logmass_range=logmass_range,
        condense=condense,
        logmass_delta=logmass_delta,
        logT_delta=logT_delta,
        logL_delta=logL_delta,
    )

    print("size regrid = ", len(et.data["M_ini"]))

    print("logM range:", logmass_range)
    print("logM:", logmass_delta)

    print("condense:", condense)
    if condense:
        print("condense deltas")
        print("logT:", logT_delta)
        print("logL:", logL_delta)

    regrid_metrics = et.grid_metrics()
    print("grid metrics: diffs [min,  max, median, mean]")
    for ckey in regrid_metrics.keys():
        print(ckey)
        print("  orig:", orig_metrics[ckey])
        print("regrid:", regrid_metrics[ckey])

    # get the new grid metrics
    # et.grid_metrics()

    # plot, plot, plot
    nls = "-"
    color = "b"
    et.plot_tracks(ax[0, 0], xval="logT", yval="logL", linestyle=nls, color=color)
    et.plot_tracks(ax[1, 1], xval="logA", yval="M_ini", linestyle=nls, color=color)
    et.plot_tracks(ax[0, 1], xval="eep", yval="M_ini", linestyle=nls, color=color)
    et.plot_tracks(ax[1, 0], xval="logT", yval="logg", linestyle=nls, color=color)

    fig.tight_layout()

    save_name = "evoltracks"
    if args.tex:
        plt.rc({"usetex": True})
    if args.savefig:
        fig.savefig("{}.{}".format(save_name, args.savefig))
    else:
        plt.show()
