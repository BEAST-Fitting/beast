"""
Make a plot of the evolutionary tracks
"""

import matplotlib.pyplot as plt

from beast.physicsmodel.model_grid import make_evoltrack_table

from beast.plotting.beastplotlib import initialize_parser


if __name__ == "__main__":
    parser = initialize_parser()
    args = parser.parse_args()

    fig, ax = plt.subplots(2, 2, figsize=(10, 10))

    # download the file live from the website
    savename = "/tmp/padova_iso.csv"
    iso_fname, oiso = make_evoltrack_table(
        "test", oet=savename, age_info=[6.0, 10.13, 0.15], z=[0.019]
    )

    # plot, plot, plot
    nls = ":"
    oiso.plot(ax[0, 0], xval="logT", yval="logL", linestyle=nls)
    oiso.plot(ax[1, 1], xval="logA", yval="M_ini", linestyle=nls)
    oiso.plot(ax[0, 1], xval="stage", yval="M_ini", linestyle=nls)
    oiso.plot(ax[1, 0], xval="logT", yval="logg", linestyle=nls)

    print("size grid = ", len(oiso.data["M_ini"]))

    grid_metrics = oiso.grid_metrics()
    for ckey in grid_metrics.keys():
        print(ckey)
        print(grid_metrics[ckey])

    fig.tight_layout()

    save_name = "evoltracks"
    if args.tex:
        plt.rc({"usetex": True})
    if args.savefig:
        fig.savefig("{}.{}".format(save_name, args.savefig))
    else:
        plt.show()
