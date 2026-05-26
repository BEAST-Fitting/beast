import numpy as np
import matplotlib.pyplot as plt

from beast.physicsmodel.stars import stellib

if __name__ == "__main__":  # pragma: no cover

    slibs = [stellib.Tlusty2025(),
             stellib.BOSZ2024(),
             stellib.Aringer2016()]
    labels = ["Tlusty2025",
              "BOSZ2025",
              "Aringer2016"]
    # slibs = [stellib.Tlusty(),
    #          stellib.Kurucz()]
    # labels = ["Tlusty",
    #           "Kurucz"]
    colors = ["r", "b", "g"]

    plt.plot([1.0, 1.0], [1.1, 1.1])
    for slib, ccol, clab in zip(slibs, colors, labels):
        slib.plot_boundary(dlogT=0.0, dlogg=0.0, color=ccol, label=clab, alpha=0.5)

    plt.xlim(3.3, 4.8)
    plt.ylim(-1.5, 6.0)
    plt.xlabel("log(Teff)")
    plt.ylabel("log(g)")
    plt.legend()
    plt.show()