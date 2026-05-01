import numpy as np
import matplotlib.pyplot as plt

from beast.physicsmodel.stars import stellib
from helpers import get_stellib_boundaries

if __name__ == "__main__":  # pragma: no cover

    slib = stellib.Tlusty2025()
    bound = get_stellib_boundaries(slib, dlogT=0.01, dlogg=0.01)
    bbox = np.array(slib.bbox())
    print(bbox)
    print(bound)

    plt.plot(slib.logT, slib.logg, "ko")
    plt.plot(bound[:, 0], bound[:, 1])
    plt.plot(bbox[:, 0], bbox[:, 1])
    plt.show()
