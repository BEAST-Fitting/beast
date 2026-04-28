import numpy as np
import matplotlib.pyplot as plt

from beast.physicsmodel.stars import stellib


def get_stellib_boundaries(s, dlogT=0.1, dlogg=0.3, closed=True):
    """
    Returns the closed boundary polygon around the stellar library with
    given margins.

    Parameters
    ----------
    s :  Stellib
        Stellar library object
    dlogT :  float
        margin in logT
    dlogg  : float
        margin in logg
    closed : bool
        if set, close the polygon

    Returns
    ------
    b : ndarray[float, ndim=2]
        (closed) boundary points: [logg, Teff]
    """
    #        use "points_inside_poly" to test wether a point is inside the limits
    #        >>> data = np.array([iso.data['logg'], iso.data['logT']]).T
    #        >>> aa = points_inside_poly(data, leftb)
    #    """
    leftb = [(k, np.max(s.logT[s.logg == k]) + dlogT) for k in np.unique(s.logg)]
    leftb += [(leftb[-1][0] + dlogg, leftb[-1][1])]
    leftb = [(leftb[0][0] - dlogg, leftb[0][1])] + leftb
    rightb = [(k, np.min(s.logT[s.logg == k]) - dlogT) for k in np.unique(s.logg)[::-1]]
    rightb += [(rightb[-1][0] - dlogg, rightb[-1][1])]
    rightb = [(rightb[0][0] + dlogg, rightb[0][1])] + rightb
    b = leftb + rightb
    if closed:
        b += [b[0]]
    return np.array(b)


if __name__ == "__main__":  # pragma: no cover

    slib = stellib.Tlusty2025()
    bound = get_stellib_boundaries(slib, dlogT=0.01, dlogg=0.01)
    bbox = np.array(slib.bbox())
    print(bbox[:, 0])

    plt.plot(slib.logT, slib.logg, "ko")
    plt.plot(bound[:, 1], bound[:, 0])
    plt.plot(bbox[:, 0], bbox[:, 1])
    plt.show()
