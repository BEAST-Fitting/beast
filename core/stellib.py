"""
Stellib class

Intent to implement a generic module to manage stellar library from various
sources.

The interpolation is impletmented from the pegase.2 fortran converted algorithm.
(this may not be pythonic though)

TODO: ugly redef of __interp__ to be cleaned up
"""
import numpy
#from numpy import interp
import pyfits
from ..external.eztables import Table
from . import grid
from ..config import __ROOT__, __WITH_C_LIBS__, __NTHREADS__

try:
    assert(__WITH_C_LIBS__), "Python code requested in config.py"
    from ..include.interp import __interp__ as __cinterp__

    def __interp__(T0, g0, T, g, dT_max=0.1, eps=1e-6):
        """
        Interpolation of the (T,g) grid at fixed Z

        Translated from Pegase.2 fortran version
        (this may not be pythonic though)

        Note: preference is always given to the temperature over
            the gravity when needed.

        Inputs:
            T0  log(Teff) to obtain
            g0  log(g) to obtain
            T   log(Teff) of the grid
            g   log(g) of the grid

        Keywords:
            dT_max  If, T2 (resp. T1) is too far from T compared to T1
                (resp. T2), i2 (resp. i1) is not used.
                (see below for namings)
            eps temperature sensitivity under which points are
                considered to have the same temperature

        Returns 4 star indexes and 4 associated weights

        if index is -1, this means the point is rejected and the
        associated weight is 0.
        """
        idx = numpy.zeros(4, dtype=numpy.int64)
        w   = numpy.zeros(4, dtype=numpy.float64)
        _T  = numpy.double(T)
        _g  = numpy.double(g)
        __cinterp__(T0, g0, idx, w, len(T), _T, _g, dT_max, eps)
        #__cinterp__(T0, g0, _T,_g, dT_max, eps)
        return idx, w

except Exception as e:
    #print "Using python code instead of c, because %s"  % e

    def __interp__(T0, g0, T, g, dT_max=0.1, eps=1e-6):
        """
        Interpolation of the (T,g) grid at fixed Z

        Translated from Pegase.2 fortran version
        (this may not be pythonic though)

        Note: preference is always given to the temperature over
            the gravity when needed.

        Inputs:
            T0  log(Teff) to obtain
            g0  log(g) to obtain
            T   log(Teff) of the grid
            g   log(g) of the grid

        Keywords:
            dT_max  If, T2 (resp. T1) is too far from T compared to T1
                (resp. T2), i2 (resp. i1) is not used.
                (see below for namings)
            eps temperature sensitivity under which points are
                considered to have the same temperature

        Returns 4 star indexes and 4 associated weights

        if index is -1, this means the point is rejected and the
        associated weight is 0.

        Naming:

        i1 = index of the star with temperature > T and gravity > g.
        Among all such stars, one chooses the one minimizing
        |Delta T|+kappa*|Delta g|.
        If no star with temperature > T and gravity > g exists, i1 = -1

        i2 = index of the star with temperature > T and gravity < g.

        i3 = index of the star with temperature < T and gravity > g.

        i4 = index of the star with temperature < T and gravity < g.

         g

        /|\
         | i3  |
         |     |  i1
         | ----x------
         |     |    i2
         |  i4 |
         |__________\ T
                /
        """
        kappa  = 0.1

        idx    = numpy.arange(len(g))
        deltag = g - g0
        deltaT = T - T0
        dist   = kappa * abs(deltag) + abs(deltaT)

        #TODO: check for update
        if dist.min() == 0:
            return (dist.argmin(), -1, -1, -1), (1., 0., 0., 0.)

        ## Looking for i_{1..4}
        # looking for i1
        ind = numpy.where( (deltag > 0.) & (deltaT > 0) )[0]
        if len(ind) == 0:
            i1 = -1
        else:
            i1  = idx[ind][dist[ind].argmin()]

        # looking for i2
        ind = numpy.where( (deltag < 0.) & (deltaT > 0) )[0]
        if len(ind) == 0:
            i2 = -1
        else:
            i2  = idx[ind][dist[ind].argmin()]

        # looking for i3
        ind = numpy.where( (deltag > 0.) & (deltaT < 0) )[0]
        if len(ind) == 0:
            i3 = -1
        else:
            i3  = idx[ind][dist[ind].argmin()]

        # looking for i4
        ind = numpy.where( (deltag < 0.) & (deltaT < 0) )[0]
        if len(ind) == 0:
            i4 = -1
        else:
            i4  = idx[ind][dist[ind].argmin()]

        if ( (i1 < 0) & (i2 < 0) & (i3 < 0) & (i4 < 0) ):
            assert(False), "Interp. Error"

        T1 = T[i1]
        T2 = T[i2]
        T3 = T[i3]
        T4 = T[i4]
        g1 = g[i1]
        g2 = g[i2]
        g3 = g[i3]
        g4 = g[i4]
        # If, T2 (resp. T1) is too far from T compared to T1
        # (resp. T2), i2 (resp. i1) is not used.
        # The same for i3 and i4.
        if ( (i1 > 0) & (i2 > 0) ):
            if (T1 < T2 - dT_max):
                i2 = -1
            elif (T2 < T1 - dT_max):
                i1 = -1

        if ( (i3 > 0) & (i4 > 0) ):
            if (T3 > T4 + dT_max):
                i4 = -1
            elif (T4 > T3 + dT_max):
                i3 = -1

        # Interpolation in the (T, g) plane between the used points
        # (at least 1, at most 4).
        # Code "0110" means that i1 = i4 = 0, i2 /=0 and i3 /= 0.
        #
        # Note: preference is always given to the temperature over
        #   the gravity when needed.

        if (i1 < 0):
            if (i2 < 0):
                if (i3 < 0):
                    if (i4 < 0):
                        #               #0000
                        assert (False), "Error"
                    else:                   # 0001
                        alpha1 = 0.
                        alpha2 = 0.
                        alpha3 = 0.
                        alpha4 = 1.
                    #endif
                elif (i4 < 0):                  # 0010
                    alpha1 = 0.
                    alpha2 = 0.
                    alpha3 = 1.
                    alpha4 = 0.
                else:                           # 0011
                    alpha1 = 0.
                    alpha2 = 0.
                    if ( abs(T3 - T4) < eps ):
                        if (g3 == g4):
                            alpha3 = 0.5
                        else:
                            alpha3 = (g0 - g4) / (g3 - g4)
                        #endif
                        alpha4 = 1. - alpha3
                    else:
                        if (T3 > T4):
                            alpha3 = 1.
                            alpha4 = 0.
                            i4 = -1
                        else:
                            alpha3 = 0.
                            i3 = -1
                            alpha4 = 1.
                        #endif
                    #endif
                #endif
            elif (i3 < 0):
                if (i4 < 0):
                    #                     #0100
                    alpha1 = 0.
                    alpha2 = 1.
                    alpha3 = 0.
                    alpha4 = 0.
                else:                         # 0101
                    alpha1 = 0.
                    if (T2 == T4):
                        alpha2 = 0.5
                    else:
                        alpha2 = (T0 - T4) / (T2 - T4)
                    #endif
                    alpha3 = 0.
                    alpha4 = 1. - alpha2
                #endif
            elif (i4 < 0):                        # 0110
                alpha1 = 0.
                if (T2 == T3):
                    alpha2 = 0.5
                else:
                    alpha2 = (T0 - T3) / (T2 - T3)
                #endif
                alpha3 = 1. - alpha2
                alpha4 = 0.
            else:                                # 0111
                # Assume that (T, g) is within the triangle i
                # formed by the three points.

                mat0 = numpy.asarray([ [ T2, T3, T4 ],
                                       [ g2, g3, g4 ],
                                       [ 1., 1.,  1.]  ])
                mat2 = numpy.asarray([ [ T0, T3, T4 ],
                                       [ g0, g3, g4 ],
                                       [ 1., 1.,  1.]  ])
                mat3 = numpy.asarray([ [ T2, T0, T4 ],
                                       [ g2, g0, g4 ],
                                       [ 1., 1.,  1.]  ])
                mat4 = numpy.asarray([ [ T2, T3, T0 ],
                                       [ g2, g3, g0 ],
                                       [ 1., 1.,  1.]  ])
                det0 = __det3x3__(mat0.ravel())
                det2 = __det3x3__(mat2.ravel())
                det3 = __det3x3__(mat3.ravel())
                det4 = __det3x3__(mat4.ravel())
                alpha1 = 0.
                alpha2 = det2 / det0
                alpha3 = det3 / det0
                alpha4 = det4 / det0

                # If (T, g) is outside the triangle formed
                # by the three used points use only two points.
                if ( (alpha2 < 0.) | (alpha2 > 1. ) | (alpha3 < 0.) | (alpha3 > 1.) | (alpha4 < 0.) | (alpha4 > 1. ) ):
                    alpha1 = 0.
                    if (T2 == T3):
                        alpha2 = 0.5
                    else:
                        alpha2 = (T0 - T3) / (T2 - T3)
                    #endif
                    alpha3 = 1. - alpha2
                    alpha4 = 0.
                    i4 = -1
                #endif
            #endif
        elif (i2 < 0):
            if (i3 < 0):
                if (i4 < 0):
                    #               #1000
                    alpha1 = 1.
                    alpha2 = 0.
                    alpha3 = 0.
                    alpha4 = 0.
                else:                      # 1001
                    if (T1 == T4):
                        alpha1 = 0.5
                    else:
                        alpha1 = (T0 - T4) / (T1 - T4)
                    # endif
                    alpha2 = 0.
                    alpha3 = 0.
                    alpha4 = 1. - alpha1
                #endif
            elif (i4 < 0):             # 1010
                if (T1 == T3):
                    alpha1 = 0.5
                else:
                    alpha1 = (T0 - T3) / (T1 - T3)
                #endif
                alpha2 = 0.
                alpha3 = 1. - alpha1
                alpha4 = 0.
            else:                     # 1011

                # Assume that (T, g) is within the triangle formed by the three points.

                mat0 = numpy.asarray([ [ T1, T3, T4 ],
                                       [ g1, g3, g4 ],
                                       [ 1., 1.,  1.]  ])
                mat1 = numpy.asarray([ [ T0, T3, T4 ],
                                       [ g0, g3, g4 ],
                                       [ 1., 1.,  1.]  ])
                mat3 = numpy.asarray([ [ T1, T0, T4 ],
                                       [ g1, g0, g4 ],
                                       [ 1., 1.,  1.]  ])
                mat4 = numpy.asarray([ [ T1, T3, T0 ],
                                       [ g1, g3, g0 ],
                                       [ 1., 1.,  1.]  ])
                det0 = __det3x3__(mat0.ravel())
                det1 = __det3x3__(mat1.ravel())
                det3 = __det3x3__(mat3.ravel())
                det4 = __det3x3__(mat4.ravel())
                alpha1 = det1 / det0
                alpha2 = 0.
                alpha3 = det3 / det0
                alpha4 = det4 / det0

                # If (T, g) is outside the triangle formed by the three used points,
                # use only two points.

                if ( (alpha1 < 0.) | (alpha1 > 1.) | (alpha3 < 0.) | (alpha3 > 1.) | (alpha4 < 0.) | (alpha4 > 1.) ):
                    if (T1 == T4):
                        alpha1 = 0.5
                    else:
                        alpha1 = (T0 - T4) / (T1 - T4)
                    #endif
                    alpha2 = 0.
                    alpha3 = 0.
                    i3 = -1
                    alpha4 = 1. - alpha1
                #endif
            #endif
        elif (i3 < 0):
            if (i4 < 0):
                #                   #1100
                if (abs(T1 - T2) < eps):
                    if (g1 == g2):
                        alpha1 = 0.5
                    else:
                        alpha1 = (g0 - g2) / (g1 - g2)
                    #endif
                    alpha2 = 1. - alpha1
                else:
                    if (T1 < T2):
                        alpha1 = 1.
                        alpha2 = 0.
                        i2 = -1
                    else:
                        alpha1 = 0.
                        i1 = -1
                        alpha2 = 1.
                    #endif
                #endif
                alpha3 = 0.
                alpha4 = 0.
            else:                       # 1101

                #Assume that (T, g) is within the triangle formed by the three points.

                mat0 = numpy.asarray([ [ T1, T2, T4 ],
                                       [ g1, g2, g4 ],
                                       [ 1., 1., 1. ]  ])
                mat1 = numpy.asarray([ [ T0, T2, T4 ],
                                       [ g0, g2, g4 ],
                                       [ 1., 1., 1. ]  ])
                mat2 = numpy.asarray([ [ T1, T0, T4 ],
                                       [ g1, g0, g4 ],
                                       [ 1., 1., 1. ]  ])
                mat4 = numpy.asarray([ [ T1, T2, T0 ],
                                       [ g1, g2, g0 ],
                                       [ 1., 1., 1. ]  ])
                det0 = __det3x3__(mat0.ravel())
                det1 = __det3x3__(mat1.ravel())
                det2 = __det3x3__(mat2.ravel())
                det4 = __det3x3__(mat4.ravel())
                alpha1 = det1 / det0
                alpha2 = det2 / det0
                alpha3 = 0.
                alpha4 = det4 / det0

                # If (T, g) is outside the triangle formed by the three used points,
                # use only two points.
                if ( (alpha1 < 0.) | (alpha1 > 1.) | (alpha2 < 0.) | (alpha2 > 1.) | (alpha4 < 0.) | (alpha4 > 1.) ):
                    if (T1 == T4):
                        alpha1 = 0.5
                    else:
                        alpha1 = (T0 - T4) / (T1 - T4)
                    #endif
                    alpha2 = 0.
                    i2 = -1
                    alpha3 = 0.
                    alpha4 = 1. - alpha1
                #endif
            #endif
        elif (i4 < 0):
            #                           # 1110
            #Assume that (T, g) is within the triangle formed by the three points.
            mat0 = numpy.asarray([ [ T1, T2, T3 ],
                                   [ g1, g2, g3 ],
                                   [ 1., 1., 1. ]  ])
            mat1 = numpy.asarray([ [ T0, T2, T3 ],
                                   [ g0, g2, g3 ],
                                   [ 1., 1., 1. ]  ])
            mat2 = numpy.asarray([ [ T1, T0, T3 ],
                                   [ g1, g0, g3 ],
                                   [ 1., 1., 1. ]  ])
            mat3 = numpy.asarray([ [ T1, T2, T0 ],
                                   [ g1, g2, g0 ],
                                   [ 1., 1., 1. ]  ])
            det0 = __det3x3__(mat0.ravel())
            det1 = __det3x3__(mat1.ravel())
            det2 = __det3x3__(mat2.ravel())
            det3 = __det3x3__(mat3.ravel())
            alpha1 = det1 / det0
            alpha2 = det2 / det0
            alpha3 = det3 / det0
            alpha4 = 0.

            # If (T, g) is outside the triangle formed by the three used points,
            # use only two points.
            if ( (alpha1 < 0.) | (alpha1 > 1.) | (alpha2 < 0.) | (alpha2 > 1.) | (alpha3 < 0.) | (alpha3 > 1.) ):
                alpha1 = 0.
                i1 = -1
                if (T2 == T3):
                    alpha2 = 0.5
                else:
                    alpha2 = (T0 - T3) / (T2 - T3)
                #endif
                alpha3 = 1. - alpha2
                alpha4 = 0.
            #endif
        #endif

        # All four points used.

        if ( (i3 >= 0) & (i4 >= 0) & (i1 >= 0) & (i2 >= 0) ):
            if (T1 != T3):
                alpha = (T0 - T3) / (T1 - T3)
            else:
                alpha = 0.5
            #endif
            if (T2 != T4):
                beta = (T0 - T4) / (T2 - T4)
            else:
                beta = 0.5
            #endif
            gprim = alpha * g1 + (1 - alpha) * g3
            gsec  = beta * g2  + (1 - beta ) * g4
            if (gprim != gsec):
                gamma = ( g0 - gsec ) / ( gprim - gsec )
            else:
                gamma = 0.5
            #endif
            alpha1 = alpha * gamma
            alpha2 = beta * ( 1 - gamma )
            alpha3 = ( 1 - alpha ) * gamma
            alpha4 = ( 1 - beta ) * ( 1 - gamma )
        #endif
        return numpy.asarray((i1, i2, i3, i4)), numpy.asarray((alpha1, alpha2, alpha3, alpha4))

    def __det3x3__(a):
        """ compute the 3x3 determinant of an array
            8 times faster than numpy.linalg.det for a matrix 3x3

        Inputs:
            a   3x3 array

        Returns the result as a float
        """
        # val  = +a[0,0] * ( a[1,1] * a[2,2] - a[2,1] * a[1,2] )
        # val += -a[0,1] * ( a[1,0] * a[2,2] - a[2,0] * a[1,2] )
        # val += +a[0,2] * ( a[1,0] * a[2,1] - a[2,0] * a[1,1] )
        val  = +a[0] * (a[4] * a[8] - a[7] * a[5])
        val += -a[1] * (a[3] * a[8] - a[6] * a[5])
        val += +a[2] * (a[3] * a[7] - a[6] * a[4])
        return val


def __interpSingle__(args):
    return numpy.asarray(interp(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9])).T


def __interpMany__(oSL, logT, logg, Z, logL, dT_max=0.1, eps=1e-06, weights=None, pool=None, nthreads=__NTHREADS__):
    if (pool is None) & (nthreads > 0):
        import multiprocessing as mp
        pool = mp.Pool(nthreads)
    if weights is None:
        seq = [ (logT[k], logg[k], Z, logL[k], oSL.Teff, oSL.logg, oSL.Z, dT_max, eps, 1.) for k in range(len(logT)) ]
    else:
        seq = [ (logT[k], logg[k], Z, logL[k], oSL.Teff, oSL.logg, oSL.Z, dT_max, eps, weights[k]) for k in range(len(logT)) ]
    if (pool is not None):
        r = pool.map( __interpSingle__, seq )
    else:
        r = map( __interpSingle__, seq )
    return numpy.vstack(r)


def interp(T0, g0, Z0, L0, T, g, Z, dT_max=0.1, eps=1e-6, weight=1.):
    """ Interpolation of the T,g grid

    Interpolate on the grid and returns star indices and
    associated weights, and Z.
    3 to 12 stars are returned.
    It calls _interp_, but reduce the output to the relevant stars.

    Inputs:
        T0  log(Teff) to obtain
        g0  log(g) to obtain
        Z0  metallicity to obtain
        L0  log(Lum) to obtain

    Keywords:
        dT_max  reject points with a temperature further than
            this value
        eps temperature sensitivity under which points are
            considered to have the same temperature

    Returns 3 to 12 star indexes and associated weights

    see _interp_

    TODO: compute new weights accounting for Z
    """
    _Z    = Z
    _Zv   = numpy.unique(_Z)
    _T    = numpy.asarray(T)
    _g    = numpy.asarray(g)
    bZ_m  = True in (_Zv == Z0)  # Z_match bool
    r     = numpy.where((_Zv < Z0))[0]
    Z_inf = _Zv[r.max()] if len(r) > 0 else -1.
    r     = numpy.where((_Zv > Z0))[0]
    Z_sup = _Zv[r.min()] if len(r) > 0 else -1.

    index   = numpy.zeros(4 * 3) - 1
    weights = numpy.zeros(4 * 3)
    Z       = numpy.zeros(4 * 3)

    if weight is None:
        weight = 1.

    if (bZ_m):
        ind         = numpy.where(_Z == Z0)
        i, w        = __interp__(T0, g0, _T[ind], _g[ind], dT_max, eps)
        index[8:]   = ind[0][i]
        weights[8:] = w
        Z[8:]       = [Z0] * 4
    else:
        if (Z_inf > 0.):
            ind         = numpy.where(_Z == Z_inf)
            i, w        = __interp__(T0, g0, _T[ind], _g[ind], dT_max, eps)
            index[:4]   = ind[0][i]
            weights[:4] = w
            Z[:4]       = [Z_inf] * 4

        if (Z_sup > 0.):
            ind          = numpy.where(_Z == Z_sup)
            i, w         = __interp__(T0, g0, _T[ind], _g[ind], dT_max, eps)
            index[4:8]   = ind[0][i]
            weights[4:8] = w
            Z[4:8]       = [Z_sup] * 4

        if ((Z_inf > 0.) & (Z_sup > 0.)):
            if ( Z_sup - Z_inf ) > 0.:
                fz = (Z0 - Z_inf) / ( Z_sup - Z_inf )
                weights[:4]  *= fz
                weights[4:8] *= ( 1. - fz )
            else:
                weights[:8]  *= 0.5

    ind = numpy.where(weights > 0)
    return index[ind].astype(int), 10 ** L0 * weight * weights[ind]  # / (weights[ind].sum()) #, Z[ind]


class Stellib(object):
    def __init__(self, *args, **kargs):
        """ Contructor """
        pass

    def _load_(self):
        pass

    def interp(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6):
        """ Interpolation of the T,g grid

        Interpolate on the grid and returns star indices and
        associated weights, and Z.
        3 to 12 stars are returned.
        It calls _interp_, but reduce the output to the relevant stars.

        Inputs:
            T0  log(Teff) to obtain
            g0  log(g) to obtain
            Z0  metallicity to obtain
            L0  log(Lum) to obtain

        Keywords:
            dT_max  reject points with a temperature further than
                this value
            eps temperature sensitivity under which points are
                considered to have the same temperature

        Returns 3 to 12 star indexes and associated weights

        see _interp_

        TODO: compute new weights accounting for Z
        """
        _Z    = self.Z
        #_Zv   = numpy.unique(_Z)
        _T    = numpy.asarray(self.grid['Teff'], dtype=numpy.double)
        _g    = numpy.asarray(self.grid['logG'], dtype=numpy.double)
        return interp(T0, g0, Z0, L0, _T, _g, _Z, dT_max=0.1, eps=1e-6)

    def interpMany(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6, weights=None, pool=None, nthreads=__NTHREADS__):
        """ run interp on a list of inputs and returns reduced
        results """

        #_t = numpy.asarray(T0)
        #_g = numpy.asarray(g0)
        r = __interpMany__(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-06, weights=weights, pool=pool, nthreads=nthreads)
        idx = numpy.unique(r[:, 0])
        # d = { idxk:0. for idxk in idx }
        d = {}
        for idxk in idx:
            d[idxk] = 0.
        for k in xrange(len(r)):
            d[ r[k, 0] ] += r[k, 1]
        del r, idx
        return numpy.asarray(d.items())

    def get_boundaries(self, dlogT=0.1, dlogg=0.3, closed=True):
        """ Returns the closed boundary polygon around the stellar library with
        given margins

        INPUTS:
            s   Stellib     Stellar library object

        KEYWORDS:
            dlogT   float       margin in logT
            dlogg   float       margin in logg
            closed  bool        if set, close the polygon

        OUTPUTS:
            b   ndarray[float, ndim=2]  (closed) boundary points: [logg, Teff]

        Note:
            use "points_inside_poly" to test wether a point is inside the limits
            >>> data = numpy.array([iso.data['logg'], iso.data['logT']]).T
            >>> aa = points_inside_poly(data, leftb)
        """
        leftb   = [(k, numpy.max(self.logT[self.logg == k]) + dlogT ) for k in numpy.unique(self.logg)]
        leftb  += [ (leftb[-1][0] + dlogg, leftb[-1][1]) ]
        leftb   = [ (leftb[0][0] - dlogg, leftb[0][1]) ] + leftb
        rightb  = [(k, numpy.min(self.logT[self.logg == k]) - dlogT ) for k in numpy.unique(self.logg)[::-1]]
        rightb += [ (rightb[-1][0] - dlogg, rightb[-1][1]) ]
        rightb  = [ (rightb[0][0] + dlogg, rightb[0][1]) ] + rightb
        b = leftb + rightb
        if closed:
            b += [b[0]]
        return numpy.array(b)

    def genQ(self, qname, r, **kwargs):
        """ Generate a composite value from a previously calculated
            interpolation
            Works on 1 desired star or a population of stars

        Inputs:
            qname   quantity name from self.grid
            r   the result from a previous interpolation

        Outputs:
            an array containing the value
        """
        return ( self.grid[qname][r[:, 0].astype(int)] * r[:, 1] ).sum()

    def genSpectrum(self, T0, g0=None, Z0=None, weights=None, **kwargs):
        """ Generate a composite sprectrum
            Does the interpolation or uses a previously calculated
            interpolation
            Works on 1 desired star or a population of stars

        Inputs:
            T0  log(Teff) of each star or a 2d-array containing
                the result from a previous interpolation
            g0  log(g) of each stars
            Z0  metallicity

            if T0 and g0 are iterable, it calls interpMany

        Keywords:
            weights individual weights of each star
            **kwargs forwarded to interp(Many)

        Outputs:
            an array containing the composite spectrum
        """
        if Z0 is not None:
            if hasattr(T0, '__iter__'):
                _r = self.interpMany(T0, g0, Z0, weights=weights, **kwargs)
            else:
                _r = numpy.asarray(self.interp(T0, g0, Z0, weight=weights, **kwargs))
        else:
            assert( T0.ndim == 2), 'error expecting 2d-array'
            _r = T0
        return ( ( (self.spectra[_r[:, 0].astype(int)].T) * _r[:, 1]) ).sum(1)


class Elodie(Stellib):
    """ Elodie 3.1 stellar library derived class """
    def __init__(self, *args, **kwargs):
        self.name = 'ELODIE v3.1 (Prugniel et al 2007, astro-ph/703658)'
        self.source = __ROOT__ + '/libs/stellib_ELODIE_3.1.fits'
        self._load_()

    def _load_(self):
        with  pyfits.open(self.source) as f:
            #load data
            self._getWaveLength_(f)
            self._getTGZ_(f)
            self._getSpectra_(f)

    def _getWaveLength_(self, f):
        self.wavelength = numpy.asarray((f[2].data[:]).tolist()).ravel()

    def _getTGZ_(self, f):
        cols = f[1].columns.names
        d = {}
        #d = { k: f[1].data.field(k) for k in cols }
        for k in cols:
            d[k] = f[1].data.field(k)
        self.grid = Table(d)
        self.grid.header['NAME'] = 'TGZ'
        del d, cols

    def _getSpectra_(self, f):
        self.spectra = f[0].data[:]

    @property
    def logg(self):
        return self.grid['logG']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']

    @property
    def NHI(self):
        return self.grid['NHI']

    @property
    def NHeI(self):
        return self.grid['NHeI']

    @property
    def NHeII(self):
        return self.grid['NHeII']


class BaSeL(Stellib):
    """ BaSeL 2.2 + Rauch stellar library derived class
        This library is used in Pegase.2
    """
    def __init__(self, *args, **kwargs):
        self.name = 'BaSeL 2.2 + Rauch (Pegase.2 version)'
        self.source = __ROOT__ + '/libs/stellib_BaSeL_2.2_Rauch.fits'
        self._load_()

    def _load_(self):
        with pyfits.open(self.source) as f:
            #load data
            self._getWaveLength_(f)
            self._getTGZ_(f)
            self._getSpectra_(f)

    def _getWaveLength_(self, f):
        self.wavelength = numpy.asarray((f[2].data[:]).tolist()).ravel()

    def _getTGZ_(self, f):
        cols = f[1].columns.names
        d = {}
        #d = { k: f[1].data.field(k) for k in cols }
        for k in cols:
            d[k] = f[1].data.field(k)
        self.grid = Table(d)
        self.grid.header['NAME'] = 'TGZ'
        del d, cols

    def _getSpectra_(self, f):
        self.spectra = f[0].data[:]

    @property
    def logg(self):
        return self.grid['logG']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']

    @property
    def NHI(self):
        return self.grid['NHI']

    @property
    def NHeI(self):
        return self.grid['NHeI']

    @property
    def NHeII(self):
        return self.grid['NHeII']


class Kurucz(Stellib):
    """
    The stellar atmosphere models by Castelli and Kurucz 2004
    """
    def __init__(self, *args, **kwargs):
        self.name = 'Kurucz 2004'
        self.source = __ROOT__ + '/libs/kurucz2004.grid.fits'
        self._load_()

    def _load_(self):
        g = grid.FileSpectralGrid(self.source)
        self.wavelength = g.lamb
        self.grid = g.grid
        self.grid.header['NAME'] = 'TGZ'
        self.spectra = g.seds

    @property
    def logT(self):
        return numpy.log10(self.grid['T0'])

    @property
    def logg(self):
        return self.grid['logG']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']

class Tlusty(Stellib):
    """
    Tlusty O and B stellar atmospheres
    """
    def __init__(self, *args, **kwargs):
        self.name = 'Tlusty'
        self.source = __ROOT__ + '/libs/tlusty.grid.fits'
        self._load_()

    def _load_(self):
        g = grid.FileSpectralGrid(self.source)
        self.wavelength = g.lamb
        self.grid = g.grid
        self.grid.header['NAME'] = 'tlusty'
        self.spectra = g.seds

    @property
    def logT(self):
        return numpy.log10(self.grid['Teff'])

    @property
    def logg(self):
        return self.grid['logG']

    @property
    def Teff(self):
        return self.grid['Teff']

    @property
    def Z(self):
        return self.grid['Z']


## utils for pegase files =======================================
#import numpy, pyfits
#
#def _readLine(f, nlines = 1, cols = None, convert=None, debug=None):
#   """
#   Read a given lines from a file stream
#   and optionaly convert field values into given format sequence
#   """
#   if nlines < 1:      return ()
#   if convert == None: convert = [str]
#   if cols == None:    cols = numpy.arange(numpy.size(convert))
#
#   data = {}
#   for k in cols: data[str(k)] = []
#   for il in range(nlines):
#       l = f.readline()
#       l = l.split()
#       for k in range(numpy.size(cols)):
#           data[str(cols[k])].append(convert[k](l[k]))
#   if nlines == 1:
#       return tuple([ data[str(k)][0] for k in cols ])
#   else:
#       return tuple([ data[str(k)] for k in cols ])
#
#
#def _readBlock(f, nvalues=0, convert=None, debug=False):
#   """
#   Read a given number of values in a file stream
#   optionaly it can also convert it into a given format
#   """
#   if nvalues == 0: return
#   data = []
#   s = len(data) # which here means 0!
#   while s < nvalues:
#       l = f.readline().split()
#       if debug: print l
#       data = numpy.hstack((data,l))
#       s = len(data)
#   if convert == None:
#       return numpy.array(data)
#   else:
#       return numpy.array([ convert(d) for d in data ])
#
#def read_pegase():
#   d = 'stel_lib_dir/'
#   name = 'stellib_BaSeL_2.2_Rauch.fits'
#
#   f_wave = open(d+'stellar_libraries_wavelengths.dat')
#   nw     = _readLine(f_wave, 1, convert = [int])[0]
#   wave   = _readBlock(f_wave, nw).astype(float)
#   f_wave.close()
#
#   f_root = open(d+'list_stellar_libraries.dat')
#   nw     = _readLine(f_root, 1, convert = [int, int])
#   l      = _readBlock(f_root, sum(nw))
#   f_root.close()
#
#   nw = len(wave)
#   def extract(fk):
#       f_fk = open(d+fk)
#       ns, z = _readLine(f_fk, 1, convert = [int, float])
#       #logT, logg, NHI, NHeI, NHeII
#       tgz   = _readBlock(f_fk, ns*5).astype(float)
#       tgz   = tgz.reshape((ns,5))
#       tgz   = numpy.vstack([ tgz.T , numpy.asarray([ z ]*ns) ]).T
#       s     = []
#       for k in range(ns):
#           s.append( _readBlock(f_fk, nw).astype(float) )
#       f_fk.close()
#
#       return tgz, numpy.asarray(s)
#
#   r = map(extract, l)
#   tgz   = numpy.vstack([ k[0] for k in r ])
#   specs = numpy.vstack([ k[1] for k in r ])
#
#   pars = ['Teff', 'logG', 'NHI', 'NHeI', 'NHeII', 'Z']
#   d    = { pars[k]:tgz[:,k] for k in range(len(pars)) }
#
#   t1    = Table(d)
#
#   t1.header['EXTNAME'] = 'TGZ'
#   t1.setUnit('Teff', 'log(K)')
#   t1.setUnit('logG', 'log(g/s**2)')
#   t1.header['COMMENT'] = 'Stellib. from Pegase.2x: BaSeL 2.2 + Rauch'
#
#   t2    = Table( {'BFIT': wave} )
#   t2.header['EXTNAME'] = 'WCA'
#   t2.setUnit('BFIT', 'AA')
#
#   pyfits.writeto(name, specs)
#   t1.write(name, append=True)
#   t2.write(name, append=True)
#
#
