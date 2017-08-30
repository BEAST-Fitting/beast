"""
Contains a class for managing paths (polylines).
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from matplotlib.path import Path
from matplotlib.nxutils import points_inside_poly


class Path(Path):
    """
    :class:`Path` represents a series of possibly disconnected,
    possibly closed, line and curve segments.

    The underlying storage is made up of two parallel numpy arrays:
      - *vertices*: an Nx2 float array of vertices
      - *codes*: an N-length uint8 array of vertex types

    These two arrays always have the same length in the first
    dimension.  For example, to represent a cubic curve, you must
    provide three vertices as well as three codes ``CURVE3``.

    The code types are:

       - ``STOP``   :  1 vertex (ignored)
           A marker for the end of the entire path (currently not
           required and ignored)

       - ``MOVETO`` :  1 vertex
            Pick up the pen and move to the given vertex.

       - ``LINETO`` :  1 vertex
            Draw a line from the current position to the given vertex.

       - ``CURVE3`` :  1 control point, 1 endpoint
          Draw a quadratic Bezier curve from the current position,
          with the given control point, to the given end point.

       - ``CURVE4`` :  2 control points, 1 endpoint
          Draw a cubic Bezier curve from the current position, with
          the given control points, to the given end point.

       - ``CLOSEPOLY`` : 1 vertex (ignored)
          Draw a line segment to the start point of the current
          polyline.

    Users of Path objects should not access the vertices and codes
    arrays directly.  Instead, they should use :meth:`iter_segments`
    to get the vertex/code pairs.  This is important, since many
    :class:`Path` objects, as an optimization, do not store a *codes*
    at all, but have a default one provided for them by
    :meth:`iter_segments`.

    .. note::

        The vertices and codes arrays should be treated as
        immutable -- there are a number of optimizations and assumptions
        made up front in the constructor that will not change when the
        data changes.

    """
    def contains_points(self, points, transform=None, **kwargs):
        """
        Returns a bool array which is *True* if the path contains the
        corresponding point.

        If *transform* is not *None*, the path will be transformed
        before performing the test.

        *radius* allows the path to be made slightly larger or
        smaller.
        """
        poly = self.to_polygons()
        gen = (points_inside_poly(points, pk) for pk in poly)
        result = np.sum(gen, axis=0)
        return result.astype(bool)
