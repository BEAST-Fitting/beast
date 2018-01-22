""" Some astro related functions """
from ..core.decorators import elementwise

import math
import numpy
from numpy import deg2rad, rad2deg, sin, cos, sqrt, arcsin, arctan2


PI = math.pi
HALFPI = PI / 2.0
D2R = PI / 180.0
R2D = 1.0 / D2R


@elementwise
def hms2deg(_str, delim=':'):
	""" Convert hex coordinates into degrees """
	if _str[0] == '-':
		neg = -1
		_str = _str[1:]
	else:
		neg = 1
	_str = _str.split(delim)
	return neg * ((((float(_str[-1]) / 60. + float(_str[1])) / 60. + float(_str[0])) / 24. * 360.))


@elementwise
def deg2dms(val, delim=':'):
    """ Convert degrees into hex coordinates """
    if val < 0:
        sign = -1
    else:
        sign = 1
    d = int( sign * val )
    m = int( (sign * val - d) * 60. )
    s = (( sign * val - d) * 60.  - m) * 60.
    return '{}{}{}{}{}'.format( sign * d, delim, m, delim, s)


@elementwise
def deg2hms(val, delim=':'):
    """ Convert degrees into hex coordinates """
    if val < 0:
        sign = -1
    else:
        sign = 1
    h = int( sign * val / 45. * 3.)   # * 24 / 360
    m = int( (sign * val / 45. * 3. - h) * 60. )
    s = (( sign * val / 45. * 3. - h) * 60.  - m) * 60.
    return '{}{}{}{}{}'.format( sign * h, delim, m, delim, s)


@elementwise
def dms2deg(_str, delim=':'):
	""" Convert hex coordinates into degrees """
	if _str[0] == '-':
		neg = -1
		_str = _str[1:]
	else:
		neg = 1
	_str = _str.split(delim)
	return (neg * ((float(_str[-1]) / 60. + float(_str[1])) / 60. + float(_str[0])))


def euler(ai_in, bi_in, select, b1950=False, dtype='f8'):
	"""
	Transform between Galactic, celestial, and ecliptic coordinates.

	INPUTS:
	long_in - Input Longitude in DEGREES, scalar or vector.
	lat_in  - Input Latitude in DEGREES
	select  - Integer (1-6) specifying type of coordinate transformation.

	select   From          To        |   select      From            To
	1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec
	2     Galactic       RA-DEC     |     5       Ecliptic      Galactic
	3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic

	Celestial coordinates (RA, Dec) should be given in equinox J2000
	unless the b1950=True keyword is set.

	OUTPUTS:
	long_out - Output Longitude in DEGREES
	lat_out  - Output Latitude in DEGREES

	KEYWORDS:
	b1950 - If this keyword is true then input and output
	     celestial and ecliptic coordinates should be given in equinox
	     B1950.

	REVISION HISTORY:
	Written W. Landsman,  February 1987
	Adapted from Fortran by Daryl Yentis NRL
	Converted to IDL V5.0   W. Landsman   September 1997
	Made J2000 the default, added /FK4 keyword  W. Landsman December 1998
	Add option to specify SELECT as a keyword W. Landsman March 2003

	Converted from IDL to numerical Python: Erin Sheldon, NYU, 2008-07-02

	"""

	# Make a copy as an array. ndmin=1 to avoid messed up scalar arrays
	ai = numpy.array(ai_in, ndmin=1, copy=True, dtype=dtype)
	bi = numpy.array(bi_in, ndmin=1, copy=True, dtype=dtype)

	twopi   = 2.0 * PI
	fourpi  = 4.0 * PI

	#   J2000 coordinate conversions are based on the following constants
	#   (see the Hipparcos explanatory supplement).
	#  eps = 23.4392911111d           Obliquity of the ecliptic
	#  alphaG = 192.85948d            Right Ascension of Galactic North Pole
	#  deltaG = 27.12825d             Declination of Galactic North Pole
	#  lomega = 32.93192d             Galactic longitude of celestial equator
	#  alphaE = 180.02322d            Ecliptic longitude of Galactic North Pole
	#  deltaE = 29.811438523d         Ecliptic latitude of Galactic North Pole
	#  Eomega  = 6.3839743d           Galactic longitude of ecliptic equator
	# Parameters for all the different conversions
	if b1950:
		#equinox = '(B1950)'
		psi    = numpy.array([ 0.57595865315, 4.9261918136,
				      0.00000000000, 0.0000000000,
				      0.11129056012, 4.7005372834], dtype=dtype)
		stheta = numpy.array([ 0.88781538514, -0.88781538514,
				      0.39788119938, -0.39788119938,
				      0.86766174755, -0.86766174755], dtype=dtype)
		ctheta = numpy.array([ 0.46019978478, 0.46019978478,
				      0.91743694670, 0.91743694670,
				      0.49715499774, 0.49715499774], dtype=dtype)
		phi    = numpy.array([ 4.9261918136,  0.57595865315,
				      0.0000000000, 0.00000000000,
				      4.7005372834, 0.11129056012], dtype=dtype)
	else:
		#equinox = '(J2000)'
		psi    = numpy.array([ 0.57477043300, 4.9368292465,
				      0.00000000000, 0.0000000000,
				      0.11142137093, 4.71279419371], dtype=dtype)
		stheta = numpy.array([ 0.88998808748, -0.88998808748,
				      0.39777715593, -0.39777715593,
				      0.86766622025, -0.86766622025], dtype=dtype)
		ctheta = numpy.array([ 0.45598377618, 0.45598377618,
				      0.91748206207, 0.91748206207,
				      0.49714719172, 0.49714719172], dtype=dtype)
		phi    = numpy.array([ 4.9368292465,  0.57477043300,
				      0.0000000000, 0.00000000000,
				      4.71279419371, 0.11142137093], dtype=dtype)

	# zero offset
	i  = select - 1
	a  = ai * D2R - phi[i]

	b = bi * D2R
	sb = sin(b)
	cb = cos(b)
	cbsa = cb * sin(a)
	b  = -stheta[i] * cbsa + ctheta[i] * sb
	w, = numpy.where(b > 1.0)
	if w.size > 0:
		b[w] = 1.0
	bo = arcsin(b) * R2D
	a  = arctan2( ctheta[i] * cbsa + stheta[i] * sb, cb * cos(a) )
	ao = ( (a + psi[i] + fourpi) % twopi) * R2D
	return ao, bo


def sphdist(ra1, dec1, ra2, dec2):
	"""measures the spherical distance between 2 points
	Inputs:
		(ra1, dec1)	in degrees
		(ra2, dec2)	in degrees
	Outputs:
		returns a distance in degrees
	"""
	dec1_r = deg2rad(dec1)
	dec2_r = deg2rad(dec2)
	return 2. * rad2deg( arcsin( sqrt( ( sin((dec1_r - dec2_r) / 2)) ** 2 + cos(dec1_r) * cos(dec2_r) * ( sin((deg2rad(ra1 - ra2)) / 2)) ** 2)))


def conesearch(ra0, dec0, ra, dec, r, outtype=0):
    """ Perform a cone search on a table
    INPUTS:
        ra0 	ndarray[ndim=1, dtype=float]	column name to use as RA source in degrees
        dec0	ndarray[ndim=1, dtype=float]	column name to use as DEC source in degrees
        ra		float                       	ra to look for (in degree)
        dec		float	                        ra to look for (in degree)
        r		float		                    distance in degrees
    KEYWORDS:
        outtype int                             0 -- minimal, indices of matching coordinates
                                                1 -- indices and distances of matching coordinates
                                                2 -- full, boolean filter and distances

    """
    @elementwise
    def getDist( pk ):
        """ get spherical distance between 2 points """
        return sphdist(pk[0], pk[1], ra, dec)

    dist = numpy.array(getDist(list(zip(ra0, dec0))))
    v = (dist <= r)

    if outtype == 0:
        return numpy.ravel(numpy.where(v))
    elif outtype == 1:
        return numpy.ravel(numpy.where(v)), dist[v]
    else:
        return v, dist


