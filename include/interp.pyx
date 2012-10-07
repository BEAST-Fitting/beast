from __future__ import division
import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float_t DTYPE_t

cdef extern from "numpy/npy_math.h":
	bint npy_isnan(DTYPE_t x)

cimport cython


@cython.boundscheck(False)
cdef DTYPE_t __det3x3__1d(np.ndarray[DTYPE_t, ndim=1] a):
	"""
	 Description: compute the 3x3 determinant of an array 
	 This is 8 times faster than numpy.linalg.det for a matrix 3x3 when in
	 python

	 Inputs:
	     np.ndarray[DTYPE_t, ndim=1] a[9]
	 Returns:
	     double val
	"""
	cdef DTYPE_t val = 0.0
	val  = +a[0]*(a[4]*a[8]-a[7]*a[5])
	val += -a[1]*(a[3]*a[8]-a[6]*a[5])
	val += +a[2]*(a[3]*a[7]-a[6]*a[4])
	return val

@cython.boundscheck(False)
cdef DTYPE_t __det3x3__(np.ndarray[DTYPE_t, ndim=2] a):
	"""
	 Name: f_det3x3
	 Description: compute the 3x3 determinant of an array 
	 This is 8 times faster than numpy.linalg.det for a matrix 3x3 when in
	 python

	 Inputs:
	     double* a (numpy.ndarray[3,3])
	 Returns:
	     double val
	"""
	cdef DTYPE_t val = 0.0
	val  = +a[0,0]*(a[1,1]*a[2,2]-a[1,2]*a[2,1])
	val += -a[1,0]*(a[0,1]*a[2,2]-a[0,2]*a[2,1])
	val += +a[2,0]*(a[0,1]*a[1,2]-a[0,2]*a[1,1])
	return val


@cython.boundscheck(False)
def __interp__(DTYPE_t T0, DTYPE_t g0, 
		np.ndarray[cython.long, ndim=1] tidx,
		np.ndarray[DTYPE_t, ndim=1] talpha,
		int n_spec,
		np.ndarray[DTYPE_t, ndim=1] T,
		np.ndarray[DTYPE_t, ndim=1] g,
		DTYPE_t dT_max=0.1, 
		DTYPE_t eps=1e-6,
		DTYPE_t kappa = 0.1 ):
	"""
	Interpolation of the (T,g) grid at fixed Z

	Translated from Pegase.2 fortran version
	(this may not be pythonic though)

	Note: preference is always given to the temperature over 
		the gravity when needed.

	Inputs:
		T0	log(Teff) to obtain
		g0	log(g) to obtain
		T	log(Teff) of the grid 
		g	log(g) of the grid 

	Keywords:
		dT_max	If, T2 (resp. T1) is too far from T compared to T1
			(resp. T2), i2 (resp. i1) is not used. 
			(see below for namings)	
		eps	temperature sensitivity under which points are
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
	#******* declarations *****
	cdef int i
	cdef long i1, i2, i3, i4
	cdef DTYPE_t alpha, beta, gamma, gprim, gsec, g1, g2, g3, g4
	cdef DTYPE_t alpha1, alpha2, alpha3, alpha4
	cdef DTYPE_t T1, T2, T3, T4
	cdef DTYPE_t delta1, delta2, delta3, delta4
	cdef DTYPE_t det0, det1, det2, det3, det4
	cdef np.ndarray[DTYPE_t, ndim=1]  mat0, mat1, mat2, mat3, mat4

	cdef DTYPE_t dist


	#******* init ******
	i1 = -1
	i2 = -1
	i3 = -1
	i4 = -1
	alpha1 = 0.
	alpha2 = 0.
	alpha3 = 0.
	alpha4 = 0.
	delta1 = 1e4
	delta2 = delta1
	delta3 = delta1
	delta4 = delta1

	for i in range(n_spec):
		dist = kappa * abs(g[i]-g0) + abs(T[i]-T0)
		if (T[i] >= T0):
			if (g[i] >= g0):
				if ( dist < delta1 ):
					i1 = i
					delta1 = dist
			else:
				if ( dist < delta2 ):
					i2 = i
					delta2 = dist


		else: #(T_spec[i] < T)
			if (g[i] >= g0):
				if (dist < delta3):
					i3 = i
					delta3 = dist
			else:
				if ( dist < delta4):
					i4 = i
					delta4 = dist
     

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
	#	the gravity when needed.

	if (i1 < 0):
		if (i2 < 0):
			if (i3 < 0):
				if (i4 < 0):
					#               #0000
					assert (False), "Error"
				else:                   #0001
					alpha1 = 0.
					alpha2 = 0.
					alpha3 = 0.
					alpha4 = 1.
				#endif
			elif (i4 < 0):                  #0010
				alpha1 = 0.
				alpha2 = 0.
				alpha3 = 1.
				alpha4 = 0.
			else:                           #0011
				alpha1 = 0.
				alpha2 = 0.
				if ( abs(T3-T4) < eps ):
					if (g3 == g4):
						alpha3 = 0.5
					else: 
						alpha3 = (g0-g4)/(g3-g4)
					#endif
					alpha4 = 1.-alpha3
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
			else:                         #0101
				alpha1 = 0.
				if (T2 == T4):
					alpha2 = 0.5
				else:
					alpha2 = (T0-T4)/(T2-T4)
				#endif
				alpha3 = 0.
				alpha4 = 1.-alpha2
			#endif
		elif (i4 < 0):                        #0110
			alpha1 = 0.
			if (T2 == T3):
				alpha2 = 0.5
			else:
				alpha2 = (T0-T3)/(T2-T3)
			#endif
			alpha3 = 1.-alpha2
			alpha4 = 0.
		else:                                #0111
			# Assume that (T, g) is within the triangle i
			# formed by the three points.

			mat0 = np.asarray([ 
			      T2, T3, T4 , 
			      g2, g3, g4 , 
			      1.,    1.,     1.  ])
			mat2 = np.asarray([ 
			      T0   , T3, T4 , 
			      g0   , g3, g4 , 
			         1.,    1.,     1.  ])
			mat3 = np.asarray([ 
			      T2, T0   , T4 , 
			      g2, g0   , g4 , 
			         1.,    1.,     1.  ])
			mat4 = np.asarray([ 
			      T2, T3, T0    , 
			      g2, g3, g0    , 
			         1.,    1.,     1.  ])
			det0 = __det3x3__1d(mat0)
			det2 = __det3x3__1d(mat2)
			det3 = __det3x3__1d(mat3)
			det4 = __det3x3__1d(mat4)
			alpha1 = 0.
			alpha2 = det2/det0
			alpha3 = det3/det0
			alpha4 = det4/det0

			# If (T, g) is outside the triangle formed 
			# by the three used points use only two points.
			if ( (alpha2 < 0.) | (alpha2 > 1. ) | (alpha3 < 0.) | (alpha3 > 1.) | (alpha4 < 0.) | (alpha4 > 1. ) ):
				alpha1 = 0.
				if (T2 == T3):
					alpha2 = 0.5
				else:
					alpha2 = (T0-T3)/(T2-T3)
				#endif
				alpha3 = 1.-alpha2
				alpha4 = 0.
				i4 = -1
			#endif
		#endif
	elif (i2 < 0): 
		if (i3 < 0):
			if (i4 < 0):
				#	            #1000
				alpha1 = 1.
				alpha2 = 0.
				alpha3 = 0.
				alpha4 = 0.
			else:                      #1001
				if (T1 == T4):
					alpha1 = 0.5
				else: 
					alpha1 = (T0-T4)/(T1-T4)
				#endif
				alpha2 = 0.
				alpha3 = 0.
				alpha4 = 1.-alpha1
			#endif
		elif (i4 < 0):			   #1010
			if (T1 == T3):
				alpha1 = 0.5
			else:
				alpha1 = (T0-T3)/(T1-T3)
			#endif
			alpha2 = 0.
			alpha3 = 1.-alpha1
			alpha4 = 0.
		else: 					  #1011

			# Assume that (T, g) is within the triangle formed by the three points.

			mat0 = np.asarray([ 
			      T1, T3, T4 , 
			      g1, g3, g4 , 
			         1.,    1.,     1.  ])
			mat1 = np.asarray([ 
			      T0   , T3, T4 , 
			      g0   , g3, g4 , 
			         1.,    1.,    1. ])
			mat3 = np.asarray([ 
			      T1, T0   , T4 , 
			      g1, g0   , g4 , 
			         1.,    1.,    1.  ])
			mat4 = np.asarray([ 
			      T1, T3, T0    , 
			      g1, g3, g0    , 
			         1.,    1.,     1.  ])
			det0 = __det3x3__1d(mat0)
			det1 = __det3x3__1d(mat1)
			det3 = __det3x3__1d(mat3)
			det4 = __det3x3__1d(mat4)
			alpha1 = det1/det0
			alpha2 = 0.
			alpha3 = det3/det0
			alpha4 = det4/det0

			# If (T, g) is outside the triangle formed by the three used points,
			# use only two points.

			if ( (alpha1 < 0.) | (alpha1 > 1.) | (alpha3 < 0.) | (alpha3 > 1.) | (alpha4 < 0.) | (alpha4 > 1.) ):
				if (T1 == T4):
					alpha1 = 0.5
				else:
					alpha1 = (T0-T4)/(T1-T4)
				#endif
				alpha2 = 0.
				alpha3 = 0. 
				i3 = -1
				alpha4 = 1.-alpha1
			#endif
		#endif
	elif (i3 < 0): 
		if (i4 < 0):
			#					#1100
			if (abs(T1-T2) < eps):
				if (g1 == g2):
					alpha1 = 0.5
				else:
					alpha1 = (g0-g2)/(g1-g2)
				#endif
				alpha2 = 1.-alpha1
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
		else:						#1101

			#Assume that (T, g) is within the triangle formed by the three points.

			mat0 = np.asarray([ 
			      T1, T2, T4 , 
			      g1, g2, g4 , 
			         1.,    1.,     1.  ])
			mat1 = np.asarray([ 
			      T0   , T2, T4 , 
			      g0   , g2, g4 , 
			         1.,    1.,    1.  ])
			mat2 = np.asarray([ 
			      T1, T0   , T4 , 
			      g1, g0   , g4 , 
			         1.,    1.,    1.  ])
			mat4 = np.asarray([ 
			      T1, T2, T0    , 
			      g1, g2, g0    , 
			         1.,    1.,     1.  ])
			det0 = __det3x3__1d(mat0)
			det1 = __det3x3__1d(mat1)
			det2 = __det3x3__1d(mat2)
			det4 = __det3x3__1d(mat4)
			alpha1 = det1/det0
			alpha2 = det2/det0
			alpha3 = 0.
			alpha4 = det4/det0


			# If (T, g) is outside the triangle formed by the three used points,
			# use only two points.

			if ( (alpha1 < 0.) | (alpha1 > 1.) | (alpha2 < 0.) | (alpha2 > 1.) | (alpha4 < 0.) | (alpha4 > 1.) ):
				if (T1 == T4):
					alpha1 = 0.5
				else:
					alpha1 = (T0-T4)/(T1-T4)
				#endif
				alpha2 = 0.
				i2 = -1
				alpha3 = 0.
				alpha4 = 1.-alpha1
			#endif
		#endif
	elif (i4 < 0):
		#							#1110
		#Assume that (T, g) is within the triangle formed by the three points.
		mat0 = np.asarray([ 
		      T1, T2, T3 , 
		      g1, g2, g3 , 
		         1.,    1.,     1.  ])
		mat1 = np.asarray([ 
		      T0   , T2, T3 , 
		      g0   , g2, g3 , 
		         1.,    1.,     1.  ])
		mat2 = np.asarray([ 
		      T1, T0   , T3 , 
		      g1, g0   , g3 , 
		         1.,    1.,     1.  ])
		mat3 = np.asarray([ 
		      T1, T2, T0    , 
		      g1, g2, g0    , 
		         1.,    1.,     1.  ])
		det0 = __det3x3__1d(mat0)
		det1 = __det3x3__1d(mat1)
		det2 = __det3x3__1d(mat2)
		det3 = __det3x3__1d(mat3)
		alpha1 = det1/det0
		alpha2 = det2/det0
		alpha3 = det3/det0
		alpha4 = 0.


		# If (T, g) is outside the triangle formed by the three used points,
		# use only two points.

		if ( (alpha1 < 0.) | (alpha1 > 1.) | (alpha2 < 0.) | (alpha2 > 1.) | (alpha3 < 0.) | (alpha3 > 1.) ):
			alpha1 = 0.
			i1 = -1
			if (T2 == T3):
				alpha2 = 0.5
			else:
				alpha2 = (T0-T3)/(T2-T3)
			#endif
			alpha3 = 1.-alpha2
			alpha4 = 0.
		#endif
	#endif

	# All four points used.

	if ( (i3 > 0) & (i4 > 0) & (i1 > 0) & (i2 > 0) ):
		if (T1 != T3):
			alpha = (T0-T3)/(T1-T3)
		else:
			alpha = 0.5
		#endif
		if (T2 != T4):
			beta = (T0-T4)/(T2-T4)
		else:
			beta = 0.5
		#endif
		gprim = alpha * g1 + (1 - alpha)*g3
		gsec  = beta * g2  + (1 - beta )*g4
		if (gprim != gsec):
			gamma = ( g0 - gsec ) / ( gprim - gsec )
		else:
			gamma = 0.5
		#endif
		alpha1 =       alpha * gamma
		alpha2 =        beta * ( 1-gamma )
		alpha3 = ( 1-alpha ) * gamma
		alpha4 = (  1-beta ) * ( 1-gamma )
	#endif
	tidx[0] = i1
	tidx[1] = i2
	tidx[2] = i2
	tidx[3] = i3

	talpha[0] = alpha1
	talpha[1] = alpha2
	talpha[2] = alpha2
	talpha[3] = alpha3
	#return np.asarray((i1,i2,i3,i4)), np.asarray((alpha1, alpha2, alpha3, alpha4))

