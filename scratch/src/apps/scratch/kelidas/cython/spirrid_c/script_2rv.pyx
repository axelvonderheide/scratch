#from scipy.stats.distributions import norm, weibull_min, uniform # import normal distribution
#import pylab as p  # import plotting tool
#from numpy import vectorize, linspace, zeros_like, sign, sum as nsum, ones, \
#                 corrcoef, sort, diff, array, loadtxt
#from numpy.random import random
#from time import clock
#from scipy.interpolate import interp1d

import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.double
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.double_t DTYPE_t
# "def" can type its arguments but not have a return type. The type of the
# arguments for a "def" function is checked at run-time when entering the
# function.
#
# The arrays f, g and h is typed as "np.ndarray" instances. The only effect
# this has is to a) insert checks that the function arguments really are
# NumPy arrays, and b) make some attribute access like f.shape[0] much
# more efficient. (In this example this doesn't matter though.)
cimport cython


cdef int Heaviside( double x ):
    ''' Heaviside function '''
    return ( np.sign( x ) + 1.0 ) / 2.0

cdef double q( double e, double la, double xi ):
    ''' Response function of a single fiber '''
    return  e / ( 1 + la ) * Heaviside( xi - e / ( 1 + la ) )

#def q_ex( e ):
#    data = loadtxt( '2rv_maple.txt', delimiter = ',' )
#    x, y = data[:, 0], data[:, 1]
#    f = interp1d( x, y, kind = 'linear' )
#    return f( e )



@cython.boundscheck(False) # turn of bounds-checking for entire function
def loop_mu_q_e(np.ndarray[DTYPE_t, ndim=1] e_arr, double d_la, double d_xi,np.ndarray[DTYPE_t, ndim=1] Theta_la, np.ndarray[DTYPE_t, ndim=1] Theta_xi, np.ndarray[DTYPE_t, ndim=1] g_la_pdf, np.ndarray[DTYPE_t, ndim=1] g_xi_pdf):
    # loop over the control variable (strain)
    # define an array of the same size as e_arr
    cdef np.ndarray mu_q_arr = np.zeros_like( e_arr, dtype=DTYPE )
    cdef int i, la, xi
    cdef double mu_q_e, dG, e, eps_, q
    for i in range(len(e_arr)):
        mu_q_e = 0.0 
        e = e_arr[i]    
        for la in range(len(Theta_la)): # loop over lambda range (array of values)
            for xi in range(len(Theta_xi)): # loop over xi range (array of values)
                dG = g_la_pdf[la] * g_xi_pdf[xi] * d_la * d_xi
                eps_ = Theta_xi[xi] - e / ( 1.0 + Theta_la[la] )
                if eps_ < 0.0 or eps_ > Theta_xi[xi]:
                    q = 0.0
                else:
                    q = e / ( 1.0 + Theta_la[la] )
                mu_q_e += q * dG
        mu_q_arr[ i ] = mu_q_e
    return mu_q_arr





#dG_la = g_la.pdf( Theta_la ) * d_la
#dG_xi = g_xi.pdf( Theta_xi ) * d_xi

#dG_grid = dG_la[:, None] * dG_xi[None, :]

#def mu_q_e( e ):
    #''' Summation / integration  over the random domain '''
    #q_e_grid = q( e, Theta_la[:, None], Theta_xi[None, :] )
    #q_dG_grid = q_e_grid * dG_grid # element by element product of two (m,m) arrays
    #return nsum( q_dG_grid ) # nsum has been imported at line 3 from numpy 

#mu_q_e_vct = vectorize( mu_q_e )

#start_time = clock()
#mu_q_arr = mu_q_e_vct( e_arr )
#print 'loop-less: elapsed time', clock() - start_time

#def mu_q_e_LHS( e ):
    #''' Summation / integration  over the random domain '''
    #q_e_grid = q( e, T_la[:, None], T_xi[None, :] )
    #return nsum( q_e_grid ) / n_int ** n_k

#mu_q_e_LHS_vct = vectorize( mu_q_e_LHS )

#start_time = clock()
#mu_q_arr_LHS = mu_q_e_LHS_vct( e_arr )
#print 'loop-less: elapsed time', clock() - start_time


#def mu_q_e_MC( e ):
    #''' Summation / integration  over the random domain '''
    #q_e_grid = q( e, T_la_MC, T_xi_MC )
    #return nsum( q_e_grid ) / n_int ** n_k

#mu_q_e_MC_vct = vectorize( mu_q_e_MC )

#start_time = clock()
#mu_q_arr_MC = mu_q_e_MC_vct( e_arr )
#print 'loop-less: elapsed time', clock() - start_time

#mu_q_ex = q_ex( e_arr )
######################
## ERROR
########################xx
#print 'reg grid error', nsum( ( mu_q_arr - mu_q_ex ) ** 2 )
#print 'LHS grid error', nsum( ( mu_q_arr_LHS - mu_q_ex ) ** 2 )
#print 'MC error', nsum( ( mu_q_arr_MC - mu_q_ex ) ** 2 )

#p.plot( e_arr, mu_q_ex, 'k-', label = 'maple' )
#p.plot( e_arr, mu_q_arr, 'b-', label = 'regular grid' )
#p.plot( e_arr, mu_q_arr_LHS, 'r-', label = 'LHS grid' )
#p.plot( e_arr, mu_q_arr_MC, 'g-', label = 'MC' )
#p.plot( Theta_xi, 0.003 * ones( n_int ), 'bo' )
#p.plot( Theta_xi * ( 1 + Theta_la ), .0028 * ones( n_int ), 'b|', label = '$\\xi(1+\lambda)$' )
#p.plot( Theta_xi[None, :] * ( 1 + Theta_la[:, None] ), 0.0028 * ones( n_int ), 'b|' )
#p.plot( T_xi, 0.002 * ones( n_int ), 'ro' )
#p.plot( T_xi * ( 1 + T_la ), 0.0018 * ones( n_int ), 'r|', label = '$\\xi(1+\lambda)$' )
#p.plot( T_xi[None, :] * ( 1 + T_la[:, None] ), 0.0018 * ones( n_int ), 'r|' )
#p.plot( T_xi_MC, 0.001 * ones( n_int ** 2 ), 'go' )
#p.plot( 0, 0, 'g|', label = '$\\xi(1+\lambda)$' )
#p.plot( T_xi_MC * ( 1 + T_la_MC ), 0.0008 * ones( n_int ** 2 ), 'g|' )
#p.legend()

#p.figure()
#p.plot( T_la[None, :] * ones( n_int )[:, None], T_xi[:, None] * ones( n_int )[None, :], 'ro', label = 'LHS grid' )
#p.plot( T_la_MC, T_xi_MC, 'gx', label = 'MC' )
#p.plot( Theta_la[None, :] * ones( n_int )[:, None], Theta_xi[:, None] * ones( n_int )[None, :], 'b.', label = 'regular grid' )


#p.show()
