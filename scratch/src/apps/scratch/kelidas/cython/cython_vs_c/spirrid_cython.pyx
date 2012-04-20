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

@cython.cdivision(True)
@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.wraparound(False)
def loop_mu_q_e(np.ndarray[DTYPE_t, ndim=1] e_arr, double d_la, double d_xi,np.ndarray[DTYPE_t, ndim=1] Theta_la, np.ndarray[DTYPE_t, ndim=1] Theta_xi, np.ndarray[DTYPE_t, ndim=1] g_la_pdf, np.ndarray[DTYPE_t, ndim=1] g_xi_pdf):
    # loop over the control variable (strain)
    # define an array of the same size as e_arr
    cdef np.ndarray mu_q_arr = np.zeros_like( e_arr, dtype=DTYPE )
    cdef int i, la, xi
    cdef double mu_q_e, dG, e, eps_, q
    for i from 0 <= i < 1000:
        mu_q_e = 0.0 
        for la from 0 <= la < 1000: # loop over lambda range (array of values)
            for xi from 0 <= xi < 1000: # loop over xi range (array of values)
                e = e_arr[i]
                q = e / ( 1.0 + Theta_la[la] )
                eps_ = Theta_xi[xi] - q
                if eps_ < 0.0 or eps_ > Theta_xi[xi]:
                    q = 0.0
                mu_q_e += q * g_la_pdf[la] * g_xi_pdf[xi] * d_la * d_xi
        mu_q_arr[ i ] = mu_q_e
    return mu_q_arr
    
    