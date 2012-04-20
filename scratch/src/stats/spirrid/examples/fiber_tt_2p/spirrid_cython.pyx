import numpy as np
cimport numpy as np
ctypedef np.double_t DTYPE_t
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def mu_q(np.ndarray[DTYPE_t, ndim=1] la_flat, double e_arr, np.ndarray[DTYPE_t, ndim=1] xi_flat):
    cdef double mu_q
    cdef double la, xi, eps = e_arr, dG, q
    cdef int i_la, i_xi
    mu_q = 0
    dG = 6.25e-06
    for i_la from 0 <= i_la <400:
        la = la_flat[ i_la ]
        for i_xi from 0 <= i_xi <400:
            xi = xi_flat[ i_xi ]
            
            # Computation of the q( ... ) function
            if eps < 0 or eps > xi:
                q = 0.0
            else:
                q = la * eps
            
                mu_q += q * dG


    return mu_q