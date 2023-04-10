cimport cython
import numpy as np
cimport numpy as np

from libc.math cimport abs, pow, log

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)

cpdef inline compute_KL_PRURm_enrichment(
    np.ndarray[np.float_t, ndim=1] S_h,
    np.ndarray[np.int_t, ndim=1] I_h,
    double alpha):

    cdef int i, n = len(S_h)
    cdef np.float_t ES, p_sum = 1.0e-16, r_sum = 1.0e-16, u_sum = 1.0e-16, Delta_sum = 0.0
    cdef np.ndarray[np.float_t, ndim = 1] p_cumsum = np.empty(n, dtype = float)
    cdef np.ndarray[np.float_t, ndim = 1] r_cumsum = np.empty(n, dtype = float)
    cdef np.ndarray[np.float_t, ndim = 1] u_cumsum = np.empty(n, dtype = float)
    cdef np.ndarray[np.float_t, ndim = 1] P = np.empty(n, dtype = float)
    cdef np.ndarray[np.float_t, ndim = 1] U = np.empty(n, dtype = float)
    cdef np.ndarray[np.float_t, ndim = 1] R = np.empty(n, dtype = float)
    cdef np.ndarray[np.float_t, ndim = 1] Delta = np.empty(n, dtype = float)
    
    if alpha == 1.0:
        
        for i in range(n):
            if I_h[i] == 1:
                p_sum += abs(S_h[i]) 
            else:
                u_sum += abs(S_h[i]) 

            r_sum += abs(S_h[i]) 
            p_cumsum[i] = p_sum   
            u_cumsum[i] = u_sum
            r_cumsum[i] = r_sum

    else:
        
        for i in range(n):
            if I_h[i] == 1:
                p_sum += pow(abs(S_h[i]), alpha) 
            else:
                u_sum += pow(abs(S_h[i]), alpha) 

            r_sum += pow(abs(S_h[i]), alpha) 
            p_cumsum[i] = p_sum
            u_cumsum[i] = u_sum
            r_cumsum[i] = r_sum
        
    for i in range(n):
        P[i] = p_cumsum[i]/p_sum
        U[i] = u_cumsum[i]/u_sum
        R[i] = r_cumsum[i]/r_sum
        Delta[i] = (P[i] * log(P[i]/R[i])) - (U[i] * log(U[i]/R[i])) - \
                   ((1. - P[i] + 1.0e-16) * log((1. - P[i] + 1.0e-16)/(1. - R[i] + 1.0e-16)) \
                    - (1. - U[i] + 1.0e-16) * log((1. - U[i] + 1.0e-16)/(1. - R[i] + 1.0e-16))) 
         
        Delta_sum += Delta[i]

    ES = 0.5*Delta_sum/n
    
    return Delta, ES