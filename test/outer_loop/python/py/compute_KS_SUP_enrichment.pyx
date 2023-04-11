cimport cython
import numpy as np
cimport numpy as np

from libc.math cimport abs, pow

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)

cpdef inline compute_KS_SUP_enrichment(
    np.ndarray[np.float_t, ndim=1] S_h,
    np.ndarray[np.int_t, ndim=1] I_h,
    double alpha):
    #np.float_t [:] S_h, # Sorted and normalized gene scores
    #np.int_t [:] I_h, # gene set membership indicator (1 = gene in set, 0 = gene not in set)
    #np.float_t alpha): # exponent

    cdef int i, n = len(S_h)
    cdef np.float_t ES, p_sum = 1.0e-16, q_sum = 1.0e-16
    cdef np.ndarray[np.float_t, ndim = 1] p_cumsum = np.empty(n, dtype = float)
    cdef np.ndarray[np.float_t, ndim = 1] q_cumsum = np.empty(n, dtype = float)
    cdef np.ndarray[np.float_t, ndim = 1] Delta = np.empty(n, dtype = float)
    
    for i in range(n):
        if I_h[i] == 1:
            p_sum += pow(abs(S_h[i]), alpha)
        else:
            q_sum += 1 

        p_cumsum[i] = p_sum
        q_cumsum[i] = q_sum
            
    for i in range(n):
        Delta[i] = p_cumsum[i]/p_sum - q_cumsum[i]/q_sum

    if abs(Delta.min()) < abs(Delta.max()):
        ES = Delta.max()
    else:
        ES = Delta.min()
    
    return Delta, ES