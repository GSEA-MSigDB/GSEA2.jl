cimport cython
import numpy as np
cimport numpy as np

from libc.math cimport abs, pow

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)

cpdef inline compute_KS_AUC_enrichment(
    np.ndarray[np.float_t, ndim=1] S_h,
    np.ndarray[np.int_t, ndim=1] I_h,
    double alpha):

    cdef int i, n = len(S_h)
    cdef np.float_t ES, p_sum = 1.0e-16, q_sum = 1.0e-16, Delta_sum = 0
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
        Delta_sum += Delta[i]

    ES = Delta_sum/n
    
    return Delta, ES