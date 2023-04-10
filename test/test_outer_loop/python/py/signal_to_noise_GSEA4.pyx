cimport cython

import numpy as np
cimport numpy as np

from libc.math cimport sqrt, abs

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)

cpdef signal_to_noise_GSEA4(
    np.int64_t[:] x, 
    np.float_t[:,:] y):
  
    cdef int i, c0, c1, size, n
    cdef float m0, m1, s0, s1, Ex_0, Ex_1, Ex2_0, Ex2_1
    cdef list zero_locs = []
    cdef list one_locs = []
    cdef list scores = []

    c0 = 0
    c1 = 0
    size = len(x)
  
    for i in range(size):
        if x[i] == 0:
            c0 += 1
            zero_locs.append(i)
        else:
            c1 += 1
            one_locs.append(i)
            
    for i in range(y.shape[0]):
        m0 = 0.0
        for j in zero_locs:
            m0 += y[i, j]
            
        m0 /= c0
            
        m1 = 0.0
        for j in one_locs:
            m1 += y[i, j]
            
        m1 /= c1
        
        s0 = 0.0
        for j in zero_locs:
            s0 += (y[i, j] - m0)**2
        
        s0 /= (c0 - 1)
        
        s1 = 0.0
        for j in one_locs:
            s1 += (y[i, j] - m1)**2
        
        s1 /= (c1 - 1)
    
        s0 = sqrt(s0)
        s1 = sqrt(s1)
        
        if s0 < 0.2 * abs(m0):
            if m0 == 0.0:
                s0 = 0.2
            else:
                s0 = 0.2 * abs(m0)
            
        if s1 < 0.2 * abs(m1):
            if m1 == 0.0:
                s1 = 0.2
            else:
                s1 = 0.2 * abs(m1)   
        
        scores.append((m1 - m0)/(s1 + s0))
    
    return scores
