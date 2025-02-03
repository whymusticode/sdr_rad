
from libc.math cimport floor
import numpy as np
cimport numpy as np
from cython cimport boundscheck, wraparound

# Disable bounds checking and wraparound to optimize for speed
@boundscheck(False)
@wraparound(False)
def resampstr(np.ndarray[np.float64_t, ndim=1] p, int m=-1):
    cdef int n = p.shape[0]
    if m == -1:
        m = n

    cdef np.ndarray[np.float64_t, ndim=1] pn = np.empty(n, dtype=np.float64)
    cdef np.ndarray[np.int32_t, ndim=1] s = np.zeros(m, dtype=np.int32)
    cdef np.ndarray[np.float64_t, ndim=1] r = np.random.rand(n)
    
    cdef double sum_p = np.sum(p)
    cdef double c = 0.0
    cdef int i, k = 0

    # Normalize p to sum to m
    for i in range(n):
        pn[i] = p[i] / sum_p * m

    # Resampling logic
    for i in range(n):
        c += pn[i]
        if c >= 1.0:
            a = floor(c)
            c -= a
            for j in range(int(a)):
                s[k] = i + 1  # MATLAB uses 1-based indexing
                k += 1
        if k < m and c >= r[k]:
            c -= 1.0
            s[k] = i + 1
            k += 1

    # Handle last element
    if s[m - 1] == 0:
        mode_val = np.argmax(np.bincount(s))
        s[m - 1] = mode_val

    return s