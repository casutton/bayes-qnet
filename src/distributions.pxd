cdef class Distribution:
    cdef double _lpdf (self, double x)
    cdef double _mean (self)
    cdef object _sample (self, int N)
    cdef _estimate (self, data)
    cdef double _dx_lpdf (self, double x)
    cdef double _quantile (self, double p)
    
