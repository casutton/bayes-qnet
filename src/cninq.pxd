
cdef class SortedDouble
cdef class OverlayIterator
cdef class Overlay

cdef class Ninq:

    cdef SortedDouble a
    cdef SortedDouble d
    cdef int Nmax

    cpdef add_birth_death (self, double t_a, double t_d)
    cpdef move_arrival (self, double a0, double a1)
    cpdef move_departure (self, double a0, double a1)
    cpdef int N (self, double t)
    cpdef OverlayIterator interval_iterator (self, double l, double r)
    cpdef int contains_arrival (self, double v)
    cpdef int contains_departure (self, double v)


cdef class Overlay:

    cdef Ninq inner
    cdef SortedDouble plus
    cdef SortedDouble minus

    cpdef move_arrival (self, double a0, double a1)
    cpdef move_departure (self, double a0, double a1)
    cpdef OverlayIterator interval_iterator (self, double l, double r)
    cdef OverlayIterator _interval_iterator (self, double l, double r)
    cpdef int N (self, double t)


cdef class SortedDouble:
    cdef double *val
    cdef int N
    cdef int capacity

    cpdef double item (self, int i)
    cpdef add_time (self, double x)
    cpdef move_time (self, double x0, double x1)
    cpdef int num_lte (self, v)
    cdef int bisect (self, double v)

 

cdef class OverlayIterator:

    cdef SortedDouble a
    cdef SortedDouble d
    cdef SortedDouble plus
    cdef SortedDouble minus
    
    cdef double r

    cdef int i_a
    cdef int i_d
    cdef int i_plus
    cdef int i_minus
    cdef int is_a
    cdef int is_d
    cdef int is_plus
    cdef int is_minus

    cdef double t0
    cdef double t1
        
    cpdef int has_next (self)
    cpdef double T0 (self)
    cpdef double T1 (self)
    cpdef int N (self)
    cpdef advance (self)
