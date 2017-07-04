cdef class Qnet

from arrivals cimport Arrivals
from arrivals cimport Event

cdef class Qnet:
    
    cdef readonly object queues     # List of Queues
    cdef readonly object templates    # List of all object that have parameters
    cdef readonly object fsm     # HMM object
    
    cdef readonly object universe # for auxiliary variables :: maps names to domains
    
    cdef readonly object qname2id
    cdef readonly object sname2id
    cdef readonly int eidx

    # private
    cdef int gibbs_resample_final (self, Arrivals arrv, Event evt) except -1
    cdef int gibbs_resample_pair (self, Arrivals arrv, Event e0, Event e1) except -1
    
    # private, used by slice sampler
    cdef double determine_upper (self, Event evt)
    cdef object dfn_range (self, Arrivals arrv, Event evt)
    cpdef double log_prob (self, Arrivals arrv) except *

