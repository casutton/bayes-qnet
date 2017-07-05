cdef class Event
cdef class Arrivals

from queues cimport Queue
from qnet cimport Qnet
cimport cython

cdef class Vars:
    cdef object map
        
cdef class Event:
    cdef readonly int eid
    cdef readonly int tid
    cdef Queue q
    cdef readonly int qid
    cdef readonly int state
    cdef public double a
    cdef public double s
    cdef public double c
    cdef public double d
    cdef public int obs_a
    cdef public int obs_d
    cdef object prev_byt
    cdef object next_byt
    cdef object prev_byq
    cdef object next_byq
    cdef object arrv

    # For G/G/k queues
    cdef public int proc 
    cdef int num_proc
    cdef double d_prev
    cdef double* d_proc_prev
    
    # used by mixture distributions, etc.
    cdef readonly Vars auxiliaries
    
    cdef set_service (self, double v)
    cdef set_departure (self, double d)

    cpdef Event duplicate (self)
    cdef Event dup_with_structure (self)
    cdef int update_from_solution (self, lp_soln) except -1
    cdef copy_dtimes (self, Event other)

    #private
    cdef init_dtimes (self)
    cdef Event _dup_inner (self, int dup_aux)

@cython.final
cdef class Arrivals:

    cdef Qnet net
    cdef object byq
    cdef object final_byq
    cdef object byt
    cdef object events
    cdef object _initial_events
    cdef object queue_orders
    cdef object ninq
    cdef object inq
    cdef object a2e_cache
    cdef object c2e_cache
    cdef int cap

    cpdef inline Event event (self, int eid)
    cpdef inline int num_events (self)

    cdef _is_initial (self, Event evt)
    cdef _is_final (self, Event evt)
