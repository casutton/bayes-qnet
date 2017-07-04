cdef class Queue

from distributions cimport Distribution
from pwfun cimport Pwfun
from arrivals cimport Arrivals
from arrivals cimport Event


cdef Pwfun _pair_proposal (Arrivals arrv, Event e0, Event e1)
cdef Pwfun _departure_proposal (Arrivals arrv, Event e0)

cdef class FactorTemplate:
    cdef object estimate (self, evtl)
    cdef object serviceLikelihood (self, Event evt)
    cdef object serviceDerivative (self, Event evt)
    cdef void resample_auxiliaries (self, Event evt)
    cdef void initialize_auxiliaries (self, Event evt)

cdef class Queue:
    cdef readonly object name
    cdef readonly FactorTemplate service

    cdef Pwfun arrivalLik (self, Arrivals arrv, Event evt)

    cdef Pwfun departureLik (self, Arrivals arrv, Event evt)

    cdef double allEventLik (self, Arrivals arrv) except *

    cdef diffListForArrival (self, Event evt, double d)
    cdef diffListForDeparture (self, Event evt, double d)
    cpdef double likelihoodDelta (self, Arrivals arrv, diffList) except *

    cpdef recompute_service (self, Event e)
    cpdef recompute_departure (self, Arrivals arrv, Event e)
    cdef double departure_of_service  (self, Event e, double s)
    
    cpdef double service_lpdf (self, Event e)
