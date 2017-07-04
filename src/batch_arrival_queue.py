import queues


class BatchArrivalQueue (Queue):

    def __repr__ (self):
        return "<QUEUE TYPE='Batch' NAME='%s'>\n<SERVICE>%s</SERVICE><BATCH_SIZE>%s</BATCH_SIZE>\n</QUEUE>" % (self.name, self.service, self.batch_size)

    def as_yaml (self):
        return "  - { name: %s, type: BATCH, service: %s, size: %s }" % (self.name, self.service.as_yaml(), self.size.as_yaml())

    def when_clear (self, Arrivals arrv):
        cdef Event last_evt = arrv.final_arrival_by_queue (self)
        return last_evt.d if last_evt else 0.0
        
    def previous_departing (self, Event e):
        return e.prev_byq # previous in arrival order
        
    cpdef recompute_service (self, Event e1):
        cdef Event e0, e2
        e0 = e1.prev_byq
        e2 = e1.next_byq
        
        if e0:
            e1.s = e1.d - max(e1.a, e0.d)
        else:
            e1.s = e1.d - e1.a
            
        if e2:
            e2.s = e2.d - max(e2.a, e1.d)
    
    cdef double departure_of_service  (self, Event e, double s):
        cdef Event byq = e.prev_byq
        cdef double d_prev = byq.d if byq else 0.0
        return s + max (e.a, d_prev)
    
    cpdef recompute_departure (self, Arrivals arrv, Event e1):
        cdef Event e_prev = e1.prev_byq
        cdef double d_prev = e_prev.d if e_prev else 0
        e1.set_departure (e1.s + max(e1.a, d_prev))
    
    def validate_service (self, evt):
        cdef Event e = evt
        d_prev = (e.prev_byq.d if e.prev_byq else 0)
        expected = e.d - max(e.a, d_prev)
        assert abs(e.s - expected) < 1e-10, "Invalid service time (expected %s was %s)\n  %s\n  %s" % (expected, e.s, e.prev_byq, e)
        if e.prev_byq:
            assert e.a >= e.prev_byq.a, "Invalid arrival times for events PREV: %s CURR: %s\n  %s\n  %s" % (e.prev_byq.a, e.a, e.prev_byq, e)
            
    cdef Pwfun arrivalLik (self, Arrivals arrv, Event evt): 
        cdef Event e_prev, e_next
        e_prev = evt.prev_byq
        e_next = evt.next_byq
        
        if e_prev:
            d_prev = e_prev.d
            L = e_prev.a
        else:
            L = 0
            d_prev = 0
            
        if e_next:
            U = min(e_next.a, evt.d)
        else:
            U = evt.d
            
        cdef ArrivalFn afn
        afn = ArrivalFn (self.service, evt, evt.d, d_prev) 
        if DEBUG: 
            print afn
            print "@ L(=%.17g)  (diff=%.17g) is %s or %s" % (L, L-afn.L(), afn.f1(L), afn.f2(L))
        
        if d_prev < L:
            return Pwfun ([L, U], [ afn.f2 ], [ afn.df2 ])
        elif U < d_prev:
            return Pwfun ([L, U], [ afn.f1 ], [ afn.df1 ])
        else:
            return Pwfun ([L, d_prev, U], [ afn.f1, afn.f2 ], [ afn.df1, afn.df2 ])
        

    cdef Pwfun departureLik (self, Arrivals arrv, Event evt):
        cdef Event e_prev, e_next
        cdef double L, A, U
        
        e_prev = evt.prev_byq
        e_next = evt.next_byq
        
        if e_prev: 
            L = max(evt.a, e_prev.d) 
        else: 
            L = evt.a            

        cdef DepartureFn dfun
        cdef FinalDepartureFn fdfun
        
        if e_next:             
            A = e_next.a
            U = e_next.d 
            dfun = DepartureFn (self.service, evt, e_next, L, A, U)
            if DEBUG: 
                print dfun
                print "@ L(=%.17g)  (diff=%.17g) is %s or %s" % (L, L-dfun.L(), dfun.f1(L), dfun.f2(L))
            
#            print "Creating departure fn:\n  %s\nL: %8.4f A: %8.4f U: %8.4f\ndfun: %s\n" % (format_six(evt), L, A, U, dfun)
            if A < L:
                return Pwfun([L, U], [dfun.f2], [dfun.df2])
            elif U < A:                
                return Pwfun([L, U], [dfun.f1], [dfun.df1])
            else:
                return Pwfun([L, A, U], [ dfun.f1, dfun.f2 ], [ dfun.df1, dfun.df2 ])
                
        else:
            fdfun = FinalDepartureFn (self.service, evt, L)
            max_x = numpy.inf
            return Pwfun ([L, max_x], [ fdfun.f ], [ fdfun.df ])
           
    cdef diffListForArrival (self, Event evt, double a):
        cdef Event e_new = evt.dup_with_structure()
        e_new.a = a
        self.recompute_service (e_new)
        return [e_new]
        
    cdef diffListForDeparture (self, Event evt, double d):
        cdef double d_prev = evt.prev_byq.d if evt.prev_byq else 0.0
        cdef Event e_new = evt.dup_with_structure()
        e_new.d = d
        e_new.s = d - max(e_new.a, d_prev)
        lst = [e_new]
        
        cdef Event byq = e_new.next_byq
        if byq:
            byq = byq.dup_with_structure()
            byq.prev_byq = e_new
            e_new.next_byq = byq
            byq.s = byq.d - max(byq.a, d)
            lst.append (byq)
            
        return lst      

    def select_job_for_service (self, evt_list):
        return evt_list[0]
