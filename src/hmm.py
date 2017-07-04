from utils import *

class HMM:
    """Simple HMM class.  States and observations all integers.
    State 0 the initial state, state N the final state."""

    def __init__ (self, a, o):
        self.a = a
        self.o = o
        self.ns = a.shape[0]
        self.n_obs = o.shape[1]
        assert a.ndim == 2
        assert o.ndim == 2
        assert a.shape[0] == a.shape[1]
        assert a.shape[0] == o.shape[0]

    def initial_state (self):
        return 0
        
    def is_final (self, state):
        return sum(self.a[state,:]) < 1e-10
        
    def num_states (self): return self.ns
        
    def successors (self, sid):
        return set([ i for i in xrange(self.ns) if self.a[sid,i] > 0 ])
    
    def possible_observations (self, sid):
        return set([ i for i in xrange(self.n_obs) if self.o[sid,i] > 0 ])

    def sample (self):
        alls = []
        obs = []
        s = self.initial_state()
        while self.is_final(s):
            s_next = roll_die (self.a[s,:])
            obs_next = self.sample_one_obs (s_next)
            obs.append ((s_next, obs_next))
            s = s_next
        return obs
        
    def sample (self, n):
        return [ self.sample() for i in xrange(n) ]

    def sample_one_obs (self, state):
        return roll_die (self.o[state,:])

    def sample_state (self, state):
        if  self.is_final(state):
            return None
        else:
            return roll_die (self.a[state,:])
        
    def conditional_state_dist (self, s_prev, s_next):
        if s_prev >= 0:
            d1 = self.a[s_prev,:]
        else:
            d1 = [1] * self.ns

        if s_next >= 0:
            d2 = self.a[:,s_next]
        else:
            d2 = [1] * self.ns

        prod = d1*d2
        return prod / sum(prod)
