import numpy
from numpy import random
import scipy.special
import netutils

from scipy.special import digamma, polygamma
from scipy.stats import distributions

cdef extern from "math.h":
    double sqrt(double x)
    double exp(double x)
    double log(double x)
    enum: INFINITY

cdef class Distribution:

    # Interface to Python.  Shouldn't need to be overridden
    
    def lpdf (self, x): return self._lpdf (x)
    def dx_lpdf (self, x): return self._dx_lpdf (x)
    def estimate (self, data): self._estimate (data)
    def mean (self): return self._mean ()
    def sample (self, N): return self._sample (N)
     
    # Override these in subclasses
    
    property parameters:
        def __get__ (self): raise NotImplementedError()
        def __set__ (self, v): raise NotImplementedError()
        
    def numParameters (self): raise NotImplementedError()
    
    cdef double _lpdf (self, double x):
        raise NotImplementedError()

    cdef double _quantile (self, double p):
        raise NotImplementedError()
    
    cdef double _dx_lpdf (self, double x):
        raise NotImplementedError()
        
    cdef _estimate (self, data):
        raise NotImplementedError()

    cdef double _mean (self):
        raise NotImplementedError()

    cdef object _sample (self, int N):
        raise NotImplementedError()
        
        
cdef class Exponential (Distribution):
    
    cdef double scale

    def __init__ (self, scale=1.0):
        self.scale = scale
        
    cdef double _lpdf (self, double x):
        return -log(self.scale) - x/self.scale

    cdef double _dx_lpdf (self, double x):
        return -1.0/self.scale
        
    cdef _estimate (self, data):
        mu = numpy.mean (data)
        if len(data) == 0:
            self.scale = 1.0
        # any less than below is unstable with python floats.
        #  once we get to using doubles consistently, we can reduce this
        elif mu < 1e-7:  
            self.scale = 1e-7
        else:
            self.scale = mu
        
    cdef double _mean (self): return self.scale
    
    def std (self): return self.scale
    
    cdef object _sample (self, int N):
        return random.exponential (scale=self.scale, size=N)

    cdef double _quantile (self, double p):
        return -self.scale * log(1.0-p)

    property parameters:
        def __get__ (self): return [self.scale]
        def __set__ (self, v): self.scale = v[0]

    def numParameters (self): return 1

    def sample_parameters (self, data):
        mu = numpy.sum (data)
        N = len(data)
        d = Gamma (N, 1/mu)
        return [ 1.0 / d.sample(1)[0] ]

    # theta_to [scale parameter] ~ IG(N, sum(data))
    def parameter_kernel (self, theta_from, theta_to, data):
        mu = numpy.sum (data)
        N = len(data)
        d = InverseGamma (N, mu)
        lp = d.lpdf(theta_to[0])
        return lp


cdef class Gamma (Distribution):
    """Gamma distribution.  The pdf of a gamma distribution with shape = a and scale = b is
    
          1/(b**a * Gamma(a)) * x^(a-1) * exp(x/b)
        
      Estimation is performed by the method of moments.  By default, a is constrained by a >= 1.
      (Otherwise, the distribution is not convex as a function of x.)  This can be overridden
      by setting self.keep_convex = 0"""
    
    cdef readonly double shape
    cdef readonly double scale
    cdef readonly double log_ba
    cdef bint keep_convex
    
    def __repr__ (self):
        return "<GAMMA SHAPE=%s SCALE=%s/>" % (self.shape, self.scale)
        
    def __init__ (self, shape, scale):
        self.shape = shape
        self.scale = scale
        self.keep_convex = 0
        
    cdef double _lpdf (self, double x):
        if self.shape == 1: # special case: handles x == 0 corretly
            return -x / self.scale - self.shape*log(self.scale)
        else:
            return (self.shape-1)*log(x) - (x / self.scale) - scipy.special.gammaln(self.shape)  -self.shape*log(self.scale) 
    
    cdef double _dx_lpdf (self, double x):
        if self.shape == 1:
            return -1 / self.scale
        else:
            return ((self.shape - 1.) / x) - (1. / self.scale)
        
    # ML, from Tom Minka
    cdef _estimate (self, data):
       cdef double mu, mu_log, a, a_old

       mu = numpy.mean(data)
       mu_log = numpy.mean (map (numpy.log, data))

       # initialization
       a = 0.5 / (log (mu) - mu_log)

       a_old = -numpy.inf
       # iterate
       while abs(a - a_old) >= 0.01:
           a_old = a
           a = 1/(1/a_old + (mu_log - log(mu) + log(a_old) - digamma(a_old)) / (a_old*a_old*(1/a_old - polygamma (1, a_old))))

       b = mu/a

       self.shape = a
       self.scale = b

    # method of moments
    cdef _estimate_mom (self, data):
        if  len(data) == 0:
            self.shape = self.scale = 1.0
            return
            
        cdef double mu, sigma, m2
        mu = numpy.mean (data)
        sigma = numpy.std (data)
        
        if sigma < 1e-10: sigma = 1e-5
        
#        m2 = sigma*sigma + mu*mu       # VX + (EX)^2 = E(X^2)
#        print "ESTIMATING GAMMA: N %d MU %s V %s" % (len(data), mu, sigma*sigma)
        self.shape = (mu*mu) / (sigma*sigma)
        self.scale = (sigma*sigma) / mu     
           
        if self.keep_convex and self.shape < 1:
            self.shape = 1.0
            self.scale = mu
            
        if self.shape == numpy.inf:
            print "ERROR: Infinite shape (mu: %s sig: %s)" % (mu, sigma)
            
        self.shape = max(1e-10, self.shape)
        self.scale = max(1e-10, self.scale)
    
    cdef object _sample (self, int N):
        return random.gamma (self.shape, size=N, scale=self.scale)

    property parameters:
        def __get__ (self): return [self.shape, self.scale]
        def __set__ (self, v): 
            self.shape, self.scale = v

    def numParameters (self): return 2

    # sample_parameters: improper prior, slice sampler
    
    def sample_parameters (self, data):
        return netutils.sample_gamma_posterior (data, self.shape, self.scale)

    def parameter_kernel (self, theta_from, theta_to, data):
        shape0, scale0 = theta_from
        shape1, scale1 = theta_to
        return netutils.gamma_transition_kernel (data, shape0, scale0, shape1, scale1)

cdef class InverseGamma (Distribution):
    """Inverse gamma distribution.  The pdf of the reciprocal of a gamma distribution with shape = a and scale = b."""
    
    cdef readonly double shape
    cdef readonly double scale
    cdef object gamma

    def __repr__ (self):
        return "<IG SHAPE=%s SCALE=%s/>" % (self.shape, self.scale)
        
    def __init__ (self, shape, scale):
        self.shape = shape
        self.scale = scale
        self.gamma = Gamma (shape, scale)
        
    cdef double _lpdf (self, double x):
        return -(self.shape+1)*log(x) - (self.scale / x) - scipy.special.gammaln(self.shape)  + self.shape*log(self.scale) 

    cdef _estimate (self, data):
       self.gamma.estimate ([ 1.0/x for x in data])
       self.shape = self.gamma.shape
       self.scale = self.gamma.scale

    cdef object _sample (self, int N):
        s1 = random.gamma (self.shape, size=N, scale=self.scale)
        return [ 1.0 / s1 for x in s1 ]

    property parameters:
        def __get__ (self): return [self.shape, self.scale]
        def __set__ (self, v): 
            self.shape, self.scale = v

    def numParameters (self): return 2

cdef double HALF_LOG_2PI = -0.918938533204673  # - log(sqrt(2*pi))    

cdef class LogNormal (Distribution):
    
    cdef readonly double meanlog
    cdef readonly double sdlog
    
    def __init__ (self, meanlog=0, sdlog=1):
        self.meanlog = meanlog
        self.sdlog = sdlog
        
    cdef double _mean (self):
        return exp(self.meanlog + self.sdlog*self.sdlog/2)
        
    def std (self):
        return sqrt( (exp(self.sdlog*self.sdlog) - 1) * exp(2*self.meanlog + self.sdlog*self.sdlog) )
        
    cdef object _sample (self, int N):
        return numpy.exp (random.normal (self.meanlog, self.sdlog, size=N))
        
    cdef double _quantile (self, double p):
        raise Exception ("NYI.")
    
    cdef double _lpdf (self, double x):
        if x == 0:
            return -numpy.inf
        else:
            return HALF_LOG_2PI - log(self.sdlog) - log(x) - ((log(x) - self.meanlog) ** 2) / (2*self.sdlog*self.sdlog)

    cdef double _dx_lpdf (self, double x):
        return -1.0/x - (log(x) - self.meanlog) / (x*self.sdlog*self.sdlog)
        
    # hack: discard 0's
    cdef compute_meansd_log (self, data):
        if len(data) == 0: 
          return (0,1)
	    
        cdef double x, mean, sd
        cdef int N = 0
	    
        mean = 0
        for x in data:
            if x < 1e-50: continue
            mean += log(x)
            N += 1
        mean = mean / N
        
        sd = 0
        for x in data:
            if x < 1e-50: continue
            sd += (log(x) - mean)*(log(x) - mean)
        sd = sqrt (sd / N)
        
        # don't allow zeros
        if sd < 1e-10: sd = 1e-10

        return mean,sd

    cdef _estimate (self, data):
        mean, sd = self.compute_meansd_log (data)
        self.meanlog = mean
        self.sdlog = sd

    def sample_parameters (self, data):
        mean, sd = self.compute_meansd_log (data)
        if sd <= 1e-5: return [mean, sd]

        # limit of a gamma-normal prior: p(mean, log sd^2) = 1
        ns = len(data)
        if ns == 0: return [0.0, 1.0]

        alpha = (ns+1) / 2
        beta = 2 / (ns*sd*sd)

        tau_new = random.gamma(alpha, beta)
        sd_new = sqrt(1/tau_new)
        mu_new = random.normal(mean, sd_new)
        return [mu_new, sd_new]

    # TODO: Check this, I suspect that it's wrong
    def parameter_kernel (self, theta0, theta1, data):
        mean, sd = self.compute_meansd_log (data)
        if sd <= 1e-5: return [mean, sd]

        # limit of a gamma-normal prior: p(mean, log sd^2) = 1
        ns = len(data)
        if ns == 0: return [0.0, 1.0]

        alpha = (ns+1) / 2
        beta = 2 / (ns*sd*sd)

        mu_new, sd_new = theta1
        lpdf_tau = Gamma(alpha, beta).lpdf(1/(sd_new**2))
        lpdf_mu = log(distributions.norm.pdf (mu_new, mean, sd_new))

        return lpdf_mu + lpdf_tau

    property parameters:
        def __get__ (self): return [self.meanlog, self.sdlog]
        def __set__ (self, v): self.meanlog, self.sdlog = v

    def numParameters (self): return 2


cdef class Dirac (Distribution):
    
    cdef double loc

    def __init__ (self, loc):
        self.loc = loc
        
    cdef _estimate (self, data):
        mu = numpy.mean (data)
        self.loc = mu
        
    cdef double _mean (self): return self.loc
    
    cdef object _sample (self, int N):
        xs = []
        for i from 0 <= i < N: xs.append(self.loc)
        return xs 

    cdef double _quantile (self, double p):
        return self.loc

    property parameters:
        def __get__ (self): return [self.loc]
        def __set__ (self, v): self.loc = v[0]

    def numParameters (self): return 1
