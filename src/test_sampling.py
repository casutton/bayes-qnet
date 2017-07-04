import sampling
import unittest
import numpy.random
import scipy.stats
import scipy.special
from math import log, exp
import pwfun
import netutils
import distributions

import sys
import pylab

def gen_log_gamma (shape, scale):
    return lambda x: log(scipy.stats.gamma.pdf(x, shape, scale=scale))
    
def gen_dlog_gamma (shape, scale):
    # gamma : G(a)/b^a * x^(a-1) * exp(x/b)
    """Returns a function that gives the derivative of the gamma log-density with respect to X."""
    def fn(x): 
        ret = ((shape-1.0)/x) - (1.0/scale)
        return ret
    return fn

class TestSampling (unittest.TestCase):  

    def test_gamma_deriv (self):
        params = [ (1.0, 0.5), (1.0, 2.0), (5.0, 0.5), (5.0, 2.5) ]
        xs = [ 0.1, 0.5, 1.0, 5.0 ]
        eps = 1e-6
        
        for shape,scale in params:
            logpdf = gen_log_gamma (shape, scale)
            deriv = gen_dlog_gamma (shape, scale)
            for x in xs:
                dx_analytic = deriv(x)
                dx_empirical = (logpdf(x+eps) - logpdf(x)) / eps
                diff = abs((dx_analytic - dx_empirical))
                self.assertTrue (diff < 1e-3, "ERROR: x %.4f, shape %.4f scale %.4f, pdf1 %.4f pdf2 %.4f analytic: %.4f empirical: %.4f" % (x, shape, scale, logpdf(x+eps), logpdf(x), dx_analytic, dx_empirical))
                
    def test_gamma_ub (self):
        """Test that UB actually puts an upper bound on the function."""
        lpdf = gen_log_gamma (2.0, 10)
        dlpdf = gen_dlog_gamma (2.0, 10)
        fu = pwfun.Pwlin (0, 30, lpdf, dlpdf)
        print fu.dump()
        
        xs = range(1,10)
        
        for x in xs:
            self.assertTrue (fu(x) >= lpdf(x), "X: %.4f  pdf(X): %.4f  UB(X): %.4f\n%s" % (x, lpdf(x), fu(x), fu.dump()))
            
            
            
    def test_gamma_via_gamma (self):
        """ Test of calibration of the check_quantiles test."""
        
        N = 1000
        
        params = [ (2.0, 5.0, 10.0), (5.0, 0.25, 0.75) ]
        for shape, scale, max_x in params:
        
            nrej = 0
            s1 = []
            while len(s1) < N:
                x = numpy.random.gamma(shape, scale) 
                if x < max_x: s1.append(x) 
                else: nrej += 1

            s2 = []
            while len(s2) < N:
                x = numpy.random.gamma(shape, scale) 
                if x < max_x: s2.append(x) 
                else: nrej += 1
            print "REJECTION: Number of rejections %d" % nrej

            s1.sort()
            s2.sort()
        
            netutils.check_quantiles (self, s1, s2, N)

            
    def test_ars_via_gamma (self):
        """ Test by sampling from a truncated gamma distribution. """
        
#        N = 10000
        N = 1000
        
        params = [ (2.0, 5.0, 10.0), (5.0, 0.25, 0.75) ]
#        params = [ (0.5, 0.75, 5.0 ) ]
        for shape, scale, max_x in params:
        
            print "==== SAMPLING : Gamma (shape: %.4f  scale: %.4f  x_max:  %.4f) ===" % (shape, scale, max_x)
            
            nrej = 0
            s1 = []
            while len(s1) < N:
                x = numpy.random.gamma(shape, scale) 
                if x < max_x: s1.append(x) 
                else: nrej += 1
            print "REJECTION: Number of rejections %d" % nrej

            lpdf = gen_log_gamma (shape, scale)
            dlpdf = gen_dlog_gamma (shape, scale)
        
            sampling.reset_statistics()
            
            s2 = sampling.arms (lpdf, dlpdf, 0.0, max_x, n=N)
            
            s1.sort()
            s2.sort()
        
            netutils.check_quantiles (self, s1, s2, N)
            
            stats = sampling.statistics()
            self.assertEquals (N, stats['values_generated'])
            self.assertEquals (0, stats['mh_rejects'])

    def test_arms_via_lognormal (self):
        """ Test by sampling from a lognormal distribution."""

        N = 5000

        params = [ (0.0, 1.0, numpy.inf), (-2.0, 0.5, numpy.inf), (2.0, 2.0, 1000.0) ]

        for mu, sig, x_max in params:

            print "==========================\nSampling lognormal mean %.4f  sd %.4f" % (mu, sig)
            
            ln = distributions.LogNormal (mu, sig)
            s_ex = []
            while len(s_ex) < N:
                x = ln.sample (1)
                if x < x_max: s_ex.append (x)

            sampling.reset_statistics()

            s_ars = sampling.arms (ln.lpdf, ln.dx_lpdf, 0.0, x_max, n=N, num_points=5)

            pylab.hold(False)
            pylab.plot (s_ars)
            pylab.savefig('iter.png')
            
            pylab.plot (s_ars[1:N-1], s_ars[2:N], 'ro')
            pylab.savefig('lag1.png')
            
            s_ex.sort()
            s_ars.sort()
            
            stats = sampling.statistics()
            print "ARMS: Num MH rejections = %d" % (stats['mh_rejects'])
            
            self.assertEquals (N, stats['values_generated'])
            self.assertTrue (0 < stats['mh_rejects'])

            netutils.check_quantiles (self, s_ex, s_ars, N)


    def test_sample_one_section(self):
        fn = lambda x: -x # LOG density
        df = lambda x: -1 # d LOG density / dx
        pwf = pwfun.Pwlin(0.0, 5.0, fn, df)
        
        N = 10000
        max_x = pwf.knots()[1]
        
        s1 = []
        while len(s1) < N:
            x = numpy.random.exponential(1)
            if x < max_x: s1.append (x)
            
        print "Exponential sampling done"
        
        s2 = []
        while len(s2) < N:
            s2.append (pwf.sample_from_section_exp (0))

        s1.sort()
        s2.sort()
        
        netutils.check_quantiles (self, s1, s2, N)


    def test_sample_one_section2 (self):
        fn = lambda x: -x # LOG density
        df = lambda x: -1 # d LOG density / dx
        pwf = pwfun.Pwlin(0.0, 5.0, fn, df)

        N = 10000
        min_x = pwf.knots()[1]
        max_x = pwf.knots()[2]

        s1 = []
        while len(s1) < N:
            x = numpy.random.exponential(1)
            if min_x < x and x < max_x: s1.append (x)

        s2 = []
        while len(s2) < N:
            s2.append (pwf.sample_from_section_exp (1))

        s1.sort()
        s2.sort()

        netutils.check_quantiles (self, s1, s2, N)
    

    def test_slice_via_normal (self):
        """ Test by sampling from a normal distribution."""

        N = 1000

        params = [ (0.0, 1.0), (-2.0, 0.5) ]

        for mu, sig in params:

            print "==========================\nSampling normal mean %.4f  sd %.4f" % (mu, sig)

            s_ex = numpy.random.normal (loc=mu, scale=sig, size=N)

            sampling.reset_statistics()

            ss = sampling.slice (lambda x: -((x - mu)**2)/(2*sig*sig), mu, N)

            pylab.hold(False)
            pylab.plot (ss)
            pylab.savefig('iter.png')

            pylab.plot (ss[1:N-1], ss[2:N], 'ro')
            pylab.savefig('lag1.png')

            s_ex.sort()
            ss.sort()

            stats = sampling.statistics()

            netutils.check_quantiles (self, s_ex, ss, N)


        def test_sample_one_section(self):
            fn = lambda x: -x # LOG density
            df = lambda x: -1 # d LOG density / dx
            pwf = pwfun.Pwlin(0.0, 5.0, fn, df)

            N = 10000
            max_x = pwf.knots()[1]

            s1 = []
            while len(s1) < N:
                x = numpy.random.exponential(1)
                if x < max_x: s1.append (x)

            print "Exponential sampling done"

            s2 = []
            while len(s2) < N:
                s2.append (pwf.sample_from_section_exp (0))

            s1.sort()
            s2.sort()

            netutils.check_quantiles (self, s1, s2, N)


        def test_sample_one_section2 (self):
            fn = lambda x: -x # LOG density
            df = lambda x: -1 # d LOG density / dx
            pwf = pwfun.Pwlin(0.0, 5.0, fn, df)

            N = 10000
            min_x = pwf.knots()[1]
            max_x = pwf.knots()[2]

            s1 = []
            while len(s1) < N:
                x = numpy.random.exponential(1)
                if min_x < x and x < max_x: s1.append (x)

            s2 = []
            while len(s2) < N:
                s2.append (pwf.sample_from_section_exp (1))

            s1.sort()
            s2.sort()

            netutils.check_quantiles (self, s1, s2, N)


    def test_sample_one_section3 (self):
        fn = lambda x: log(2) + 2*x # LOG density
        df = lambda x: 2 # d LOG density / dx
        
        pwf = pwfun.Pwlin(0.0, 5.0, fn, df)

        N = 10000
        max_x = pwf.knots()[1]

        s1 = []
        M = exp(fn(max_x))
        while len(s1) < N:
            x = numpy.random.uniform(0, max_x)
            unif = numpy.random.rand()
            if unif < exp(fn(x)) / M: s1.append (x)

        print "Exponential sampling done"

        s2 = []
        while len(s2) < N:
            s2.append (pwf.sample_from_section_exp (0))

        s1.sort()
        s2.sort()

        netutils.check_quantiles (self, s1, s2, N)
        
    def test_roll_die (self):
        p = [ 0.5, 1, 0.5 ]
        N = 10000
        smpl = [ sampling.roll_die_unnorm (p) for i in xrange(N) ]
        x0 = sum ([ i==0 for i in smpl])
        x1 = sum ([ i==1 for i in smpl])
        x2 = sum ([ i==2 for i in smpl])
        self.assertEquals (N, x0+x1+x2)
        self.assertTrue (abs(x0 - N/4) < (2*N/4*0.25*0.75))
        self.assertTrue (abs(x1 - N/2) < (N*0.5*0.5))
        self.assertTrue (abs(x2 - N/4) < (2*N/4*0.25*0.75))
        
    def test_pwfun_ars (self):
        N = 10000
        fn = pwfun.Pwfun ([0, 1, 100], [ f1, f2 ], [ df1, df2 ])        
        x_ars = [ sampling.arms_pwfun(fn) for i in xrange(N) ]
        x_true = [ sample_f() for i in xrange(N) ]
        x_ars.sort()
        x_true.sort()
        netutils.check_quantiles (self, x_ars, x_true, N)
    
    def test_inf_gaussian (self): 
        fn = lambda x: -0.5*x*x
        dfn = lambda x: -x
        flnorm = pwfun.Pwfun ([0, numpy.inf], [fn], [dfn])

        N = 10000
        
        smp = [ sampling.arms_pwfun(flnorm) for i in xrange(N) ]
        smp.sort()
        
        true_sample = []
        while len(true_sample) < N:
            x = numpy.random.standard_normal()
            if x > 0: true_sample.append (x)
        true_sample.sort()
        
        netutils.check_quantiles (self, smp, true_sample, len(smp))
        
    
    def test_sample_inf_section(self):
        fn = lambda x: -x # LOG density
        df = lambda x: -1 # d LOG density / dx
        pwf = pwfun.Pwlin(0.0, numpy.inf, fn, df)

        print pwf.dump()
        
        self.assertEquals (3, len(pwf.knots()))
        self.assertEquals (numpy.inf, pwf.knots()[2])
        
        N = 100
        xmin = pwf.knots()[1]
        s1 = [ xmin+numpy.random.exponential(1) for i in xrange(N) ]

        s2 = []
        while len(s2) < N:
            s2.append (pwf.sample_from_section_exp (1))

        s1.sort()
        s2.sort()

        netutils.check_quantiles (self, s1, s2, N)

    def test_slice_ln (self):         
        ln = distributions.LogNormal (-2, 2)
        N = 10000
#        smp = sampling.slice (ln.lpdf, 100, lower=50, upper=200, N)
        smp = sampling.slice (ln.lpdf, 1.0, lower=2.672846e-05, upper=6.852485e+02, N=N)
        theoretical = [ 0.00000000, 0.01042964, 0.02514132, 0.04741574, 0.08153734, 0.13533528, 0.22462885, 0.38627761, 0.72850737, 1.75611350, numpy.inf ]
        netutils.check_quantiles (self, smp, theoretical, N)
        print "Please ensure there were no warnings from _slice about excessive steps."
        
    def test_stupid_reject (self):
        fn = lambda x: -2*x
        df = lambda x: -2
        pwf = pwfun.Pwfun([0.0, 5.0], [fn], [df])
        N = 100

        s1 = []
        while len(s1) < N:
            x = numpy.random.exponential (0.5)
            if x < 5.0:
                s1.append(x)

        s2 = []
        while len(s2) < N:
            s2.append (sampling.rejection (pwf))

        netutils.check_quantiles (self, s1, s2, N)


        
def sample_f ():
    u = numpy.random.rand()
    if u < 0.25:
        return numpy.random.rand()
    else:
        return 1 + numpy.random.exponential (scale=2)
        
def f1 (x):
    if 0 <= x <= 1:
        return log(0.25)
    else:
        return 0

def df1(x): return 0
    
def f2 (x):
    if x > 1:
        return log(0.75*0.5) - 0.5*(x-1)
    else:
        return 0
        
def df2 (x):
    return -0.5
    
if __name__ == "__main__":
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
            suite = unittest.TestLoader().loadTestsFromName("test_sampling.TestSampling.%s" % (test_name,))
            unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()
    
