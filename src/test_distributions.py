import unittest
import distributions
import numpy
from numpy import random
import pwfun

import netutils
import sampling

class TestDistributions (unittest.TestCase):  

    def test_exponential (self):
        exp = distributions.Exponential (scale=0.5)
        xs = exp.sample (100000)
        mu = numpy.mean (xs)
        self.assertTrue (abs(mu - 0.5) < 0.01, "Mean off: was %.4f expected 0.5" % mu)

    def test_exponential_estimate (self):
        xs = [ 1.25, 1.5, 1.75 ]
        exp = distributions.Exponential ()
        exp.estimate (xs)
        self.assertEquals (1.5, exp.mean(), 5)

    def test_gamma_sample (self):
        exp1 = distributions.Exponential (1.5)
        gamm = distributions.Gamma (shape=2, scale=1.5)
        N = 1000
        smp1 = gamm.sample (N)
        smp2a = exp1.sample (N)
        smp2b = exp1.sample (N)
        smp2 = numpy.array(smp2a) + numpy.array(smp2b)
        smp1.sort()
        smp2.sort()
        netutils.check_quantiles (self, smp1, smp2, N)

    def test_gamma_estimate (self):
        gamm = distributions.Gamma (shape=3.5, scale=0.75)
        N = 100000
        smp = gamm.sample (N)
        gamm.estimate (smp)
        self.assertTrue (abs(3.5 - gamm.shape) < 0.1, "Estimation wrong: expected %.4f was %.4f" % (3.5, gamm.shape))
        self.assertTrue (abs(0.75 - gamm.scale) < 0.1, "Estimation wrong: expected %.4f was %.4f" % (0.75, gamm.shape))

    def test_gamma_pdf (self):
        gamm = distributions.Gamma (shape=3, scale=0.5)
        self.assertAlmostEqual (-1.227411, gamm.lpdf(2.0), 5)

    def test_lognormal_estimate (self):
        ln = distributions.LogNormal()
        smp = ln.sample (1000)
        ln.estimate (smp)
        self.assertTrue (abs(ln.meanlog) < 0.1)
        self.assertTrue (abs(1 - ln.sdlog) < 0.1)

    def test_lognormal_estimate2 (self):
        ln = distributions.LogNormal(meanlog=-1, sdlog=2)
        smp = ln.sample (1000)
        ln.estimate (smp)
        self.assertTrue (abs(-1 - ln.meanlog) < 0.1, "Error in mean: Expected %s, got %s" % (-1, ln.meanlog))
        self.assertTrue (abs(2 - ln.sdlog) < 0.1, "Error in sd: expected %s, got %s" % (2, ln.sdlog))

    def test_lognormal_pdf (self):
        ln = distributions.LogNormal(meanlog=1, sdlog=0.5)
        self.assertAlmostEqual (ln.lpdf(1),  -2.225791, 5)
        
        ln2 = distributions.LogNormal(meanlog=0.9, sdlog=0.5)
        ps = [ -18.436309, -6.779582, -3.695824, -2.321165, -1.616153, -1.245645, -1.065585, -1.003215 ]
        xs = [ 0.10, 0.35, 0.60, 0.85, 1.10, 1.35, 1.60, 1.85 ]
        for p,x in zip(ps,xs):
            self.assertAlmostEqual (ln2.lpdf(x), p, 4)

    def test_lognormal_boundary (self):
        ln = distributions.LogNormal(meanlog=1, sdlog=0.5)
        print ln.lpdf(0)
        self.assertTrue (numpy.isinf (ln.lpdf(0)))

    def test_exponential_sample_parameters (self):
        sampling.set_seed(3242)
        for rep in range(10):
            mu = random.exponential(1.0)
            f = distributions.Exponential (mu)
            x = f.sample(10000)            
            params = [ f.sample_parameters(x)[0] for i in xrange(10000) ] 
            self.assertTrue (abs(mu - numpy.mean(params)) < 0.03, "Mismatch: MU %s params %s" % (mu, numpy.mean(params)))

    def test_ln_sample_parameters (self):
        sampling.set_seed(3242)
        for rep in range(10):
            mu = random.uniform(-1, 1)
            sd = random.uniform(0.5, 1.25)
            print rep, mu, sd
            f = distributions.LogNormal (mu, sd)
            x = f.sample(1000)
            print numpy.mean(map(numpy.log, x))
            params = [ f.sample_parameters(x) for i in xrange(10000) ] 
            mu1 = [ p[0] for p in params ]
            sd1 = [ p[1] for p in params ]
            print numpy.mean(mu1)
            print numpy.mean(sd1)
            self.assertTrue (abs(mu - numpy.mean(mu1)) < 0.1, "Mismatch: MU %s params %s" % (mu, numpy.mean(mu1)))
            self.assertTrue (abs(sd - numpy.mean(sd1)) < 0.1, "Mismatch: std %s params %s" % (sd, numpy.mean(sd1)))


    def test_gamma_sample_parameters (self):
        sampling.set_seed(3242)
        for rep in range(1):
            shape = random.uniform(0.5, 3.0)
            scale = random.uniform(0.0, 10.0)
            print "REP", rep, shape, scale

            f = distributions.Gamma (shape, scale)
            x = f.sample(1000)            
            params = [ f.sample_parameters(x) for i in xrange(1000) ] 
            shape1 = [ p[0] for p in params ]
            scale1 = [ p[1] for p in params ]
            for p in params: print "P", " ".join (map(str,p)), p[0]*p[1]

            self.assertTrue (abs(scale - numpy.mean(scale1)) < 0.03, "Mismatch: MU %s params %s" % (scale, numpy.mean(scale1)))
            self.assertTrue (abs(shape - numpy.mean(shape1)) < 0.03, "Mismatch: SHAPE %s params %s" % (shape, numpy.mean(shape1)))

    def test_exponential_kernel (self):
        sampling.set_seed(3242)
        N = 100
        for rep in range(10):
            # pick an Exponential distribution, arbitrarily
            f = distributions.Exponential (1.0)
            self.do_test_kernel (f, N)

    def test_gamma_kernel (self):
        sampling.set_seed(3242)
        N = 100
        for rep in range(10):
            f = distributions.Gamma (2.0, 5.0)
            self.do_test_kernel (f, N)
            f = distributions.Gamma (5.0, 0.75)
            self.do_test_kernel (f, N)

    def test_ln_kernel (self):
        sampling.set_seed(3242)
        N = 100
        for rep in range(10):
            f = distributions.LogNormal (0.5, 0.25)
            self.do_test_kernel (f, N)
            f = distributions.LogNormal (0.0, 1.0)
            self.do_test_kernel (f, N)

    def do_test_kernel (self, f, N):
        data = f.sample (10)
        
        # find special point mu1, check E_mu T(mu1 <== mu)  = p(mu1)
        np = len(f.parameters)
        s_mu = numpy.zeros((np, N))
        for i in range(N):
            s_mu[:,i] = f.sample_parameters(data)
        theta_mean = s_mu.sum(axis=1) / N
        theta_std = 0.2 * theta_mean
        print "Data MEAN %.5f   SD: %.5f" % (numpy.mean(data), numpy.std(data))
        print "Bayes sample MEAN %s  SD: %s" % (theta_mean, theta_std)

        # assumes some kind of CLT-like concentration has happend
        k01 = f.parameter_kernel(theta_mean - 2*theta_std, theta_mean, data)
        k10 = f.parameter_kernel(theta_mean, theta_mean - 2*theta_std, data)
        print "K(mu-1SD: %s ==> mu: %s) %.5f" % (theta_mean - 2*theta_std, theta_mean, k01)
        print "K(mu ==> mu-1SD)", k10
        self.assertTrue (k01 >= k10)

        k21 = f.parameter_kernel(theta_mean + 2*theta_std, theta_mean, data)
        k12 = f.parameter_kernel(theta_mean, theta_mean + 2*theta_std, data)
        print "K(mu+1SD ==> mu)", k21
        print "K(mu ==> mu+1SD)", k12
        self.assertTrue (k21 >= k12)
#        self.assertTrue (abs(k01 - k21) <= 0.01)

    def test_ig_pdf (self):
        ig = distributions.InverseGamma (shape=5.0, scale=2.0)
        self.assertAlmostEquals (-5.896807, ig.lpdf(0.1), 5)

import sys
def main():
    if len(sys.argv) > 1:
        for test_name in sys.argv[1:]:
           suite = unittest.TestLoader().loadTestsFromName("test_distributions.TestDistributions.%s" % (test_name,))
           unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.main()

if __name__ == "__main__":
    main()


