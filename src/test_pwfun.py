#!/usr/bin/python

import unittest
import mytime
from math import exp, log, sqrt
import numpy

import pwfun
import distributions
import netutils

class TestPwfun (unittest.TestCase):  

    def test_exponential (self):
        Zs = self.fexp1.integrate_sections()
        
        self.assertEqual (3, len(Zs))
        self.assertAlmostEqual (Zs[0], 0.6321206, 5)
        self.assertAlmostEqual (Zs[1], 0.05850982, 5)
        self.assertAlmostEqual (Zs[2], 0.0007851141, 5)

    def test_exp_time (self):
        N = 1000
        tmr = mytime.timeit()
        for i in xrange(N): Z = self.fexp1.integrate()
        tmr.total ("Integrated %d three-piece fns" % N)
        
    def test_pwlin_bound (self):
        xs = [-1, -0.5, 0, 1, 2, 4]
        for x in xs:
            self.assertTrue (self.flin(x) >= -x*x*x*x)
        
    def test_pwlin_knots (self):
        knots_expected = [ -5, -1, 0.25, 0.75, 5 ]
        self.assertEquals (5, len(self.fsq.knots()))        
        for k1,k2 in zip(self.fsq.knots(), knots_expected):
            self.assertAlmostEqual (k1, k2, 10)
    
    def test_pwlin_values (self):
        self.assertAlmostEqual (self.fsq(-4), -12, 10)        
        self.assertAlmostEqual (self.fsq(-2.5), -6, 10)
        self.assertAlmostEqual (self.fsq(-1), 0, 10)
        self.assertAlmostEqual (self.fsq(0), 0, 10)    
        self.assertAlmostEqual (self.fsq(0.5), -0.25, 10)    
        self.assertAlmostEqual (self.fsq(0.75), -0.5, 10)
        self.assertAlmostEqual (self.fsq(1.5), -2, 10)
        self.assertAlmostEqual (self.fsq(3), -5, 10)
            
    def test_pwlin_integral (self):
        print (self.fsq.dump())
        Z1 = self.fsq.integrate_sections ()
        z_exp = [ -32, 0, -0.125, -20.1875 ]
        self.assertEquals (len(Z1), len(z_exp))
        for x1,x2 in zip(Z1,z_exp):
            self.assertAlmostEqual (x2, x1, 5)    

    def test_pwlin_expint (self):
        Z1 = self.floglin.integrate_sections ()
        Z2 = self.fsq.integrate_exp_sections ()
        for x1,x2 in zip(Z1,Z2):
            self.assertAlmostEqual (x1, x2, 5)
        
    def test_argmin_knottiness (self):
        fsq = pwfun.Pwlin (-5, 5, lambda x: -x*x, lambda x: -2*x)
        self.check_knots (fsq)

        fsq.argmin (0)
        self.check_knots (fsq)

        fsq.argmin (0.5)
        self.check_knots (fsq)
        
    def test_argmin_knot_past_end (self):
        """Checks that Pwlin.argmin dosen't break in the edge case that the new point is past all existing tangent points"""
        fsq = pwfun.Pwlin (-5, 5, lambda x: -x*x, lambda x: -2*x)
        fsq.argmin (4)
        self.check_knots(fsq)

    def test_argmin_knot_before_beginning (self):
        """Checks that Pwlin.argmin dosen't break in the edge case that the new point is past all existing tangent points"""
        fsq = pwfun.Pwlin (-5, 5, lambda x: -x*x, lambda x: -2*x)
        fsq.argmin (-4.5)
        self.check_knots(fsq)

    def check_knots (self, pwlin):
        x_tan = pwlin.tangent_points()
        ds = pwlin.section_derivs()
        hs = pwlin.section_heights()
        knots = pwlin.knots()
        
        for i in xrange(len(x_tan)-1):
            x1 = x_tan[i]
            x2 = x_tan[i+1]
            k = knots[i+1]
            z1 = hs[i] + ds[i]*(k - x1)
            z2 = hs[i+1] + ds[i+1]*(k - x2)
            self.assertTrue (x1 <= k, "WEIRDNESS: x1 %.4f k %.4f\n%s" % (x1, k, pwlin.dump()))
            self.assertTrue (k <= x2)
            self.assertAlmostEqual (z1, z2, 5)
            
    def test_pwfun_multiply (self):        
        fn1 = self.f_sq_log
        fn2 = self.f_abs
        f12 = pwfun.multiply (fn1, fn2)
        
        self.assertEquals (f12.knots(), [ -5.0, -1.0, 0.0, 1.0, 5.0 ])
        
        self.check_pwfun_product (fn1, fn2, f12, [ -2, -0.5, 0.5, 2.0 ])

    def test_pwfun_multiply2 (self):        
        f1 = self.f_sq_log
        f2 = pwfun.Pwfun ([-2.0, 0.0, 1.0], self.f_abs.fns, self.f_abs.derivs)
        f12 = pwfun.multiply (f1, f2)
        
        self.assertEquals (f12.knots(), [ -2.0, -1.0, 0.0, 1.0 ])
        
        self.check_pwfun_product (f1, f2, f12, [ -2, -1.5, -0.5, 0.5] )

        
    def test_pwfun_multiply3 (self):
        k1 = [ 3., 7., 9. ]
        k2 = [ 0., 1., 5. ]
        f1 = pwfun.Pwfun (k1, self.f_sq_log.fns[0:2], self.f_sq_log.derivs[0:2])
        f2 = pwfun.Pwfun (k2, self.f_sq_log.fns[0:2], self.f_sq_log.derivs[0:2])
        f12 = pwfun.multiply (f1, f2)
        self.assertEquals (f12.knots(), [ 3., 5. ])
        self.check_pwfun_product (f1, f2, f12, [3.0, 3.5, 4.0, 4.5 ])
        
    def check_pwfun_product (self, fn1, fn2, f12, points):
        for x in points:
            self.assertAlmostEquals (f12(x), fn1(x) * fn2(x), 5, "X: %.4f f1(x) = %.4f  f2(x) = %.4f  f12(x) = %.4f" % (x, fn1(x), fn2(x), f12(x)))
            delta = 1e-5
            dx_emp = (f12(x+delta) - f12(x)) / delta
            dx_analytic = f12.fprime (x)
            self.assertAlmostEquals (dx_emp, dx_analytic, 2)
        
    def test_inf_pwfun (self):
        fn = lambda x: -0.5*x*x
        dfn = lambda x: -x
        fexp = pwfun.Pwlin (0, numpy.inf, fn, dfn)
        print fexp.dump()
        self.assertEquals(3, len(fexp.knots()))
        
        pts = [0., 0.5, 1.0, 1.5, 2.0, 2.5]
        for x in pts:
            self.assertTrue (fn(x) <= fexp(x), "Upper bound doesn't hold: %s fn(%s) = %s  [real: %s]" % (x, x, fexp(x), fn(x)))

        fexp.argmin(3.0)
#        fexp.argmin(997.0)
        print fexp.dump()
        
        self.assertEquals(4, len(fexp.knots()))
        
        for x in pts:
            self.assertTrue (fn(x) <= fexp(x), "Upper bound doesn't hold: %s fn(%s) = %s  [real: %s]" % (x, x, fexp(x), fn(x)))
        
        for i in xrange(1000):
            u = numpy.random.uniform (0, 10)
            fexp.argmin(u)

        self.assertAlmostEquals (0.5*sqrt(2*numpy.pi), sum(fexp.integrate_exp_sections()), 4)


    def test_expinf_pwfun (self):
        gamm = distributions.Gamma(shape=2, scale=0.5)
        fexp = pwfun.Pwlin (0, numpy.inf, gamm.lpdf, gamm.dx_lpdf)
        print fexp.dump()
        self.assertEquals(3, len(fexp.knots()))

    def test_extreme_pwlin (self):
        exp = distributions.Exponential(scale=0.1)
        fexp = pwfun.Pwlin (0, 500, exp.lpdf, exp.dx_lpdf)
        cs = fexp.integrate_exp_sections()
        self.assertAlmostEquals (1.0, cs[0], 10)
        self.assertAlmostEquals (0.0, cs[1], 10)
        
    def test_extreme_sample (self):
        fn = lambda x: 5*x
        dfn = lambda x: 5.0
        fexp = pwfun.Pwlin (0, 1000, fn, dfn)
        print fexp.dump()
        
        k = fexp.knots()
        for i in xrange(100):
            x = fexp.sample_from_section_exp(1)
            self.assertTrue (x >= k[1], "x=%s not in %s..%s" % (x, k[1], k[2]))
            self.assertTrue (x <= k[2], "x=%s not in %s..%s" % (x, k[1], k[2]))

    def test_extreme_sample2 (self):
        fn = lambda x: -5*x
        dfn = lambda x: -5.0
        fexp = pwfun.Pwlin (0, 1000, fn, dfn)
        print fexp.dump()

        k = fexp.knots()
        for i in xrange(100):
            x = fexp.sample_from_section_exp(1)
            self.assertTrue (x >= k[1], "x=%s not in %s..%s" % (x, k[1], k[2]))
            self.assertTrue (x <= k[2], "x=%s not in %s..%s" % (x, k[1], k[2]))
            
    def test_sample_exp (self):
        fn = lambda x: -x
        dfn = lambda x: -1.0
        fexp = pwfun.Pwlin (0, 10, fn, dfn)
        fexp.argmin (1.0)
        
        knt = fexp.knots()
        k0 = knt[1]
        k1 = knt[2]
        
        N = 1000
        sampled = [ fexp.sample_from_section_exp (1) for i in xrange (N) ]
        
        # compare to rejection
        gold = []
        while len (gold) < N:
            x = numpy.random.exponential ()
            if (x > k0) and (x < k1): gold.append (x)
            
        sampled.sort()
        gold.sort()

        print "sampled: min %s med %s max %s" % (sampled[0], sampled[N/2], sampled[-1])
        print "gold:    min %s med %s max %s" % (gold[0], gold[N/2], gold[-1])
        
        netutils.check_quantiles (self, sampled, gold, N)

        
    def setUp (self):
        fn1 = lambda x: exp(-x)
        fn2 = lambda x: exp(-2*x)
        fn3 = lambda x: exp(-3*x)
        self.fexp1 = pwfun.Pwfun ([ 0, 1, 2, 3 ], [ fn1, fn2, fn3 ])
        
        xs = [-1, 0, 2, 5]
        derivs = [ 0.5, 0, -1.5 ]
        heights = [0, 2, -3]
        self.flin = pwfun.Pwlin (xs[0], xs[-1], lambda x: -x*x*x*x, lambda x: -4*x*x*x)

        self.fsq = pwfun.Pwlin (-5, 5, lambda x: -x*x, lambda x: -2*x)
        self.fsq.argmin (0)
        self.fsq.argmin (0.5)

        sq_knots = [-5, -1.0, 0.25, 0.75]
        sq_tan = [-2, 0.0, 0.5, 1] 
        sq_derivs = [4, -0.0, -1.0, -2]
        sq_heights = [-4, -0.0, -0.25, -1]
                
        fns = [ genexpfn(i, sq_tan, sq_derivs, sq_heights) for i in xrange(len(heights)) ]
        self.floglin = pwfun.Pwfun (sq_knots, fns)

        k1 = [ -5., -1., 1., 5. ]
        fns1 = [ lambda x: x*x, \
                 lambda x: 0.5*(x-1), \
                 lambda x: log(x) ]
        derivs1 = [ lambda x: 2*x, \
                    lambda x: 0.5, \
                    lambda x: 1/x ]
        self.f_sq_log = pwfun.Pwfun (k1, fns1, derivs1)
        
        k2 = [ -5., 0., 5.]
        fn2 = [ lambda x: -x, lambda x: x ]
        dx2 = [ lambda x: -1, lambda x: 1 ]
        self.f_abs = pwfun.Pwfun (k2, fn2, dx2)


def genfn(i, xs, derivs, heights): 
    return lambda x: heights[i] + derivs[i]*(x-xs[i])

def genexpfn(i, xs, derivs, heights): 
    return lambda x: exp(heights[i] + derivs[i]*(x-xs[i]))


if __name__ == "__main__":
    unittest.main()
#    test_name = "test_sample_exp"
#    suite = unittest.TestLoader().loadTestsFromName("test_pwfun.TestPwfun.%s" % (test_name,))
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
