import unittest
import stupidlp
import cvxopt.solvers

class TestStupidLP (unittest.TestCase):  
    
    def test_qp1 (self):
        qp = stupidlp.QP()
        
        qp.add_objsq (1.0, "V1")
        qp.add_objsq (1.5, "V2")
        
        qp.add_objlin (3., "V1")
        qp.add_objlin (-4., "V2")
        
        soln = qp.solve()
        
        self.assertAlmostEquals (soln["V1"], -3., 5)
        self.assertAlmostEquals (soln["V2"], 2.6667, 3)

    def test_qp2 (self):
        qp = stupidlp.QP()

        # minimize variance of a 
        qp.add_objsq (-1.0, "P2")
        qp.add_objsq_cross (-4.0, "P2", "P3")
        qp.add_objlin (1., "P2")

        qp.add_le0 ([-1.0], ["P1"])
        qp.add_le0 ([-1.0], ["P2"])
        qp.add_le0 ([-1.0], ["P3"])
        
        qp.add_eq ([1.0, 1.0, 1.0 ], ["P1", "P2", "P3" ], 1.0)
        qp.add_eq ([0., 1., 2.], ["P1","P2","P3"], 1./3)
        
        print qp
        
#        cvxopt.solvers.options['debug'] = 1
        
        soln = qp.solve()
        
        self.assertAlmostEquals (soln['P1'], 0.83333327337120411, 5)
        self.assertAlmostEquals (soln['P2'], 1.19924258508573e-07, 5)
        self.assertAlmostEquals (soln['P3'], 0.1666666067045374, 5)


    def test_qp3 (self):
        qp = stupidlp.QP()

        # minimize variance of a 
        qp.add_objsq (-1.0, "P2")
        qp.add_objsq_cross (-4.0, "P2", "P3")
        qp.add_objlin (1., "P2")

        qp.add_le0 ([-1.0], ["P1"])
        qp.add_le0 ([-1.0], ["P2"])
        qp.add_le0 ([-1.0], ["P3"])

        qp.add_eq ([1.0, 1.0, 1.0 ], ["P1", "P2", "P3" ], 1.0)        
        qp.add_eq ([0., 1., 2.], ["P1","P2","P3"], 1./3)
        qp.add_eq ([1.], ["P2"], 0.1)
        soln = qp.solve()
        
        expected = {'P2': 0.10000000000000001, 'P3': 0.11666666666666667, 'P1': 0.78333333333333321}
        for vname in soln:
            self.assertAlmostEquals (expected[vname], soln[vname], 5)

    def test_qp4 (self):
        qp = stupidlp.QP()
        
        qp.add_objsq (1., "X")
        qp.add_objsq (1., "Y")
        
        qp.add_eq ([1., 1.], ["X", "Y"], 1.0)
        
        qp.add_var_le ("X", 0.2)

        soln = qp.solve()
        print soln

        expected = {'Y': 0.80000000642238134, 'X': 0.19999999357761869}
        for vname in soln:
            self.assertAlmostEquals (expected[vname], soln[vname], 5)

    # Below fails because CVXOPT doesn't like the fact that the equality constraint matrix is singular
    def ignore_test_qp5 (self):
        qp = stupidlp.QP()

        # minimize variance of a, constraining P2 == 0.1
        qp.add_objsq (-1.0, "P2")
        qp.add_objsq_cross (-4.0, "P2", "P3")
        qp.add_objlin (1., "P2")

        qp.add_le0 ([-1.0], ["P1"])
        qp.add_le0 ([-1.0], ["P2"])
        qp.add_le0 ([-1.0], ["P3"])

        qp.add_eq ([1.0, 1.0, 1.0 ], ["P1", "P2", "P3" ], 1.0)        
        qp.add_eq ([0., 1., 2.], ["P1","P2","P3"], 1./3)
        qp.add_identity ("P2", 0.1)
        
        soln = qp.solve()

        expected = {'P3': 0.11666666666666667, 'P1': 0.78333333333333321}
        for vname in soln:
            self.assertAlmostEquals (expected[vname], soln[vname], 5)


        
if __name__ == "__main__":
    unittest.main() 
#    test_name = "test_qp5"
#    suite = unittest.TestLoader().loadTestsFromName("test_stupidlp.TestStupidLP.%s" % (test_name,))
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
