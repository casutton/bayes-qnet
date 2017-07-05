# Simple linear programming wrapper. Not really general

# Completely stubbing everything out to reduce dependences
# that aren't needed to reproduce paper. July 2017

class LP:    pass

class QP:    pass

# from cvxopt.base import matrix, spmatrix
# from cvxopt import solvers, lapack
# import mytime

# solvers.options['maxiters'] = 500

# import pyglpk
# import os

# DEBUG = 0

# class LP:
    
#     def __init__ (self):
#         self.vars = []
#         self.constraints = []
#         self.fname = None
        
#     def add_vars (self, text):
#         self.vars.append (text)
        
#     def set_objective (self, text):
#         self.objective = text
        
#     def add_constraint (self, text):
#         self.constraints.append (text)
        
#     def solve (self):
#         self.fname = os.tempnam (os.getcwd(), "lp-model")
#         f = open (self.fname, 'w')
#         for var in self.vars:
#             f.write (var)
#             f.write ("\n")
#         f.write (self.objective)
#         f.write ("\n")
#         for ctr in self.constraints:
#             f.write (ctr)
#             f.write ("\n")
#         f.close()
        
#         print "Attempting to solve LP: ", self.fname
#         return pyglpk.solve_lp (self.fname)
        
#     def cleanup (self):
#         if self.fname:
#             os.remove (self.fname)
#             self.fname = None
        
# class QP:
    
#     def __init__ (self):
#         self.identities = dict()
#         self.ctrs_eq = []
#         self.b = []
#         self.ctrs = []
#         self.le_rhs = []
#         self.objsq = []
#         self.objlin = []
                       
#     def lookup_var (self, vname):
#         if vname in self.vars:
#             return self.vars[vname]
#         else:
#             vi = len(self.vars)
#             self.vars[vname] = vi
#             return vi
                        
#     def add_objsq (self, coeff, var):
#         self.objsq.append ( ( coeff, var, var ) )

#     def add_objsq_cross (self, coeff, v1, v2):
#         self.objsq.append ( ( coeff, v1, v2 ) )

#     def add_objlin (self, coeff, var):
#         self.objlin.append ( ( coeff, var ) )
        
#     def add_le0 (self, coeff, vars):
#         assert len(coeff) == len(vars)
#         eqn = ( coeff, vars[:] )
        
#         if eqn in self.ctrs:
#             raise Exception ("Equation %s already in QP: %s" % (eqn, self.ctrs))

#         self.ctrs.append (eqn)
#         self.le_rhs.append (0.0)

#     def add_var_le (self, var, value):
#         self.ctrs.append ( ([1.], [var]) )
#         self.le_rhs.append (value)

#     def add_var_ge (self, var, value):
#         self.ctrs.append ( ([-1.], [var]) )
#         self.le_rhs.append (-value)

#     def add_le (self, v1, v2):
#         self.ctrs.append ( ([1., -1.], [v1, v2]) )
#         self.le_rhs.append (0)
        
#     def add_eq (self, coeff, vars, rhs):
#         assert len(coeff) == len(vars)
#         self.ctrs_eq.append (( coeff, vars[:] ))
#         self.b.append (rhs)
    
#     def add_identity (self, var, value):
#         self.identities[var] = value
        
#     def collapse_identities (self):
#         self.collapse_identities_for (self.ctrs, self.le_rhs)
#         self.collapse_identities_for (self.ctrs_eq, self.b)
        
#         lin_new = dict()
#         sqnew = []
#         for c,v in self.objlin:
#             if not v in self.identities:
#                 lin_new[v] = c
                
#         for c, v1, v2 in self.objsq:
#             if (v1 in self.identities) and (v2 in self.identities):
#                 # drop it
#                 pass
#             elif (v1 in self.identities):
#                 c_old = lin_new.get (v2, 0)
#                 lin_new[v2] = c_old + c
#             elif (v2 in self.identities):
#                 c_old = lin_new.get (v1, 0)
#                 lin_new[v1] = c_old + c
#             else:
#                 # neither fixed
#                 sqnew.append ( (c,v1,v2) )

#         self.objsq = sqnew
#         self.objlin = [ (c,v) for v,c in lin_new.iteritems() ]
        
#         print self
            
#     def collapse_identities_for (self, all_lhs, all_rhs):
#         for i in xrange(len(all_lhs)):
#             w,lhs = all_lhs[i]
#             for j in xrange(len(lhs)):
#                 vi = lhs[j]
#                 if vi in self.identities:
#                     lhs[j] = None
#                     all_rhs[i] -= self.identities[vi]
#         j = 0
#         while j < len(all_lhs):
#             w,lhs = all_lhs[j]
#             if all(map(lambda x: x is None, lhs)):
#                 del all_lhs[j]
#                 del all_rhs[j]
#             j = j + 1

#     def construct_var_map (self):
#         self.vars = dict()
#         for w,lhs in self.ctrs:
#             for v in lhs:
#                 self.lookup_var (v)
#         for w,lhs in self.ctrs_eq:
#             for v in lhs:
#                 self.lookup_var (v)
#         for w,v1,v2 in self.objsq:
#             self.lookup_var(v1)
#             self.lookup_var(v2)
#         for w,v1 in self.objlin:
#             self.lookup_var(v1)
            
#     def construct_problem (self):
#         self.collapse_identities()
#         self.construct_var_map()
        
#         nv = len(self.vars)
#         neq = len(self.ctrs_eq)
#         n_ineq = len(self.ctrs)
        
#         print "QP: Num variables = %d  Equality constraints = %d  Inequality constraints %d" % (nv, neq, n_ineq)
        
#         # build P
#         x = []
#         I = []
#         J = []
#         for c,v1,v2 in self.objsq:
#             idx1 = self.vars[v1]
#             idx2 = self.vars[v2]
#             idx_min = min(idx1,idx2)
#             idx_max = max(idx1,idx2)
#             x.append (c)
#             I.append (idx_max)
#             J.append (idx_min)
#         P = spmatrix(x, I, J, (nv, nv))
            
#         q = matrix(0.0, (nv,1))
#         for c,v in self.objlin: 
#             idx = self.vars[v]
#             q[idx] = c
        
#         # build A (inequality constraints)
#         G = self.matrix_from_hash (self.ctrs, n_ineq, nv)
#         A = self.matrix_from_hash (self.ctrs_eq, neq, nv)

#         b = matrix (0.0, (A.size[0],1))
#         for i in xrange(len(self.b)):
#             b[i] = self.b[i]

#         h = matrix (0.0, (G.size[0],1))
#         for i in xrange(len(self.le_rhs)):
#             h[i] = self.le_rhs[i]
        
#         if DEBUG:
# #            S = matrix(0.0, (min(A.size),1))
# #            lapack.gesvd(A,S)
# #            r = 0
# #            for i in S: 
# #                if abs(i) > 1e-10:
# #                    r = r+1
# #            print "Rank(A) approx ",r

#             print "P: ", P
#             print "q: ", q
#             print "A: ", A
#             print "b: ", b
#             print "G: ", G            
            
# #        A = matrix(0.0, (0,nv))
# #        G = matrix(0.0, (0,nv))
# #        b = matrix(0.0, (0,1))

#         return P,q,A,b,G,h

#     def __repr__ (self): 
#         s = "min "
#         for c,n1,n2 in self.objsq:
#             s = s + " + (%.10f * %s * %s)" % (c, n1, n2) 
#         for c,n1 in self.objlin:
#             s = s + " + (%.10f * %s)" % (c, n1) 
#         s = s + "\nsubject to\n  "
#         s = s + "%d\n  " % len(self.ctrs)
#         for ((c_list,v_list), rhs) in zip(self.ctrs, self.le_rhs):            
#             for i in xrange(len(c_list)):
#                 v = v_list[i]
#                 c = c_list[i]
#                 s = s + " + (%.10f * %s)" % (c,v)
#             s = s + " <=  %s\n  " % rhs
#         ri = 0
#         for c_list,v_list in self.ctrs_eq:            
#             for i in xrange(len(c_list)):
#                 if v_list[i]:
#                     v = v_list[i]
#                     c = c_list[i]
#                     s = s + " + (%.10f * %s)" % (c,v)
#             b = self.b[ri]
#             s = s + " == %.10f" % b
#             s = s + "\n  "
#             ri = ri + 1
#         return s
        
#     def unpack_soln (self, vec):    
#         soln = dict()
#         for name,vi in self.vars.iteritems():
#             soln[name] = vec[vi]
#         return soln
        
#     def solve (self):
#         tmr = mytime.timeit()
        
#         P,q,A,b,G,h = self.construct_problem()
#         tmr.tick ("QP: Constructing problem")

#         soln = solvers.coneqp(P, q, G, h, None, A, b, None, kktsolver="chol2")
#         tmr.tick ("QP: solving problem")
        
#         if soln['status'] != 'optimal':
#             print "Error in QP: ", soln['status']
#         return self.unpack_soln (soln['x'])
        

#     def matrix_from_hash (self, hash, nr, nc):
#         x = []
#         I = []
#         J = []
    
#         ri = 0
#         for c_list,v_list in hash:
#             for i in xrange(len(c_list)):
#                 vname = v_list[i]
#                 if vname:
#                     vidx = self.vars [vname] 
#                     c = c_list[i]
#                     print "(%d, %s, %s)" % (ri, vidx, c)
#                     I.append (ri)
#                     J.append (vidx)
#                     x.append (c)
#             ri = ri + 1  
        
#         print "I:", I
#         print "J:", J
#         print (nr,nc)
#         A = spmatrix (x, I, J, (nr,nc))
#         return A
    
