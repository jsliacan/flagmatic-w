from flagmatic.all import *

p = GraphProblem(3, forbid="3:121323")
c = GraphBlowupConstruction("2:12")
p.set_extremal_construction(c)
p.solve_sdp()
p.make_exact()

assert(p._exact_Q_matrices[0] == matrix([[1/2,-1/2],[-1/2,1/2]]))