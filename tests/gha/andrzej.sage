from flagmatic.all import *

p = GraphProblem(5, forbid="3:121323", density="5:1223344551")
c = GraphBlowupConstruction("5:1223344551")
p.set_extremal_construction(c)
p.solve_sdp()
p.make_exact()

assert(p._bound == 24/625)