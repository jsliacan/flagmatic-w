from flagmatic.all import *

p = GraphProblem(5, forbid="3:121323", density="5:1223344551")
p.set_extremal_construction(GraphBlowupConstruction("5:1223344551"))
p.solve_sdp()
p.make_exact()
p.verify_robust_stability("3:12")
p.verify_perfect_stability()