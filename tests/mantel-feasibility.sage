from flagmatic.all import *

p = GraphProblem(3, mode="feasibility")
p.add_assumption("0:", [("3:121323(0)",-1)], 0)
p.add_assumption("0:", [("2:12(0)",1)], 1/2)
#p.add_assumption("0:", [("3:12(0)",1)], 1)
p.solve_sdp()
