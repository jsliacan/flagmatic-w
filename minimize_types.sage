from flagmatic.all import *
import random


from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

N = 6
density = [("3:", 1), ("5:12131415232425343545", 1)] #[("4:", 1)]
forbid = None #"5:12131415232425343545"
label = r"Z??G`@?@wrDSLGQoigbKO]CA?^{VDsjIqehgmK[EM[OzIqCyegO|FO_^{?_?"
            
with suppress_stdout():
    P = GraphProblem(N, density=density, forbid=forbid, minimize=True)
    P.solve_sdp()
    bound = eval(P._sdp_solver_output.split("\n")[-8].split(" ")[-1]) #29/13**3
    cand_types = P._types
    
tol = 1e-10

removed_types = []
tested_types = []

orig_n_types = len(cand_types)

print("Starting with "+str(orig_n_types)+" types")

while True:
    current_types = [t for t in cand_types if t not in removed_types]
    types_to_consider = [t for t in current_types if t not in tested_types]
    
    if len(types_to_consider) == 0: break
    
    next_type = random.choice(types_to_consider)
    next_types = [t for t in current_types if t != next_type]
    
    with suppress_stdout():
        P = GraphProblem(N, density=density, forbid=forbid, types=next_types, minimize=True)
        P.solve_sdp()
        sol = eval(P._sdp_solver_output.split("\n")[-8].split(" ")[-1])

    if sol > bound-tol:
        removed_types.append(next_type)
        tested_types = []
        
        print("removed type "+str(next_type)+" down to "+str(len(current_types)-1))
        
    else:
        tested_types.append(next_type)
        print("couldn't remove "+str(next_type))
        
current_types = [t for t in cand_types if t not in removed_types]

with suppress_stdout():
    C = GraphBlowupConstruction(GraphFlag(Graph(label).complement()))
    P = GraphProblem(N, density=density, forbid=forbid, types=current_types, minimize=True)
    P.set_extremal_construction(C)
    P.solve_sdp()
    reduced_bound = eval(P._sdp_solver_output.split("\n")[-8].split(" ")[-1])

print("finish with "+str(len(current_types))+" out of "+str(orig_n_types)+" getting "+str(reduced_bound)+" instead of "+str(bound)+":")
for t in current_types:
    print(t)
    

        
P.make_exact_plus(limit_denominator=100)
P.write_certificate("c_3,5_compact.cert")