from flagmatic.all import *
from sage.all import *

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

tol = 1e-10

if sys.argv[1] == 'c_4':
    N = 8
    density = [("4:", 1), ("4:121314232434", 1)]
    forbid = None
    C = None
    certificate = "c_4.js"
    reduce_types = False

elif sys.argv[1] == 'c_5':
    N = 8
    density = [("5:", 1), ("5:12131415232425343545", 1)]
    forbid = None
    C = None
    certificate = "c_5.js"
    reduce_types = False

elif sys.argv[1] == 'c_4,5':
    N = 8
    density = [("4:", 1), ("5:12131415232425343545", 1)]
    forbid = None
    C = None
    certificate = "c_45.js"
    reduce_types = False

elif sys.argv[1] == 'c_3,4':
    N = 6
    density = [("3:", 1), ("4:121314232434", 1)]
    forbid = None
    C = GraphBlowupConstruction(GraphFlag(graphs.SchlaefliGraph()))
    certificate = "c_34.js"
    reduce_types = True

elif sys.argv[1] == 'c_3,5':
    N = 6
    density = [("5:", 1), ("3:121323", 1)]
    forbid = None
    C = GraphBlowupConstruction(GraphFlag(graphs.SchlaefliGraph().complement()))
    certificate = "c_35.js"
    reduce_types = True

elif sys.argv[1] == 'g_4,5':
    N = 7
    density = "4:"
    forbid = "5:12131415232425343545"
    C = GraphBlowupConstruction(GraphFlag(Graph(r"LJ]lmZRnn]]\v[")), no_symmetry=True)
    certificate = "g_45.js"
    reduce_types = True

elif sys.argv[1] == 'g_6,3':
    N = 7
    density = ("6:", 1)
    forbid = "3:121323"
    C = GraphBlowupConstruction("g:12131415162728292a373b3c3d484b4e4f595c5e5g6a6d6f6g7e7f7g8c8d8g9b9d9fabacaebgcfde")
    certificate = "g_63.js"
    reduce_types = True
    
elif sys.argv[1] == 'g_4,6':
    N = 8
    density = "4:"
    forbid = "6:121314151623242526343536454656"
    C = None
    certificate = "g_46.js"
    reduce_types = False
    
elif sys.argv[1] == 'g_4,7':
    N = 8
    density = "4:"
    forbid = "7:121314151617232425262734353637454647565767"
    C = None
    certificate = "g_47.js"
    reduce_types = False

elif sys.argv[1] == 'g_5,4':
    N = 8
    density = "5:"
    forbid = "4:121314232434"
    C = None
    certificate = "g_55.js"
    reduce_types = False

elif sys.argv[1] == 'g_5,5':
    N = 8
    density = "5:"
    forbid = "5:12131415232425343545"
    C = None
    certificate = "g_55.js"
    reduce_types = False
    
elif sys.argv[1] == 'g_5,6':
    N = 8
    density = "5:"
    forbid = "6:121314151623242526343536454656"
    C = None
    certificate = "g_56.js"
    reduce_types = False
    
elif sys.argv[1] == 'g_5,7':
    N = 8
    density = "5:"
    forbid = "6:121314151617232425262734353637454647565767"
    C = None
    certificate = "g_57.js"
    reduce_types = False


if reduce_types:
    with suppress_stdout():
        P = GraphProblem(N, density=density, forbid_induced=forbid, minimize=True)
        # P.set_extremal_construction(C)
        P.solve_sdp(solver="sdpa_dd")
        bound = eval(P._sdp_solver_output.split('\n')[-10].split(' = ')[1])
        # bound = eval(P._sdp_solver_output.split("\n")[-8].split(" ")[-1]) #29/13**3
        cand_types = P._types

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
            P = GraphProblem(N, density=density, forbid_induced=forbid, types=next_types, minimize=True)
            P.solve_sdp(solver="sdpa_ddx")
            sol = eval(P._sdp_solver_output.split('\n')[-10].split(' = ')[1])
            # sol = eval(P._sdp_solver_output.split("\n")[-8].split(" ")[-1])

        if sol > bound-tol:
            removed_types.append(next_type)
            tested_types = []

            print("removed type "+str(next_type)+" down to "+str(len(current_types)-1)+" and getting "+str(sol)+" instead of "+str(bound))

        else:
            tested_types.append(next_type)
            print("couldn't remove "+str(next_type)+" getting "+str(sol)+" instead of "+str(bound))

    current_types = [t for t in cand_types if t not in removed_types]

    with suppress_stdout():
        P = GraphProblem(N, density=density, forbid_induced=forbid, types=current_types, minimize=True)
        if C is not None: P.set_extremal_construction(C)
        P.solve_sdp(solver="sdpa_qd")
        reduced_bound = eval(P._sdp_solver_output.split('\n')[-10].split(' = ')[1])
        # reduced_bound = eval(P._sdp_solver_output.split("\n")[-8].split(" ")[-1])

    print("finished with "+str(len(current_types))+" out of "+str(orig_n_types)+" getting "+str(reduced_bound)+" instead of "+str(bound)+":")
    for t in current_types:
        print(t)

    P.make_exact(2000)
    P.write_certificate(certificate)
    
else:
    P = GraphProblem(N, density=density, forbid_induced=forbid, minimize=True)
    if C is not None: P.set_extremal_construction(C)
    P.solve_sdp(solver="csdp")
    P.make_exact(2^30)
    P.write_certificate(certificate)