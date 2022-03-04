from flagmatic.all import *

def tuple_to_string(data):
    return ','.join([str(x) for x in data])

def write_problem(I, J, M, C, c, sense, problem_name):
    
     with open("/home/sage/"+problem_name, 'w') as f:
        f.write(str(I))
        f.write(" ")
        f.write(str(J))
        f.write(" ")
        f.write(sense)
        f.write("\n")

        f.write(" ".join([str(m) for m in M]))
        f.write("\n")

        f.write(" ".join([str(x) for x in c]))
        f.write("\n")

        for i in range(I):
            for j in range(J):
                f.write(str(i))
                f.write(" ")
                f.write(str(j))
                f.write(" ")
                f.write(" ".join([tuple_to_string(x) for x in C[i][j]]))
                f.write("\n")


if __name__=='__main__':
    import sys

    _, problem, N = sys.argv

    N = int(N)

    if problem == 'Mantel':
        p = GraphProblem(N, forbid="3:121323")
        sense = "min"

    elif problem == 'Turan_C5':
        p = GraphProblem(N, forbid=(3,3), density="5:1223344551")
        sense = "min"

    elif problem == 'Turan_Conjecture':
        p = ThreeGraphProblem(N, forbid=(4,4))
        sense = "min"

    elif problem == 'Thomason_K2K3':
        p = GraphProblem(N, density=[("2:", 4/7),("3:121323", 3/7)], minimize=True)
        sense = "max"
        
    elif problem == 'Thomason_K2K4':
        p = GraphProblem(N, density=[("2:", 1),("4:121314232434", 1)], minimize=True)
        sense = "max"

    elif problem == 'AsymThomason_K2K3':
        p = GraphProblem(N, density=[("2:", 1),("3:121323", 1)], minimize=True)
        sense = "max"
                                                                                                                       
    elif problem == 'Thomason_K3K4':
        p = GraphProblem(N, density=[("3:", 1),("4:121314232434", 1)], minimize=True)
        sense = "max"
                                                                                                                       
    elif problem == 'Thomason_K3K4_25s26a27s26':
        p = GraphProblem(N, density=[("3:", 25/26),("4:121314232434", 27/26)], minimize=True)
        sense = "max"
                                                                                                                       
    elif problem == 'Thomason_K3K4_27s26a25s26':
        p = GraphProblem(N, density=[("3:", 27/26),("4:121314232434", 25/26)], minimize=True)
        sense = "max"

    elif problem == 'Thomason_K3K5':
        p = GraphProblem(N, density=[("3:", 1),("5:12131415232425343545", 1)], minimize=True)
        sense = "max"

    elif problem == 'Thomason_K3K6':
        p = GraphProblem(N, density=[("3:", 1),("6:121314151623242526343536454656", 1)], minimize=True)
        sense = "max"

    elif problem == 'Thomason_K3K7':
        p = GraphProblem(N, density=[("3:", 1),("7:121314151617232425262734353637454647565767", 1)], minimize=True)
        sense = "max"

    elif problem == 'Thomason_K3K8':
        p = GraphProblem(N, density=[("3:", 1),("8:12131415161718232425262728343536373845464748565758676878", 1)], minimize=True)
        sense = "max"
        
    elif problem == 'Thomason_K4K5':
        p = GraphProblem(N, density=[("4:", 1),("5:12131415232425343545", 1)], minimize=True)
        sense = "max"

    elif problem == 'Thomason_K3K3':
        p = GraphProblem(N, density=[("3:", 1),("3:121323", 1)], minimize=True)
        sense = "max"

    elif problem == 'Thomason_K4K4':
        p = GraphProblem(N, density=[("4:", 1),("4:121314232434", 1)], minimize=True)
        sense = "max"

    elif problem == 'Thomason_K3K4_axis':
        p = GraphProblem(N, density=[("4:", 1)], forbid="3:121323", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K4K3_axis':
        p = GraphProblem(N, density=[("3:", 1)], forbid="4:121314232434", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K4K5_axis':
        p = GraphProblem(N, density=[("5:", 1)], forbid="4:121314232434", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K5K4_axis':
        p = GraphProblem(N, density=[("4:", 1)], forbid="5:12131415232425343545", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K6K4_axis':
        p = GraphProblem(N, density=[("4:", 1)], forbid="6:121314151623242526343536454656", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K7K4_axis':
        p = GraphProblem(N, density=[("4:", 1)], forbid="7:121314151617232425262734353637454647565767", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K6K5_axis':
        p = GraphProblem(N, density=[("5:", 1)], forbid="6:121314151623242526343536454656", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K6K6_axis':
        p = GraphProblem(N, density=[("5:", 1)], forbid="6:121314151623242526343536454656", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K7K5_axis':
        p = GraphProblem(N, density=[("5:", 1)], forbid="7:121314151617232425262734353637454647565767", minimize=True)
        sense = "max"
        
    elif problem == 'Thomason_K7K7_axis':
        p = GraphProblem(N, density=[("7:", 1)], forbid="7:121314151617232425262734353637454647565767", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K4K4_axis':
        p = GraphProblem(N, density=[("4:", 1),("4:121314232434", 1)], forbid="4:121314232434", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K5K5_axis':
        p = GraphProblem(N, density=[("5:", 1),("5:12131415232425343545", 1)], forbid="5:12131415232425343545", minimize=True)
        sense = "max"

    elif problem == 'Thomason_K3eK3e':
        p = GraphProblem(N, density=[("4:1223", 1),("4:12132314", 1)], minimize=True)
        sense = "max"

    elif problem == 'Thomason_K5K5':
        p = GraphProblem(N, density=[("5:", 1),("5:12131415232425343545", 1)], minimize=True)
        sense = "max"

    I, J, M, C, c = p.get_sdp()
    write_problem(I, J, M, C, c, sense, "%s_%d" % (problem, N))
