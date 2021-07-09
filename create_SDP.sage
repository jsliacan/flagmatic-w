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

    elif problem == 'Thomason_K3K4':
        p = GraphProblem(N, density=[("3:",1),("4:121314232434", 1)], minimize=True)
        sense = "max"

    elif problem == 'Thomason_K4K4':
        p = GraphProblem(N, density=[("4:",1),("4:121314232434", 1)], minimize=True)
        sense = "max"

    I, J, M, C, c = p.get_sdp()
    write_problem(I, J, M, C, c, sense, "/home/sage/out/%s_%d" % (problem, N))