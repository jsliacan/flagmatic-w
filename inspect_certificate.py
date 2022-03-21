"""

Copyright (c) 2011, E. R. Vaughan. All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import fractions
import getopt
import itertools
import json
import re
import sys
import numpy
import subprocess

HELP_TEXT = """
Usage: inspect_certificate.py CERTIFICATE [BOUND TAU_GRAPH B_GRAPH CERTIFICATE_TAU [CERTIFICATE_B]]  OPTIONS
Possible options:
--help                       Display this message.
--admissible-graphs          Display the admissible graphs.
--flags                      Display the types and flags.
--r-matrices                 Display the R matrices.
--qdash-matrices             Display the Q' matrices.
--q-matrices                 Compute and display the Q matrices.
--pair-densities             Compute and display the flag pair densities.
--verify-bound               Verify the bound.
--sharp-graphs               Display the admissible graphs that are sharp.
--flag-algebra-coefficients  Display each admissible graph's flag algebra coefficient.
--stability                  Verify stability as well.
--log                        Log the bound.
"""

try:
    from sage.all import *
    using_sage = True
except ImportError:
    using_sage = False

should_print_help = False

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "", ["help", "admissible-graphs",
                                                      "flags", "r-matrices", "qdash-matrices",
                                                      "pair-densities", "q-matrices", "verify-bound",
                                                      "sharp-graphs", "flag-algebra-coefficients", "stability", "log"])

except getopt.GetoptError:
    should_print_help = True
    opts, args = ((), ())

if len(args) < 1:
    should_print_help = True

action = ""

for o, a in opts:
    if o == "--help":
        should_print_help = True
    elif o == "--admissible-graphs":
        action = "print admissible graphs"
    elif o == "--flags":
        action = "print flags"
    elif o == "--r-matrices":
        action = "print r_matrices"
    elif o == "--qdash-matrices":
        action = "print qdash_matrices"
    elif o == "--pair-densities":
        action = "print pair densities"
    elif o == "--q-matrices":
        action = "print q_matrices"
    elif o == "--verify-bound":
        action = "verify bound"
    elif o == "--sharp-graphs":
        action = "print sharp graphs"
    elif o == "--flag-algebra-coefficients":
        action = "print flag algebra coefficients"
    elif o == "--stability":
        action = "verify stability"
    elif o == "--log":
        action = "log"
        
if should_print_help:
    print HELP_TEXT
    sys.exit(0)

certificate_filename = args[0]
lower_bound = ""
tau = ""
B = ""
certificate_filename_tau = ""

if action == "verify stability":
    if using_sage == False:
        print("This step needs Sage. Please run 'sage -python inspect_certificate.py ...'.")
        sys.exit()
        
    if len(args) < 4:
        raise ValueError("If '--stability' option chosen, need at least 4 arguments: certificate, bound, tau, B.\n")

    # target bound
    lb = map(int, args[1].split("/"))
    lower_bound = Integer(lb[0])/lb[1] #fractions.Fraction(lb[0],lb[1])
    print("Lower bound:", lower_bound)
    # tau graph
    tau = str(args[2])
    print("Type tau:", tau)
    # B graph
    B = str(args[3])
    print("Graph B:", B)

    if tau != "1:" and len(args) < 5:
        raise ValueError("With '--stability' option and tau is not '1:', second certificate is needed.\n")
    # cert_tau
    if tau != "1:":
        certificate_filename_tau = str(args[4])
        print("Certificate when tau is forbidden:", certificate_filename_tau)

    
try:
    if certificate_filename[-3:] == ".gz":
        import gzip
        certf = gzip.open(certificate_filename)
    elif certificate_filename[-4:] == ".bz2":
        import bz2
        certf = bz2.BZ2File(certificate_filename)
    else:
        certf = open(certificate_filename)
except IOError:
    sys.stdout.write("Could not open certificate.\n")
    sys.exit(1)
print

certificate = json.load(certf)

print('Problem is "{}".'.format(certificate["description"]))
print('Claimed bound is {}.'.format(certificate["bound"]))

minimize = "minimize" in certificate["description"]

if using_sage:
    x = polygen(QQ)
else:
    x = None

if "field" in certificate.keys():

    if certificate["field"] != "RationalField()":

        print('Field used is "{}".'.format(certificate["field"]))

        if not using_sage:
            print("This script must be run using Sage's python, as it does not use the rational field.")
            sys.exit(1)

        try:
            K = sage_eval(certificate["field"], locals={'x': x})
        except AttributeError:
            K = RationalField()
        x = K.gen()

if action == "":
    sys.exit(0)

if "3-graph" in certificate["description"]:
    edge_size = 3
elif "2-graph" in certificate["description"]:
    edge_size = 2
else:
    print("Unsupported graph kind.")
    sys.exit(1)

oriented = "oriented" in certificate["description"]

def make_number(ch):
    try:
        CH = ch.capitalize()
        return "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ".index(CH)
    except ValueError:
        print("Failed to parse...")


def kill_twins(g):
    """
    Identifies all twins in g and returns the resulting graph.

    INPUT:
    g - graph to be de-twinned; in SAGE format
    """
    twins = None
    twin_free = False
    while not twin_free:
        found_twins = False
        twin = None
        gvertices = g.vertices()
        for i in range(g.num_verts()-1):
            for j in range(i+1,g.num_verts()):
                if g.vertex_boundary({gvertices[i]}) == g.vertex_boundary({gvertices[j]}):
                    twin = vertices[i]
                    found_twins = True
                    break
            if found_twins:
                break
        if found_twins:
            g.delete_vertex(twin)
        else:
            twin_free = True
    return g

    
class Flag(object):

    edge_size = edge_size
    oriented = oriented

    __slots__ = ("n", "t", "edges")

    def __init__(self, s=None):

        if not s:
            self.n = 0
            self.t = 0
            self.edges = ()
            return

        m = re.match(r'(\d+):(\d*)(?:\((\d+)\)|)', s)
        edges_str = ""
        if m:
            n = int(m.group(1))
            t = int(m.group(3)) if m.group(3) else 0
            if not (0 <= n <= 9 and 0 <= t <= n):
                raise ValueError

        elif s[1] == ':': # only get here if parsing above failed
            n = int(make_number(s[0]))
            t = 0
            if s[-1] == ')':
                t = int(make_number(s[-2]))
                edges_str = s[2:-3]
            else:
                edges_str = s[2:]
            
        else:
            raise ValueError("Wrong representation.")
        

            
        edges = []
        if m:
            for i in range(0, len(m.group(2)), self.edge_size):
                edges.append(tuple(map(int, m.group(2)[i:i + self.edge_size])))
        else:
            for i in range(0, len(edges_str), self.edge_size):
                edges.append(tuple(map(make_number, edges_str[i:i+self.edge_size])))
                
        self.n = n
        self.t = t
        self.edges = tuple(edges)

    def __repr__(self):
        s = "{}:{}".format(self.n, "".join(str(x) for e in self.edges for x in e))
        if self.t > 0:
            s += "({})".format(self.t)
        return s

    def __eq__(self, other):
        return self.n == other.n and self.t == other.t and self.edges == other.edges

    def __ne__(self, other):
        return not self.__eq__(other)

    def induced_subgraph(self, S):
        if not all(0 <= x <= self.n for x in S):
            raise ValueError
        good_edges = [e for e in self.edges if all(x in S for x in e)]
        p = [0 for _ in range(self.n + 1)]
        for i, x in enumerate(S):
            p[x] = i + 1
        h = Flag()
        h.n = len(S)
        h.t = 0
        if self.oriented:
            h.edges = tuple(sorted([tuple([p[x] for x in e]) for e in good_edges]))
        else:
            h.edges = tuple(sorted([tuple(sorted([p[x] for x in e])) for e in good_edges]))
        return h

    def minimal_isomorph(self):
        min_edges = self.edges
        for p in itertools.permutations(range(self.t + 1, self.n + 1), self.n - self.t):
            pt = tuple(range(1, self.t + 1)) + p
            if self.oriented:
                edges = tuple(sorted([tuple([pt[x - 1] for x in e]) for e in self.edges]))
            else:
                edges = tuple(sorted([tuple(sorted([pt[x - 1] for x in e])) for e in self.edges]))
            if edges < min_edges:
                min_edges = edges
        h = Flag()
        h.n = self.n
        h.t = self.t
        h.edges = min_edges
        return h

    def induced_flag(self, tv, ov):
        type_vertices = list(tv)
        other_vertices = list(ov)
        edges = []
        for e in self.edges:
            if any(not (x in tv or x in ov) for x in e):
                continue
            edge = []
            for x in e:
                if x in tv:
                    edge.append(1 + type_vertices.index(x))
                else:
                    edge.append(len(tv) + 1 + other_vertices.index(x))
            edges.append(tuple(edge))
        h = Flag()
        h.n = len(tv) + len(ov)
        h.t = len(tv)
        h.edges = tuple(edges)
        return h

    def delete_vertex(self, vi):
        # delete vi without normalizing edges
        if self.t > 0:
            raise ValueError("Can only delete vertices of unlabelled graphs.\n")
        if not 1 <= vi <= n:
            raise ValueError("Vertex index (1,...,n) must be in that range.\n")

        ledges = list(self.edges)
        newledges = list()
        for e in ledges:
            if vi in e:
                continue
            else:
                newledges.append(e)

        self.n = self.n-1
        self.edges = tuple(newledges)


                
admissible_graphs = [Flag(s) for s in certificate["admissible_graphs"]]
types = [Flag(s) for s in certificate["types"]]
flags = [[Flag(s) for s in f] for f in certificate["flags"]]

if action == "print admissible graphs":

    print("There are {} admissible graphs:".format(len(admissible_graphs)))
    for i, g in enumerate(admissible_graphs):
        print("{}. {}".format(i + 1, g))
    sys.exit(0)

if action == "print flags":
    for ti, _type in enumerate(types):
        nf = len(flags[ti])
        print("Type {} ({}) has {} flags:".format(ti + 1, _type, nf))
        for i in range(nf):
            print("   {}. {}".format(i + 1, flags[ti][i]))
    sys.exit(0)


def stringify(a):
    if a is None:
        return "Not used."
    if isinstance(a, list):
        return map(stringify, a)
    s = str(a)
    try:
        n = int(s)
        return n
    except ValueError:
        return s


if action == "print r_matrices":

    for ti, _type in enumerate(types):
        print("R matrix for type {} ({}):".format(ti + 1, _type))
        print("{}".format(stringify(certificate["r_matrices"][ti])))
    sys.exit(0)

if action == "print qdash_matrices":

    for ti, _type in enumerate(types):
        print("Q' matrix for type {} ({}):".format(ti + 1, _type))
        print("{}".format(stringify(certificate["qdash_matrices"][ti])))
    sys.exit(0)


print("Computing Q matrices...")

Qs = []

for ti, _type in enumerate(types):

    if certificate["qdash_matrices"][ti] is None:
        Qs.append(None)
        continue

    if using_sage:
        QD = [[sage_eval(str(s), locals={'x': x}) for s in row] for row in certificate["qdash_matrices"][ti]]
    else:
        QD = [[fractions.Fraction(s) for s in row] for row in certificate["qdash_matrices"][ti]]

    if certificate["r_matrices"][ti] is None:
        Qs.append(QD)
        continue

    if using_sage:
        R = [[sage_eval(str(s), locals={'x': x}) for s in row] for row in certificate["r_matrices"][ti]]
    else:
        R = [[fractions.Fraction(s) for s in row] for row in certificate["r_matrices"][ti]]

    nq = len(QD)
    nf = len(flags[ti])

    Q = [[0 for j in range(i, nf)] for i in range(nf)]
    for l in range(nq):
        for k in range(l, nq):
            qlk = QD[l][k - l]
            if qlk != 0:
                for i in range(nf):
                    for j in range(i, nf):
                        Q[i][j - i] += R[i][l] * R[j][k] * qlk
                        if k != l:
                            Q[i][j - i] += R[i][k] * R[j][l] * qlk
    Qs.append(Q)


if action == "print q_matrices":

    for ti, _type in enumerate(types):
        print("Q matrix for type {} ({}):".format(ti + 1, _type))
        print("{}".format(stringify(Qs[ti])))
    sys.exit(0)

print("Computing pair densities...")

pair_densities = {}
n = certificate["order_of_admissible_graphs"]

for ti, _type in enumerate(types):
    nf = len(flags[ti])
    s = _type.n
    if len(set(flag.n for flag in flags[ti])) != 1:
        raise ValueError("Flags for a given type must all have the same order.")
    m = flags[ti][0].n
    for g in admissible_graphs:
        pairs_found = [[0 for j in range(k, nf)] for k in range(nf)]
        total = 0
        for p in itertools.permutations(range(1, n + 1)):
            if p[s] > p[m]:
                continue
            is_good_permutation = True
            for c in range(s, n):
                if c in (s, m, 2 * m - s):
                    continue
                if p[c - 1] > p[c]:
                    is_good_permutation = False
                    break
            if not is_good_permutation:
                continue
            total += 1

            it = g.induced_subgraph(p[:s])
            if it != _type:
                continue

            f1 = g.induced_flag(p[:s], p[s: m]).minimal_isomorph()
            f2 = g.induced_flag(p[:s], p[m: 2 * m - s]).minimal_isomorph()

            if1 = flags[ti].index(f1)
            if2 = flags[ti].index(f2)

            if if1 <= if2:
                pairs_found[if1][if2 - if1] += 1
            else:
                pairs_found[if2][if1 - if2] += 1

        for k in range(nf):
            for j in range(k, nf):
                pf = pairs_found[k][j - k]
                if pf > 0:
                    if j == k:
                        pf *= 2
                    if using_sage:
                        pair_densities[(ti, admissible_graphs.index(g), k, j)] = Integer(pf) / (total * 2)
                    else:
                        pair_densities[(ti, admissible_graphs.index(g), k, j)] = fractions.Fraction(pf, total * 2)

if action == "print pair densities":

    for i, g in enumerate(admissible_graphs):
        print("Pair densities for admissible graph {} ({}):".format(i + 1, g))

        for ti, _type in enumerate(types):
            print("   Non-zero densities for type {} ({}):".format(ti + 1, _type))

            for key, d in pair_densities.items():
                if key[:2] == (ti, i):
                    print("      Flags {} and {} ({} and {}): {}".format(key[2] + 1, key[3] + 1,)
                                                                         flags[ti][key[2]], flags[ti][key[3]], d)
    sys.exit(0)

print("Computing bound...")

if certificate["admissible_graph_densities"] and isinstance(certificate["admissible_graph_densities"][0], list):
    graph_densities = certificate["admissible_graph_densities"]
    density_coefficients = certificate["density_coefficients"]
else:
    graph_densities = [certificate["admissible_graph_densities"]]
    density_coefficients = [1]

bounds = [0 for _ in certificate["admissible_graphs"]]
for di, sdc in enumerate(density_coefficients):
    if using_sage:
        dc = sage_eval(str(sdc), locals={'x': x})
    else:
        dc = fractions.Fraction(sdc)
    for i, s in enumerate(graph_densities[di]):
        if using_sage:
            bounds[i] += dc * sage_eval(str(s), locals={'x': x})
        else:
            bounds[i] += dc * fractions.Fraction(s)

for key in pair_densities.keys():
    t, i, j, k = key
    if Qs[t] is None:  # check that type is used
        continue
    d = pair_densities[key]
    v = Qs[t][j][k - j]
    if minimize:
        v *= -1
    if j == k:
        bounds[i] += d * v
    else:
        bounds[i] += d * v * 2

bound = min(bounds) if minimize else max(bounds)

print("Bound is {}.".format(bound))

if action == "print sharp graphs":

    print("Sharp graphs are:")
    for i, g in enumerate(admissible_graphs):
        if bounds[i] == bound:
            print("{}. {}".format(i + 1, g))
    sys.exit(0)

if action == "print flag algebra coefficients":

    print("There are {} admissible graphs:".format(len(admissible_graphs)))
    for i, g in enumerate(admissible_graphs):
        print("{}. ({}) : {}".format(i + 1, g, bounds[i]))
    sys.exit(0)




if action == "log":
    try:
        f = open("bound.txt", 'w')
        f.write(str(bound))
        print("Bound logged to 'bound.txt'.")
        f.close()
    except IOError:
        print("Couldn't open 'bound.txt' file for writing.")
        sys.exit()
        


if action == "verify stability":

    print("\nChecking conditions for robust stability:")
    
    # Check that each forbidden graph is twin-free
    # --------------------------------------------
    claim1 = False # assume to begin with
    
    # read forbidden graphs from certificate (in description)
    forbidden_graphs = list()
    begin = False
    descr = certificate["description"]
    if "forbid" in descr:
        descr_list = descr.strip().split(" ")
        for i in range(len(descr_list)):
            item = descr_list[i]
            if item == "forbid":
                begin = True
                continue
            if begin == True:
                if item[-1] == ";" or item[-1] == ",":
                    item = item[:-1]
                if item[1] == ":":
                    forbidden_graphs.append(Flag(item))


    #check twin-freeness
    twins_exist = False
    twins_graph = None # forbidden graph B that has twins
    twins_vx = None # 2-tuple of x,y twins in B
    for g in forbidden_graphs:
        
        if twins_exist: # only continue if no twins found in previous graphs
            break
        
        vg = range(1,g.n+1) # vertices of g
        eg = g.edges # tuple of tuples
        
        for x in vg:
            
            if twins_exist: # only continue if no twins found for previous x vertices
                break
            
            Nx = list()
            for (a,b) in eg:
                if a == x:
                    Nx.append(b)
                elif b == x:
                    Nx.append(a)
                    
            for y in range(x+1,g.n+1):
                Ny = list()
                if not ((x,y) in eg or (y,x) in eg):
                    for (a,b) in eg or (b,a) in eg:
                        if a == y:
                            Ny.append(b)
                        elif b == x:
                            Ny.append(b)
                            
                    # compare with Nx
                    if set(Nx) & set(Ny):
                        twins_exist = True
                        twins_graph = g
                        twins_vx = (x,y)
                        break

        
    if twins_exist:
        print("\033[31m[FAIL] \033[mAt least one of the forbidden graphs is NOT twin-free.")
        print("       Witness:", str(twins_graph)+". Consider vertices", twins_vx[0], "and", str(twins_vx[1])+ ".\n")
    else:
        claim1 = True
        print("\033[32m[OK]   \033[mAll forbidden graphs are twin-free.")


    # Check that the upper bound is tight
    # -----------------------------------
    claim2 = False
    if bound == lower_bound:
        claim2 = True
        print("\033[32m[OK]   \033[mLower bound and upper bound match", str(bound)+".")
    else:
        print("\033[31m[FAIL] \033[mLower bound is", lower_bound, "and upper bound is", str(bound)+".")

    
    # Check that |tau| <= N-2
    # -----------------------
    tgraph = Flag(tau)
    claim3 = (tgraph.n <= admissible_graphs[0].n)

    if claim3:
        print("\033[32m[OK]   \033[mTau is on", tgraph.n, "vertices, which is at most", admissible_graphs[0].n-2,"=",admissible_graphs[0].n,"- 2.")
    else:
        print("\033[31m[FAIL] \033[mThe order of tau is wrong:", str(tgraph.n)+". It must be at most", str(admissible_graphs[0].n)+".")



    # Check that all sharp graphs have a strong hom into the construction graph B
    # (strong hom = preserves edges & nonedges; no need injective)
    # ---------------------------------------------------------------------------
    claim4 = False
    
    # sharp graphs
    sharp_graphs = list()
    for i, g in enumerate(admissible_graphs):
        if bounds[i] == bound:
            sharp_graphs.append(g)

    # convert sharp graphs to Sage form
    sharp_graphs_sage = list()
    for sg in sharp_graphs:
        g_str = sg.__repr__()
        g = Graph(sg.n)
        for i in range(2,len(g_str),2):
            g.add_edge((int(g_str[i])-1, int(g_str[i+1])-1))
        sharp_graphs_sage.append(g)

    # IDEA: contract twins until no more twins; then check if subgraph of B
    if B == "g:12131415162728292a373b3c3d484b4e4f595c5e5g6a6d6f6g7e7f7g8c8d8g9b9d9fabacaebgcfde": # comparing literally to save time!!!
        C = graphs.ClebschGraph() # sage graph exists in Sage, they have tools
    else:
        C = Graph(make_number(B[0]))
        for i in range(2,len(B),2):
            C.add_edge((make_number(B[i])-1, make_number(B[i+1])-1))
            
            
    # deal with same neighbourhoods
    failed = list()
    num_graphs_left = len(sharp_graphs_sage)
    for sg in sharp_graphs_sage:
        twins = None
        twin_free = False
        while not twin_free:
            found_twins = False
            twin = None
            gvertices = sg.vertices()
            for i in range(sg.num_verts()-1):
                for j in range(i+1,sg.num_verts()):
                    if sg.vertex_boundary({gvertices[i]}) == sg.vertex_boundary({gvertices[j]}):
                        twin = gvertices[i]
                        found_twins = True
                        break
                if found_twins:
                    break
            if found_twins:
                sg.delete_vertex(twin)
            else:
                twin_free = True

        # Check if this contracted sharp graph has induced+injective hom into B
        h = C.subgraph_search(sg, induced=True)
        if h == None:
            failed.append(sg)
            num_graphs_left -= 1

    if not failed:
        claim4 = True
        print("\033[32m[OK]   \033[mAll sharp graphs admit strong hom into B.")
    else:
        print("\033[31m[FAIL] \033[mNOT all sharp graphs admit strong homomorphism into B.")
        print("       e.g. no strong hom from", failed[0].edges(labels=None), "into B", "("+B+").")



    # There is exactly one strong homomorphism from tau to B (up to automorph of tau)
    # -------------------------------------------------------------------------------
    claim5 = False
    strong_hom = None # will store the unique strong hom if found

    # !!!
    # treat CLEBSCH graph separately, too big and Sage has tools
    Tg = Flag(tau)
    Tg = Tg.minimal_isomorph()
    Tau = Flag("6:1223344551")
    Tau = Tau.minimal_isomorph()
    if B == "g:12131415162728292a373b3c3d484b4e4f595c5e5g6a6d6f6g7e7f7g8c8d8g9b9d9fabacaebgcfde" and Tg == Tau: # comparing B literally!!!
        Bgraph = graphs.ClebschGraph()
        Tgraph = Graph(Tau.n)
        for e1,e2 in Tg.edges:
            Tgraph.add_edge((e1-1,e2-1))
            
        autB = Bgraph.automorphism_group()
        card_autB = autB.cardinality()
        count_T_in_B = Bgraph.subgraph_search_count(Tgraph, induced=True)
        if count_T_in_B == card_autB: # there's exactly 1 strong hom (up to automorph grp of tau)
            Tcopy_in_B = Bgraph.subgraph_search(Tgraph, induced=True)
            strong_hom = Tcopy_in_B.vertices()
            claim5 = True
            print("\033[32m[OK]   \033[mThere is exactly 1 strong homomorphism from tau into B.")
        else:
            print("\033[31m[FAIL] \033[mThe number of strong homomorphisms from tau to B is wrong.")

    else:
        Bg = Flag(B)
        Bgraph = Bg.minimal_isomorph()
        Tgraph = Tg.minimal_isomorph()
        # set-up
        otuples = Tuples(range(1, Bgraph.n+1), Tgraph.n)
        possible_edges_t = Combinations(range(1,Tgraph.n+1), 2) # edge size = 2, assume
        coTgraph = Flag(Tgraph.__repr__()) # easiest way to make a copy
        coTgraph.edges = tuple([tuple(x) for x in possible_edges_t if tuple(x) not in Tgraph.edges])
        possible_edges_B = Combinations(range(1, Bgraph.n+1), 2)
        coBgraph = Flag(Bgraph.__repr__())
        coBgraph.edges = tuple([tuple(x) for x in possible_edges_B if tuple(x) not in Bgraph.edges])

        if Bgraph.n > 0:
            # make use of Sage functions
            sageB = Graph(Bgraph.n)
            for e in Bgraph.edges:
                sageB.add_edge(e)
            Baut_group_order = sageB.automorphism_group().order()
            strong_hom_count = 0
            for tpl in otuples:
                # for each map into B, check if it induces T
                edge_missing = False
                for edge in Tgraph.edges:
                    if edge_missing == True:
                        break
                    imedge1 = (tpl[edge[0]-1],tpl[edge[1]-1])
                    imedge2 = (tpl[edge[1]-1],tpl[edge[0]-1])
                    if imedge1 in Bgraph.edges or imedge2 in Bgraph.edges:
                        continue
                    else:
                        edge_missing = True
                        break
                if edge_missing==True:
                    continue # go to next perm
                coedge_missing = False
                for coedge in coTgraph.edges:
                    if coedge_missing == True:
                        break
                    imcoedge1 = (tpl[coedge[0]-1],tpl[coedge[1]-1])
                    imcoedge2 = (tpl[coedge[1]-1],tpl[coedge[0]-1])
                    if imcoedge1 in coBgraph.edges or imcoedge2 in coBgraph.edges:
                        continue
                    else:
                        coedge_missing = True
                        break

                if coedge_missing or edge_missing:
                    continue # this wasn't a strong hom embedding of tau into B
                else:
                    strong_hom = tpl
                    strong_hom_count += 1

            if strong_hom_count == Baut_group_order: # there's exactly 1 strong hom (up to automorph grp of tau)
                claim5 = True
                strong_hom = [x-1 for x in strong_hom]
                print("\033[32m[OK]   \033[mThere is exactly 1 strong homomorphism from tau into B.")
            else:
                print("\033[31m[FAIL] \033[mThe number of strong homomorphisms from tau to B is wrong.")



    # Different vertices of B attach differently to an embedding of tau in it
    # -----------------------------------------------------------------------
    claim6 = False
    
    found_identical_neighbourhoods = False
    witness_neighbourhood = None
    
    # take strong_hom from previous search
    if strong_hom == None:
        print("\033[31m[FAIL] \033[mThere is no strong homomorphism from tau into B. There should be exactly 1.")

    else:
        strong_hom = set(strong_hom)

        # prepare B and Tau
        Tg = Flag(tau)
        Tg = Tg.minimal_isomorph()
        Tau = Flag("6:1223344551")
        Tau = Tau.minimal_isomorph()

        # treat CLEBSCH case separately, Sage has tools
        if B == "g:12131415162728292a373b3c3d484b4e4f595c5e5g6a6d6f6g7e7f7g8c8d8g9b9d9fabacaebgcfde" and Tg == Tau:
            Bgraph = graphs.ClebschGraph()
        else:
            Bg = Flag(B)
            Bg = Bg.minimal_isomorph()
            Bg = Bg.__repr__()
            Bgraph = Graph(make_number(Bg[0]))
            for i in range(2,len(B),2):
                Bgraph.add_edge((make_number(Bg[i])-1, make_number(Bg[i+1])-1))

        # compute neighbourhoods for vertices in B intersected with embedding of Tau (by strong_hom)
        attachments_in_T = list()
        for v in Bgraph.vertices():
            vneigh =  Bgraph.vertex_boundary({v})
            attachments_in_T.append(strong_hom.intersection(vneigh))

        witness = list()
        for i in range(len(attachments_in_T)-1):
            for j in range(i+1,len(attachments_in_T)):
                if attachments_in_T[i] == attachments_in_T[j]:
                    witness.append((i,j))
                    witness.append(list(attachments_in_T[i]))

        if not witness:
            claim6 = True
            print("\033[32m[OK]   \033[mDifferent vertices of B attach differently to an embedding of tau in B.")
        else:
            print("\033[31m[FAIL] \033[mThere are two vertices", str(witness[0][0])+",",witness[0][1],"that attach to Tau in the same way:", str(witness[1])+".")


    # Check if forbidding tau improves the bound
    # ------------------------------------------
    claim7 = False
    
    print("Verifying that forbidding tau improves the bound...")

    if tau == "1:":
        bound_tau = 0
    else:
        command = "sage -python inspect_certificate.py "+certificate_filename_tau+" --log"
        try:
            f = subprocess.call(command, shell=True)
        except ValueError:
            print("Ooops! Things went wrong! Bound probably not written into 'bound.txt'.")

        bfile = open("bound.txt", 'r')
        bound_tau = Rational(bfile.readlines()[0])

    if minimize == False and bound_tau < bound:
        claim7 = True
        print("\033[32m[OK]   \033[m"+str(bound), "= lambda(Forb(\cal F)) > lambda(Forb(\cal F and tau)) =", str(bound_tau)+".")
    elif minimize == True and bound_tau > bound:
        claim7 = True
        print("\033[32m[OK]   \033[m"+str(bound), "= lambda(Forb(\cal F)) < lambda(Forb(\cal F and tau)) =", str(bound_tau)+".")
    else:
        print("\033[31m[FAIL] \033[m"+str(bound), "= lambda(Forb(\cal F)) while lambda(Forb(\cal F and tau)) =", str(bound_tau)+".")

        
    # PERFECT STABILITY

    # only hope for perfect stability, if problem robustly stable
    if not (claim1 and claim2 and claim3 and claim4 and claim5 and claim6 and claim7):
        raise ValueError("Robust stability has not been verified. Please do that first.")
        sys.exit()

    print("Robust stability verified!")
    
    # Verifying CLAIM 3: rk(Q_tau) = dim(Q_tau)-1
    # ---------------------------------------------------
    claimQ = False

    Tgraph = Flag(tau)
    Tgraph = Tgraph.minimal_isomorph()
    itau = types.index(Tgraph)
    dminus = numpy.array(Qs[itau]).shape[0]-1
    Q_tau = [[0 for x in range(dminus+1)] for y in range(dminus+1)]
    k = 0
    for i in range(dminus+1):
        for j in range(i,dminus+1):
            Q_tau[i][j] = Qs[itau][i][j-i]
            Q_tau[j][i] = Q_tau[i][j]
        
    Q_tau = matrix(QQ, Q_tau)
    rk =  rank(Q_tau)

    if rk == dminus:
        claimQ = True
        print("\033[32m[OK]   \033[mMatrix Q_tau has rk(Q_tau) = dim(Q_tau)-1 =", str(rk)+".")
        print("Perfect stability verified! Done.")
    else:
        print("\033[31m[FAIL] \033[mMatrix Q_tau has rank", rk, "and dim", str(dminus+1)+".")



    if not claimQ:
        # Check if forbidding B improves the bound (only if claimQ fails)
        claimB = False

        try:
            certB_filename = args[5].strip()
        except IOError:
            print("You did not provide a certificate file to verify perfect stability. Giving up...")
            sys.exit()

        command = "sage -python inspect_certificate.py "+certB_filename+" --log"
        try:
            f = subprocess.call(command, shell=True)
        except ValueError:
            print("Ooops! Things went wrong! Bound probably not written into 'bound.txt'.")

        try:
            bfile = open("bound.txt", 'r')
            bound_tau = Rational(bfile.readlines()[0])
            if minimize == False and bound_tau < bound:
                claimB = True
                print("\033[32m[OK]   \033[m"+str(bound), "= lambda(Forb(\cal F)) > lambda(Forb(\cal F and B)) =", str(bound_tau)+".")
            elif minimize == True and bound_tau > bound:
                claimB = True
                print("\033[32m[OK]   \033[m"+str(bound), "= lambda(Forb(\cal F)) < lambda(Forb(\cal F and B)) =", str(bound_tau)+".")
            else:
                print("\033[31m[FAIL] \033[m"+str(bound), "= lambda(Forb(\cal F)) while lambda(Forb(\cal F and B)) =", str(bound_tau)+".")

        except ValueError:
            print("Couldn't open file bound.txt or read the bound.")
