claim4 = False

# sharp graphs
sharp_graphs = list()
f = open("sg73N8.txt", 'r')
for line in f:
    sharp_graphs.append(line.split(' ')[1])
    
# convert sharp graphs to Sage form
sharp_graphs_sage = list()
for sg in sharp_graphs:
    g_str = sg.__repr__()
    g = Graph(8)
    for i in range(2,len(g_str),2):
        g.add_edge((int(g_str[i])-1, int(g_str[i+1])-1))
    sharp_graphs_sage.append(g)

# IDEA: contract twins until no more twins; then check if subgraph of B
Bgraph = Flag(B)
Bgraph = Bgraph.minimal_isomorph()
Gclebsch = Flag("g:12131415162728292a373b3c3d484b4e4f595c5e5g6a6d6f6g7e7f7g8c8d8g9b9d9fabacaebgcfde") # Clebsch graph
Gclebsch = Gclebsch.minimal_isomorph()
if Bgraph == Gclebsch:
    C = graphs.ClebschGraph()
    "was here"
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
    print "\033[32m[OK]   \033[mAll sharp graphs admit strong hom into B."
else:
    print "\033[31m[FAIL] \033[mNOT all sharp graphs admit strong homomorphism into B."
    print "       e.g. no strong hom from", failed[0].edges(labels=None), "into B", "("+B+")."



# There is exactly one strong homomorphism from tau to B (up to automorph of tau)
# -------------------------------------------------------------------------------
claim5 = False
strong_hom = None # will store the unique strong hom if found
"""
Tg = Flag(tau)
Bg = Flag(B)
Gclebsch = Flag("g:12131415162728292a373b3c3d484b4e4f595c5e5g6a6d6f6g7e7f7g8c8d8g9b9d9fabacaebgcfde") # Clebsch graph
Tg = Tgraph.minimal_isomorph()
Tau = Flag("6:1223344551")
Tau = Tau.minimal_isomorph()
Bg = Bg.minimal_isomorph()
Gclebsch = Gclebsch.minimal_isomorph()
"""
Bgraph = graphs.ClebschGraph()
Tgraph = Graph(6)
Tgraph.add_edge((0,1))
Tgraph.add_edge((1,2))
Tgraph.add_edge((2,3))
Tgraph.add_edge((3,4))
Tgraph.add_edge((4.0))

autB = Bgraph.automorphism_group()
card_autB = autB.cardinality()
count_T_in_B = Bgraph.subgraph_search_count(Tgraph, induced=True)
if count_T_in_B == card_autB: # there's exactly 1 strong hom (up to automorph grp of tau)
    Tcopy_in_B = Bgraph.subgraph_search(Tgraph)
    strong_hom = Tcopy.vertices()
    claim5 = True
    print "\033[32m[OK]   \033[mThere is exactly 1 strong homomorphism from tau into B."
else:
    print "\033[31m[FAIL] \033[mThe number of strong homomorphisms from tau to B is wrong."

