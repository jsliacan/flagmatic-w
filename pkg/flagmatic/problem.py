# -*- coding: utf-8 -*-
"""

flagmatic 2

Copyright (c) 2012, E. R. Vaughan. All rights reserved.

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


Further development of Flagmatic is supported by ERC.
http://cordis.europa.eu/project/rcn/104324_en.html
"""

import gzip, json, os, sys
import numpy
import itertools
import pexpect
import sage.all

from sage.structure.sage_object import SageObject
from sage.rings.all import Integer, Rational, QQ, ZZ, RDF
from sage.functions.other import floor
from sage.matrix.all import matrix, identity_matrix, block_matrix, block_diagonal_matrix
from sage.modules.misc import gram_schmidt
from sage.misc.misc import SAGE_TMP
#from sage.combinat.all import Permutations, Combinations, Tuples
from sage.matrix.constructor import ones_matrix, vector
from copy import copy

from .hypergraph_flag import make_graph_block, print_graph_block
from .flag import *
from .three_graph_flag import *
from .graph_flag import *
from .oriented_graph_flag import *
from .multigraph_flag import *
from .construction import *
from .blowup_construction import *

# pexpect in Sage 4.8 has a bug, which prevents it using commands with full paths.
# So for now, CSDP has to be in a directory in $PATH.

cdsp_cmd = "csdp"
sdpa_cmd = "sdpa"
sdpa_dd_cmd = "sdpa_dd"
sdpa_qd_cmd = "sdpa_qd"
dsdp_cmd = "dsdp"

def process_products_mp(tg, flag, n, flag_cls, graphs):
    graph_block = make_graph_block(graphs, n)
    
    s = tg.n
    m = (n + s) / 2
    
    flags_block = make_graph_block(flag, m)
    rarray = flag_cls.flag_products(graph_block, tg, flags_block, None)
    
    return rarray

def generate_flags_mp(flag_cls, m, tg, forbidden_edge_numbers, forbidden_graphs, forbidden_induced_graphs):
    return flag_cls.generate_flags(m, tg, forbidden_edge_numbers=forbidden_edge_numbers, forbidden_graphs=forbidden_graphs, forbidden_induced_graphs=forbidden_induced_graphs)

def zero_eigenvectors_mp(construction, types, flags):
    return construction.zero_eigenvectors(types, flags)

def compute_densities_mp(g, dg):
    dv = 0
    for h, coeff in dg:
        if h.n == g.n:
            if g == h:
                dv += coeff
        else:
            dv += coeff * g.subgraph_density(h)
    return dv

def block_structure(M):
    """
    Given a matrix, this function returns a tuple. The first entry is the number of
    row subdivisions that the matrix has. The second entry is a list of the sizes of the
    row subdivisions, and the third entry is a list containing the rows at which each
    subdivision begins. (Note that column subdivisions are ignored.)
    """
    row_div = M.subdivisions()[0]
    div_offsets = [0] + row_div
    div_sizes = row_div + [M.nrows()]
    num_blocks = len(div_sizes)

    for i in range(1, num_blocks):
        div_sizes[i] -= div_sizes[i - 1]

    return num_blocks, div_sizes, div_offsets


def safe_gram_schmidt(M):
    """
    Performs Gram Schmidt orthogonalization using Sage functions. The returned matrix
    will have the same subdivisions are the input matrix, and it will be sparse if and
    only if the input matrix is sparse.

    In addition, the new matrix.gram_schmidt method appears to have certain issues, at
    least under Sage 4.8, when dealing with number fields (in particular the ones used
    by the 'maxs3' and 'maxs4' problems.

    So, this function uses the old gram_schmidt function, whenever the base ring of the
    matrix is not the rational field.
    """

    subdivision = M.subdivisions()
    BR = M.base_ring()
    sparse = M.is_sparse()

    if BR == QQ:
        M, mu = M.gram_schmidt()
        # .gram_schmidt doesn't appear to preserve sparsity, so recreate matrix from rows.
        M = matrix(QQ, M.rows(), sparse=sparse)

    else:  # .gram_schmidt is broken for some number fields in 4.8.
        rows, mu = gram_schmidt(M.rows())
        M = matrix(BR, rows, sparse=sparse)

    M.subdivide(subdivision)
    return M


def LDLdecomposition(M):  # TODO: does this handle matrices with zero eigenvalues?
    MS = M.parent()
    D = MS.matrix()
    if M.is_zero():
        D.set_immutable()
        return D, D
    L = copy(MS.identity_matrix())
    for i in xrange(M.nrows()):
        for j in xrange(i):
            L[i, j] = (Integer(1) / D[j, j]) * (M[i, j] - sum(L[i, k] * L[j, k] * D[k, k] for k in xrange(j)))
        D[i, i] = M[i, i] - sum(L[i, k] ** 2 * D[k, k] for k in xrange(i))
    L.set_immutable()
    D.set_immutable()
    return L, D


class Problem(SageObject):
    r"""
    This is the principal class of flagmatic. Objects of this class represent Turán-type
    problems.
    """

    def __init__(self, flag_cls, order=None, forbid_induced=None, forbid=None,
                 forbid_homomorphic_images=False, density=None, minimize=False,
                 type_orders=None, types=None, max_flags=None, compute_products=True,
                 mode="plain"):
        
        r"""
        Creates a new Problem object. Generally it is not necessary to call this method
        directly, as Problem objects are more easily created using the helper functions:

        sage: problem = GraphProblem()
        sage: problem = ThreeGraphProblem()
        sage: problem = OrientedGraphProblem()
        sage: problem = TwoMultigraphProblem()
        sage: problem = ThreeMultigraphProblem()

        If this method is called directory, then a class that inherits from Flag should be
        provided, for example:

        sage: problem = Problem(GraphFlag)

        INPUT:
        
        - order: order of admissible graphs
        - forbid_induced: list of graph-strings specifying forbidden induced subgraphs
        - forbid: list of graph-srings specifying forbidden subgraphs
        - forbid_homomorphic_images:
        - density: density quantum graph (linear combination of graphs)
        - minimize: set to True if minimization problem, else maximization
        - type_orders: list of orders for types to be used
        - types: list of types to be used
        - max_flags:
        - compute_products: set to True if need to compute flag products
        - mode: plain/optimization/feasibility (see Flagmatic documentation)
        """

        self._flagmatic_version = "2.0"

        if issubclass(flag_cls, Flag):
            self._flag_cls = flag_cls
        else:
            raise ValueError

        self._stable = False
        self._robustly_stable = False
        self._perfectly_stable = False
        
        self._n = 0
        self._mode = "plain"
        
        self._field = QQ
        self._approximate_field = RDF

        self._forbidden_edge_numbers = []
        self._forbidden_graphs = []
        self._forbidden_induced_graphs = []

        self._assumptions = []
        self._assumption_flags = []
        
        self.state("specify", "yes")
        self.set_objective(minimize=minimize)

        # following wouldn't need to be stored
        # but want them when problem is re-run for stability purposes
        self._forbid_homomorphic_images = forbid_homomorphic_images
        self._max_flags = max_flags
        self._compute_products = compute_products
        self._type_orders = type_orders
        self._types_from_input = types
        
        if density is None:
            self.set_density(flag_cls.default_density_graph())
        else:
            self.set_density(density)
        

        if not forbid_induced is None:
            self.forbid_induced(forbid_induced)

        if not mode is None:
            if mode == "optimization": self._mode = mode
            elif mode == "feasibility": self._mode = mode
            elif not (mode == "plain"):
                raise ValueError
                
            
            
        if not forbid is None:
            self.forbid(forbid)

        if forbid_homomorphic_images:
            self.forbid_homomorphic_images()

            
        if not order is None:
            self.generate_flags(order, type_orders=type_orders, types=types, max_flags=max_flags, compute_products=compute_products)
            

    def state(self, state_name=None, action=None):
        r"""
        Keeps track of which things have been done. To get a list of all the states, enter

        sage: problem.state()

        If the argument state_name is supplied, then the value of the state called
        state_name is returned. This will be either "yes", "no" or "stale". "yes" means
        that the thing has been done, "no" means that it has not. "stale" means that it
        has been done, but subsequent actions mean that it needs to be re-done.

        If an action is supplied, it will change the value of the state state_name. The
        action must be one of "yes", "no" or "stale". It is not recommended that the
        value of any states be changed by the user.

        Setting ``force`` to True allows states to be set irrespective of other
        states.

        """

        # We use tuples here because there is an order to the states.
        state_list = [
            ("specify", {
                "requires": [],
                "depends": []
            }),
            ("set_objective", {
                "requires": [],
                "depends": []
            }),
            ("compute_flags", {
                "requires": [],
                "depends": ["specify"]
            }),
            ("set_construction", {
                "requires": ["compute_flags"],
                "depends": ["set_objective"]
            }),
            ("add_zero_eigenvectors", {
                "requires": ["set_construction"],
                "depends": []
            }),
            ("compute_block_bases", {
                "requires": ["compute_flags"],
                "depends": []
            }),
            ("compute_flag_bases", {
                "requires": ["set_construction"], # makes block bases if required
                "depends": ["add_zero_eigenvectors"]
            }),
            ("compute_products", {
                "requires": ["compute_flags"],
                "depends": []
            }),
            ("set_active_types", {
                "requires": ["compute_flags"],
                "depends": []
            }),
            ("set_block_matrix_structure", {
                "requires": ["compute_flags"],
                "depends": ["set_active_types"]
            }),
            ("write_sdp_input_file", {
                "requires": ["set_block_matrix_structure"],
                "depends": ["set_objective", "set_active_types"]
            }),
            ("write_sdp_initial_point_file", {
                "requires": ["set_block_matrix_structure"],
                "depends": ["set_objective", "set_active_types"]
            }),
            ("run_sdp_solver", {
                "requires": ["write_sdp_input_file"],
                "depends": ["write_sdp_input_file", "write_sdp_initial_point_file"]
            }),
            ("read_solution", {
                "requires": ["run_sdp_solver"],
                "depends": []
            }),
            ("check_solution", {
                "requires": ["read_solution"],
                "depends": []
            }),
            ("transform_solution", {
                "requires": ["compute_flag_bases", "check_solution"],
                "depends": []
            }),
            ("add_sharp_graphs", {
                "requires": ["set_construction"],
                "depends": []
            }),
            ("make_exact", {
                "requires": ["check_solution"],
                "depends": ["transform_solution", "add_sharp_graphs"]
            }),
            ("meet_target_bound", {
                "requires": ["transform_solution", "make_exact"],
                "depends": []
            }),
            ("check_exact", {
                "requires": ["make_exact"],
                "depends": ["meet_target_bound"]
            }),
            ("diagonalize", {
                "requires": ["meet_target_bound", "check_exact"],
                "depends": []
            }),
            ("write_certificate", {
                "requires": ["check_exact"],
                "depends": ["diagonalize"]
            })
        ]

        state_names = [s[0] for s in state_list]
        states = dict(state_list)

        if not hasattr(self, "_states"):
            self._states = dict((sn, "no") for sn in state_names)

        if state_name is None:
            width = max(len(sn) for sn in state_names)
            for sn in state_names:
                sys.stdout.write(("%%-%ds : %%s\n" % width) % (sn, self._states[sn]))
            return

        if not state_name in state_names:
            raise ValueError("unknown state.")

        if action == "yes" or action == "force_yes":
            if action == "yes" and not all(self._states[sn] == "yes" for sn in states[state_name]["requires"]):
                raise NotImplementedError("not ready for this yet!")
            self._states[state_name] = "yes"
            for sn in states:
                if state_name in states[sn]["depends"] and self._states[sn] == "yes":
                    self.state(sn, "stale")

        elif action == "stale":
            self._states[state_name] = "stale"
            for sn in states:
                if state_name in states[sn]["depends"]:
                    self.state(sn, "stale")

        elif action == "ensure_yes":
            if self._states[state_name] != "yes":
                raise NotImplementedError("not ready for this yet!")

        elif action == "ensure_yes_or_stale":
            if self._states[state_name] not in ["yes", "stale"]:
                raise NotImplementedError("not ready for this yet!")

        elif action is None:
            pass

        else:
            raise ValueError("unknown action.")

        return self._states[state_name]


    @property
    def flag_cls(self):
        return self._flag_cls


    @property
    def order(self):
        return self._n

    # Deprecated - use order instead.
    @property
    def n(self):
        return self._n

    # TODO: sanity checking of type orders

    def generate_flags(self, order, type_orders=None, types=None, max_flags=None, compute_products=True):
        r"""
        Generates the types and flags that will be used in the problem.

        INPUT:

         - ``order`` -- an integer. The order for the unlabelled flags (the admissible
           graphs).

         - ``type_orders`` -- (default: None) None, or a list of integers. If a list of
           integers is given, then types of these orders will be generated. If None is
           given then types of all orders less than ``order``, and congruent modulo 2 to
           ``order`` will be generated.

         - ``types`` - (default: None) A list of flags to use as types.

         - ``max_flags`` -- (default: None) None, or an integer. If an integer is given,
           then it provides an upper bound on the number of flags each type can have. If a
           type has more than this many flags, it will be removed.

         - ``compute_products`` -- (default: True) Boolean. If True then the flag products
           will be computed. For some large problems this may take a long time. If False,
           then the flag products must be computed later using the ``compute_products``
           method.
        """

        n = order
        if type_orders is None:
            orders = [(s, floor((n + s) / 2)) for s in range(n % 2, n - 1, 2)]
        else:
            orders = []
            for to in type_orders:
                if type(to) is tuple:
                    orders.append(to)
                else:
                    orders.append((to, floor((n + to) / 2)))

        self.state("compute_flags", "yes")
        self._n = n

        sys.stdout.write("Generating graphs...\n")
        self._graphs = self._flag_cls.generate_graphs(n, forbidden_edge_numbers=self._forbidden_edge_numbers,
                                                      forbidden_graphs=self._forbidden_graphs, forbidden_induced_graphs=self._forbidden_induced_graphs,
                                                      use_mp=True, show_progress=True)
        sys.stdout.write("Generated %d graphs.\n" % len(self._graphs))

        for g in self._graphs:    # Make all the graphs immutable
            g.set_immutable()

        sys.stdout.write("Computing density...\n")
        self._compute_densities()

        sys.stdout.write("Generating types and flags...\n")
        self._types = []
        self._flags = []

        allowed_types = []
        if types:
            for h in types:
                if isinstance(h, str):
                    h = self._flag_cls(h)
                if not isinstance(h, self._flag_cls):
                    raise ValueError
                allowed_types.append(h)

        for s, m in orders:

            these_types = self._flag_cls.generate_graphs(s, forbidden_edge_numbers=self._forbidden_edge_numbers,
                                                         forbidden_graphs=self._forbidden_graphs,
                                                         forbidden_induced_graphs=self._forbidden_induced_graphs)

            if types:
                these_types = [h for h in these_types if h in allowed_types]

            sys.stdout.write("Generated %d types of order %d, " % (len(these_types), s))

            import multiprocessing as mp
            arguments = [(self._flag_cls, m, tg, self._forbidden_edge_numbers, self._forbidden_graphs, self._forbidden_induced_graphs) for tg in these_types]
            p = mp.Pool()
            these_flags = p.starmap(generate_flags_mp, arguments)
            p.close()
            # these_flags = []
            # for tg in these_types:
            #     these_flags.append(self._flag_cls.generate_flags(m, tg, forbidden_edge_numbers=self._forbidden_edge_numbers,
            #                                                      forbidden_graphs=self._forbidden_graphs,
            #                                                      forbidden_induced_graphs=self._forbidden_induced_graphs))
            sys.stdout.write("with %s flags of order %d.\n" % (sum([len(L) for L in these_flags]), m))

            self._types.extend(these_types)
            self._flags.extend(these_flags)

        num_types = len(self._types)

        if not max_flags is None:
            bad_indices = [i for i in range(num_types) if len(self._flags[i]) > max_flags]
            if len(bad_indices) > 0:
                good_indices = [i for i in range(num_types) if not i in bad_indices]
                self._types = [self._types[i] for i in good_indices]
                self._flags = [self._flags[i] for i in good_indices]
                sys.stdout.write("Removed types %s as they have too many flags.\n" % bad_indices)

        num_types = len(self._types)  # may have changed!

        self._active_types = range(num_types)

        for ti in range(num_types):              # Make everything immutable!
            self._types[ti].set_immutable()
            for g in self._flags[ti]:
                g.set_immutable()

        if compute_products:
            self.compute_products()


    @property
    def graphs(self):
        r"""
        Read-only. A (copy) of the list of admissible graphs. Modifying the list will have
        no effect on the Problem.

        """
        return copy(self._graphs)

    @property
    def types(self):
        """
        Read-only. A (copy) of the list of types. Modifying the list will have
        no effect on the Problem.

        """
        return copy(self._types)

    @property
    def flags(self):
        """
        Read-only. A (copy) of the flags, as a list of lists. The flags are listed by
        type, in the same order as the types appear in the types list. Modifying the lists
        will have no effect on the Problem.

        """
        return copy(self._flags)

    @property
    def density_graphs(self):
        """
        Read-only. A (copy) of the list of density graphs. Modifying the lists
        will have no effect on the Problem.

        """
        return copy(self._density_graphs)

    def set_objective(self, minimize=False):
        r"""
        Sets the Problem to be a "maximization" or a "minimization" problem.

        INPUT:

         - ``minimize`` - boolean (default: False) sets whether the problem is a
           minimization problem (or a maximization problem).

        In a maximization problem, the objective is to find a good (i.e. low) upper bound
        on a density. In a minimization problem the objective is to find a good (i.e.
        high) lower bound on a density.
        """
        if not type(minimize) is bool:
            raise ValueError

        self.state("set_objective", "yes")

        self._minimize = minimize


    def clear_densities(self):
        
        self._density_graphs = []
        self._active_densities = []
        self._density_coeff_blocks = []
        
        self._compute_densities()


        
    def _compute_densities(self):

        import multiprocessing as mp
        p = mp.Pool()
        self._densities = []
        for dg in self._density_graphs:
            arguments = [(g, dg) for g in self._graphs]
            density_values = p.starmap(compute_densities_mp, arguments)
            # density_values = []
            # for g in self._graphs:
            #     dv = 0
            #     for h, coeff in dg:
            #         if h.n == g.n:
            #             # comparison will be fast, as both g and h should have
            #             # _certified_minimal_isomorph set to True
            #             if g == h:
            #                 dv += coeff
            #         else:
            #             dv += coeff * g.subgraph_density(h)
            #     density_values.append(dv)
            self._densities.append(density_values)
        p.close()
            
    def set_density(self, *args):

        self.state("set_objective", "yes")

        flattened_args = []
        for arg in args:
            if isinstance(arg, list):
                flattened_args.extend(arg)
            else:
                flattened_args.append(arg)

        density_graphs = []

        for arg in flattened_args:

            if isinstance(arg, str) and "." in arg:
                arg = tuple(map(int, arg.split(".")))

            if isinstance(arg, tuple):

                if len(arg) != 2:
                    raise ValueError

                if arg[0] in ZZ:

                    k, ne = arg
                    if k < self._flag_cls().r:
                        raise ValueError
                    max_e = self._flag_cls.max_number_edges(k)
                    if not ne in range(max_e + 1):
                        raise ValueError

                    # Don't forbid anything - if we do, we'll have to keep list
                    # updated whenever forbidden things change. Also it seems the most
                    # appropriate thing to do.
                    graphs = self._flag_cls.generate_graphs(k)
                    for g in graphs:
                        if g.ne == ne:
                            density_graphs.append((g, Integer(1)))
                    continue

                else:
                    h, coeff = arg

            else:
                h, coeff = arg, Integer(1)

            if isinstance(h, str):
                h = self._flag_cls(h)

            if not isinstance(h, self._flag_cls):
                raise ValueError

            h = copy(h)
            h.make_minimal_isomorph()
            density_graphs.append((h, coeff))

        if len(density_graphs) == 0:
            raise ValueError

        # Note that this function only sets one of the densities.
        self._density_graphs = [density_graphs]
        self._active_densities = [0]
        self._density_coeff_blocks = [[0]]

        if self.state("compute_flags") == "yes":
            self._compute_densities()

    def _forbid(self, h, induced):

        self.state("specify", "yes")

        if isinstance(h, str) and "." in h:
            h = tuple(map(int, h.split(".")))

        if isinstance(h, tuple):
            k, ne = h
            if k < self._flag_cls().r:
                raise ValueError
            max_e = self._flag_cls.max_number_edges(k)
            if not ne in range(max_e + 1):
                raise ValueError
            if induced:
                self._forbidden_edge_numbers.append((k, ne))
                sys.stdout.write("Forbidding %d-sets from spanning exactly %d edges.\n" % (k, ne))
            else:
                for i in range(ne, max_e + 1):
                    self._forbidden_edge_numbers.append((k, i))
                sys.stdout.write("Forbidding %d-sets from spanning at least %d edges.\n" % (k, ne))
            return

        if isinstance(h, str):
            h = self._flag_cls(h)

        if not isinstance(h, self._flag_cls):
            raise ValueError

        if induced:
            self._forbidden_induced_graphs.append(copy(h))
            self._forbidden_induced_graphs.sort(key=lambda g: (g.n, g.ne))
            sys.stdout.write("Forbidding %s as an induced subgraph.\n" % h)
        else:
            self._forbidden_graphs.append(copy(h))
            self._forbidden_graphs.sort(key=lambda g: (g.n, g.ne))
            sys.stdout.write("Forbidding %s as a subgraph.\n" % h)

    def forbid(self, *args):
        r"""
        Sets the problem to be in the theory of graphs that do not contain the given subgraphs.

        INPUT:

         - arguments must be Flags, strings, or lists of these things. Flags must be of the
           appropriate class for the problem. Strings will be interpreted as representing flags
           of the appropriate class.

        """
        for h in args:
            if isinstance(h, list):
                for x in h:
                    self._forbid(x, False)
            else:
                self._forbid(h, False)


    def forbid_induced(self, *args):
        r"""
        Sets the problem to be in the theory of graphs that do not have contain the given
        graphs as induced subgraphs.

        INPUT:

         - arguments must be Flags, strings, or lists of these things. Flags must be of the
           appropriate class for the problem. Strings will be interpreted as representing flags
           of the appropriate class.

        """
        for h in args:
            if isinstance(h, list):
                for x in h:
                    self._forbid(x, True)
            else:
                self._forbid(h, True)


    def forbid_homomorphic_images(self):
        r"""
        Restricts the problem to be in the theory of graphs that do not contain homomorphic images
        of graphs already specified using ``forbid``. For certain problems this will make the
        computation simpler, without affecting the result.
        """
        self._forbid_homomorphic_images = True
        
        L = sum([g.homomorphic_images() for g in self._forbidden_graphs], [])
        LM = self._flag_cls.minimal_by_inclusion(L)
        if len(LM) == 0:
            return
        #sys.stdout.write("Forbidding")
        for g in LM:
            #sys.stdout.write(" %s" % repr(g))
            self._forbid(g, False)
        #sys.stdout.write("\n")

        # TODO: warn if already solved
    

    def add_assumption(self, typegraph, lincomb, const=0, equality=False):
    
        """
        Convert assumption from the general form:
        [linear combination of flags on one type] >= c   OR
        [linear combination of flags on one type] == c 

        into an assumptions of the form
        [linear combination of flags on one type] >= 0

        INPUT:
        
        - typegraph: # it is the common type or the entire assumption,
                     # e.g. "3:121323" for labelled triangle
        
        - lincomb: # this is the linear combination of terms (flag,
                   # coef) as a list, i.e. LHS of the assumption
        
        - const: # RHS of the assumption (should be some rational in
                 # form a/b or a)
        
        - equality: # whether the assumption is equality True or
                    # inequality False; default is False

        EXAMPLE:
         
        problem = GraphProblem(4, mode="optimization")
        problem.add_assumption("0:", [("2:12(0)", 1)], 1/2, equality=True)
        
        """


        if self._mode == "plain":
            sys.stdout.write("\nCannot add assumptions in 'plain' mode.\n")
            sys.stdout.write("Change mode?\n")
            sys.stdout.write("\t0\t stay in 'plain' mode\n")
            sys.stdout.write("\t1\t change to 'optimization' mode\n")
            sys.stdout.write("\t2\t change to 'feasibility' mode\n")
            choice = raw_input("Enter your choice now and press return: ")

            if choice == '1':
                self._mode = "optimization"
            elif choice == '2':
                self._mode = "feasibility"
            elif choice == '0':
                sys.stdout.write("Staying in 'plain mode'. Cannot add assumptions!\n")
                return
            else:
                raise ValueError

        # PARSING ASSUMPTIONS FROM INPUT
        
        if self._flag_cls().r == 2:

            if self._flag_cls().oriented:
                
                try:
                    tg = OrientedGraphFlag(typegraph)
                    tg.make_minimal_isomorph()

                    cst = Rational(const)
                    eq = equality
                    indep = False

                    lcomb = [[OrientedGraphFlag(g), Rational(c)] for g,c in lincomb]
                    for term in lcomb: term[0].make_minimal_isomorph()  # convert flag to the one Flagmatic knows

                    # if RHS nonzero, add a type to the LHS with -const coefficient (works with '0:' type as well)
                    if cst:
                        n = max([term[0].n for term in lcomb])
                        if (tg.n == 0) and (n == self._n): # then convert constant into a family H
                            for H in self._graphs:
                                is_in_lcomb = False
                                for term in lcomb:
                                    if H == term[0]:
                                        is_in_lcomb = True
                                        term[1] += -cst
                                if not is_in_lcomb:
                                    lcomb.append([H,-cst])
                        else:
                            fg = OrientedGraphFlag(tg._repr_() + "("+str(tg.n)+")")
                            lcomb.append([fg, -cst])

                    cst = 0 # make RHS zero. (not necessary though, not used again)

                    to_throw_out = [0 for term in lcomb]
                    for i in range(len(lcomb)):
                        if lcomb[i][1] == 0:
                            to_throw_out[i] = 1

                    counter = 0
                    for i in range(len(to_throw_out)):
                        if to_throw_out[i] == 1:
                            lcomb.pop(i-counter) # remove terms with coeff=0
                            counter += 1

                    lcomb = [tuple(x) for x in lcomb] # make [graph, coeff] into (graph, coeff)

                except ValueError:
                    print("You are trying to feed 'add_assumption()' unhealthy things!")

            else:

                try:
                    tg = GraphFlag(typegraph)
                    tg.make_minimal_isomorph()

                    cst = Rational(const)
                    eq = equality
                    indep = False

                    lcomb = [[GraphFlag(g), Rational(c)] for g,c in lincomb]
                    for term in lcomb: term[0].make_minimal_isomorph()  # convert flag to the one Flagmatic knows

                    # if RHS nonzero, add a type to the LHS with -const coefficient (works with '0:' type as well)
                    if cst:
                        n = max([term[0].n for term in lcomb])
                        if (tg.n == 0) and (n == self._n): # then convert constant into a family H
                            for H in self._graphs:
                                is_in_lcomb = False
                                for term in lcomb:
                                    if H == term[0]:
                                        is_in_lcomb = True
                                        term[1] += -cst
                                if not is_in_lcomb:
                                    lcomb.append([H,-cst])
                        else:
                            fg = GraphFlag(tg._repr_() + "("+str(tg.n)+")")
                            lcomb.append([fg, -cst])

                    cst = 0 # make RHS zero. (not necessary though, not used again)

                    to_throw_out = [0 for term in lcomb]
                    for i in range(len(lcomb)):
                        if lcomb[i][1] == 0:
                            to_throw_out[i] = 1

                    counter = 0
                    for i in range(len(to_throw_out)):
                        if to_throw_out[i] == 1:
                            lcomb.pop(i-counter) # remove terms with coeff=0
                            counter += 1

                    lcomb = [tuple(x) for x in lcomb] # make [graph, coeff] into (graph, coeff)

                except ValueError:
                    print("You are trying to feed 'add_assumption()' unhealthy things!")
                    
                
        elif self._flag_cls().r == 3:

            try:
                tg = ThreeGraphFlag(typegraph)
                tg.make_minimal_isomorph()

                cst = Rational(const)
                eq = equality
                indep = False
                
                lcomb = [[ThreeGraphFlag(g), Rational(c)] for g,c in lincomb]
                for term in lcomb: term[0].make_minimal_isomorph()  # convert flag to the one Flagmatic knows

                
                # if RHS nonzero, add a type to the LHS with -const coefficient (works with '0:' type as well)
                if cst:
                    n = max([term[0].n for term in lcomb])
                    if (tg.n == 0) and (n == self._n): # then convert constant into a family H
                        for H in self._graphs:
                            for term in lcomb:
                                if H == term[0]:
                                    term[1] -= cst
                                else:
                                    lcomb.append([H,-cst])
                    else:
                        fg = ThreeGraphFlag(tg._repr_() + "("+str(tg.n)+")")
                        lcomb.append([fg, -cst])
                        # make RHS zero. (not necessary though, not used again)

                cst = 0

                to_throw_out = [0 for term in lcomb]
                for i in range(len(lcomb)):
                    if lcomb[i][1] == 0:
                        to_throw_out[i] = 1

                for i in to_throw_out:
                    if i == 1:
                        lcomb.pop(i) # remove terms with coeff=0

                lcomb = [tuple(x) for x in lcomb] # make [graph, coeff] into (graph, coeff)
                
            except ValueError:
                print("You are trying to feed 'add_assumption()' unhealthy things!")
            

        # translate the assumption to the simple ones and add them one by one

        if equality: # assumption is equality

            indep = True # in this case assumption does not to go the objective function
            minus_lcomb = [(g,-c) for g,c in lcomb]

            self._add_assumption(tg, lcomb, independent=indep)
            self._add_assumption(tg, minus_lcomb, independent=indep)


        else: # assumption is already inequality

            if self._mode == "optimization":
                self._add_assumption(tg, lcomb, independent=True)
            elif self._mode == "feasibility":
                self._add_assumption(tg, lcomb, independent=indep)

        
    def _add_assumption(self, tg, terms, independent=False):


        if self._mode == "feasibility" and (not self._assumptions):
            sys.stdout.write( "In feasibility mode now: not using density graphs.\n")
            self.clear_densities()
        elif self._mode == "feasibility":
            pass
        elif self._mode == "optimization":
            pass
        else:
            raise ValueError("Something is wrong!\n")
            
        self.state("set_objective", "yes")

        # treat 0 type separately
        # for the optimization mode
        """
        if tg.n == 0:
            num_densities = 1
        """
        
        m = self.n - max([t[0].n for t in terms]) + tg.n

        assumption_flags = self._flag_cls.generate_flags(m, tg, forbidden_edge_numbers=self._forbidden_edge_numbers, forbidden_graphs=self._forbidden_graphs, forbidden_induced_graphs=self._forbidden_induced_graphs)

        num_densities = len(assumption_flags)
        sys.stdout.write("Added %d quantum graphs.\n" % num_densities)
        
        num_graphs = len(self._graphs)
        quantum_graphs = [[Integer(0) for i in range(num_graphs)] for j in range(num_densities)]
        
        assumption_flags_block = make_graph_block(assumption_flags, m)
        graph_block = make_graph_block(self._graphs, self.n)
        
        for i in range(len(terms)):
            fg = terms[i][0]
            flags_block = make_graph_block([fg], fg.n)
            rarray = self._flag_cls.flag_products(graph_block, tg, flags_block, assumption_flags_block)
            
            for row in rarray:
                gi = row[0]
                j = row[1]  # always 0
                k = row[2]
                value = Integer(row[3]) / Integer(row[4])
                quantum_graphs[k][gi] += value * terms[i][1]

        self._assumptions.append((tg, terms))
        self._assumption_flags.append(assumption_flags)

        num_previous_densities = len(self._density_graphs)
        
        for qg in quantum_graphs:
            dg = []
            for gi in range(num_graphs):
                if qg[gi] != 0:
                    dg.append((self._graphs[gi], qg[gi]))
            self._density_graphs.append(dg)

        new_density_indices = range(num_previous_densities, num_previous_densities + len(quantum_graphs))
        self._active_densities.extend(new_density_indices)
        
        if not independent: # never happens in optimization mode
            # make assumptions look like one big assumption (some coeffs must be > 0)
            if not self._density_coeff_blocks:
                self._density_coeff_blocks.append(new_density_indices)
            else:
                self._density_coeff_blocks[0].extend(new_density_indices)

        sys.stdout.write("Re-computing densities...\n")
        sys.stdout.flush()
        self._compute_densities()


    def set_inactive_types(self, *args):
        r"""
        Specifies that the Q matrices for certain types should be zero matrices.

        INPUT:

        - arguments should be integers, specifying the indices of types in ``types``
          that should be marked as being "inactive".
        """
        self.state("set_active_types", "yes")

        num_types = len(self._types)

        for arg in args:
            ti = int(arg)
            if not ti in range(num_types):
                raise ValueError
            if ti in self._active_types:
                self._active_types.remove(ti)
            else:
                sys.stdout.write("Warning: type %d is already inactive.\n" % ti)


    def set_inactive_densities(self, *args):
        r"""
        Specifies that the coefficients of certain densities should be zero.
        
        INPUT:
        
        - arguments should be integers, specifying the indices of densities that should be
          marked as being "inactive".
        """
        for arg in args:
            di = int(arg)
            if not di in range(len(self._density_graphs)):
                raise ValueError
            if di in self._active_densities:
                self._active_densities.remove(di)
            else:
                sys.stdout.write("Warning: density %d is already inactive.\n" % di)


    def set_approximate_field(self, field):
        r"""
        Specifies the approximate field used when reading in and transforming the solution
        matrices.

        INPUT:

        - ``field`` - a field that uses finite precision, or in other words, a field that uses
          floating point arithmetic.

        EXAMPLES:

        sage: problem.set_approximate_field(RealField(113))

        This specifies that quadruple precision should be used. The default approximate field is
        the real dense field, RDF(), which uses 53 bits of precision.
        """
        if not field.is_field():
            raise ValueError("not a field.")

        if field.is_exact():
            raise ValueError("field must be floating point.")

        self._approximate_field = field


    # TODO: sanity check target_bound ?

    def set_extremal_construction(self, construction=None, field=None, target_bound=None):
        r"""
        Sets the extremal construction. This will be used to determine sharp graphs and forced
        zero eigenvectors.

        INPUT:

         - ``construction`` - a Construction object (default: None). None can be specified, in
           which case sharp graphs and forced zero eigenvectors can be added manually by using
           ``add_sharp_graphs`` and ``add_zero_eigenvectors``.

         - ``field`` - a field object (default: None). If ``construction`` is None, then this
           argument can be used to specify the exact field to use for computing the bound. This
           argument must be None if ``construction`` is not None, as the field will be taken from
           the Construction object.

         - ``target_bound`` - a number (default: None). If ``construction`` is None, then this
           argument can be used to specify the bound that should be aimed for. This argument must
           be None if ``construction`` is not None, as the target bound will be taken from
           the Construction object.
        """
        num_types = len(self._types)

        if construction is None:

            if not field.is_field():
                raise ValueError("not a valid field.")

            if not field.is_exact():
                raise ValueError("field must be an exact field (not floating point).")

            self.state("set_construction", "yes")

            self._construction = None
            self._field = field
            sys.stdout.write("Field is \"%s\" with embedding x=%s.\n" %
                            #(str(self._field), self._field.gen_embedding().n(digits=10)))
                             (str(self._field), self._field.gen().n(digits=10)))
            self._target_bound = target_bound
            sys.stdout.write("Set target bound to be %s (%s).\n" %
                (self._target_bound, self._target_bound.n(digits=10)))

            self._zero_eigenvectors = []
            for ti in range(num_types):
                M = matrix(self._field, 0, len(self._flags[ti]), sparse=True)
                M.set_immutable()
                self._zero_eigenvectors.append(M)

            self._sharp_graphs = []
            self._sharp_graph_densities = []

            return

        if not isinstance(construction, Construction):
            raise ValueError("not a valid construction.")

        if not field is None:
            raise ValueError("field should be None if construction is given.")

        if not target_bound is None:
            raise ValueError("target_bound should be None if construction is given.")

        self.state("set_construction", "yes")

        num_graphs = len(self._graphs)
        num_densities = len(self._densities)

        self._construction = construction
        self._field = construction.field

        sys.stdout.write("Determining which graphs appear in construction...\n")
        
        # TODO: mp subgraph densities
        sharp_graphs = construction.subgraph_densities(self._n)
        target_densities = [0 for j in range(num_densities)]
        self._sharp_graphs = []
        self._sharp_graph_densities = []

        for pair in sharp_graphs:
            g, den = pair
            # TODO: possibly mp this loop
            for gi in range(num_graphs):
                if g.is_labelled_isomorphic(self._graphs[gi]):
                    self._sharp_graphs.append(gi)
                    self._sharp_graph_densities.append(den)
                    for j in range(num_densities):
                        target_densities[j] += self._densities[j][gi] * den
                    break
            else:
                sys.stdout.write("Warning: non-admissible graph %s appears in construction!\n" % g)

        # set target_bound to equal the maximum - probably this will always be what is wanted...
        # print(target_densities)
        self._target_bound = max(target_densities)

        sys.stdout.write("Density of construction is %s.\n" % self._target_bound)

        from tqdm import tqdm
        import multiprocessing as mp
        arguments = [(construction, self._types[ti], self._flags[ti]) for ti in range(len(self._types))]
        p = mp.Pool()
        self._zero_eigenvectors = p.starmap(zero_eigenvectors_mp, tqdm(arguments))
        p.close()
        
        for ti in range(len(self._types)):
            sys.stdout.write("Found %d zero eigenvectors for type %d.\n" % (self._zero_eigenvectors[ti].nrows(), ti))
                
        # self._zero_eigenvectors = []
        # 
        # for ti in range(len(self._types)):
        # 
        #     self._zero_eigenvectors.append(construction.zero_eigenvectors(self._types[ti], self._flags[ti]))
        # 
        #     sys.stdout.write("Found %d zero eigenvectors for type %d.\n" % (
        #         self._zero_eigenvectors[ti].nrows(), ti))
        
        for ti in range(len(self._types)):
            self._zero_eigenvectors[ti].set_immutable()


    # TODO: reinstate the per-block option?

    def add_zero_eigenvectors(self, ti, eigenvectors, use_bases=False):
        r"""
        Adds a zero eigenvector. This method is necessary when not all the zero eigenvectors
        can be determined from the construction.

        INPUT:

         - ``ti`` - integer specifying which type the eigenvector is for.

         - ``eigenvectors`` - a vector, or matrix, of zero eigenvector(s) to add for type
           ``ti``. If adding more than one vector, the vectors can be given as the rows of a
           matrix. (Alternatively, they can be added one at a time with multiple calls to
           this method.)

         - ``use_bases`` - specifies that the vector is given in the basis of the Q' matrix, as
           opposed to the standard basis. The vector will be tranformed before being added to the
           zero eigenvectors.
        """
        self.state("add_zero_eigenvectors", "yes")

        if use_bases:
            self.state("compute_flag_bases", "ensure_yes_or_stale")
            NZ = (self._flag_bases[ti].T).solve_left(eigenvectors)
        else:
            NZ = eigenvectors

        self._zero_eigenvectors[ti] = self._zero_eigenvectors[ti].stack(NZ)
        self._zero_eigenvectors[ti].set_immutable()

    # TODO: is good idea to assume zero densities?

    def add_sharp_graphs(self, *args):
        r"""
        Adds a sharp graph. This method is necessary when not all the sharp graphs can be
        determined from the construction.

        INPUT:

         - arguments are integers specifying the indices of the graphs to set as sharp.
        """
        self.state("add_sharp_graphs", "yes")

        num_graphs = len(self._graphs)

        for arg in args:
            si = int(arg)
            if not si in range(num_graphs):
                raise ValueError
            if not si in self._sharp_graphs:
                self._sharp_graphs.append(si)
                self._sharp_graph_densities.append(Integer(0))
            else:
                sys.stdout.write("Warning: graph %d is already marked as sharp.\n" % si)

    def change_solution_bases(self, use_blocks=True):
        r"""
        Transforms the solution's Q matrices, so that they are (hopefully) positive definite. A
        construction should have been set previously, and this will be used to determine forced
        zero eigenvectors.

        This method is called from ``make_exact`` by default, and so is not normally explicitly
        invoked.

        INPUT:

         - ``use_blocks`` - Boolean (default: True). Specifies whether to apply an additional
           change of basis so that the matrices have a block structure with two blocks. This uses
           the invariant anti-invariant idea of Razborov.
        """

        if self.state("compute_flag_bases") != "yes":
            self.compute_flag_bases(use_blocks)

        self.state("transform_solution", "yes")

        num_types = len(self._types)

        sys.stdout.write("Transforming matrices")

        self._sdp_Qdash_matrices = []

        for ti in range(num_types):

            B = self._flag_bases[ti]
            if B.nrows() > 0:
                row_div = B.subdivisions()[0]
                M = B * self._sdp_Q_matrices[ti] * B.T
                M.subdivide(row_div, row_div)
                # zero out bits that should be zero. Note the copy() seems to be needed.
                M = block_diagonal_matrix([copy(M.subdivision(i,i)) for i in range(len(row_div) + 1)])
                M.set_immutable()
                self._sdp_Qdash_matrices.append(M)
            else:
                self._sdp_Qdash_matrices.append(matrix(self._approximate_field, 0, 0))
            sys.stdout.write(".")
            sys.stdout.flush()

        sys.stdout.write("\n")

    def compute_block_bases(self):
        r"""
        Computes a basis for the solution's Q matrices, that will give them a block structure with
        two blocks. This uses the invariant anti-invariant idea of Razborov. This method is not
        normally explicitly invoked. By default, ``change_problem_bases`` and
        ``change_solution_bases`` call it.
        """
        self.state("compute_block_bases", "yes")

        self._block_bases = []

        for ti in range(len(self._types)):

            B = self._flag_cls.flag_basis(self._types[ti], self._flags[ti])
            num_blocks, block_sizes, block_offsets = block_structure(B)
            sys.stdout.write("Type %d (%d flags) blocks: %s \n" % (ti, len(self._flags[ti]), block_sizes))
            self._block_bases.append(B)

    def compute_flag_bases(self, use_blocks=True, keep_rows=False, use_smaller=False):
        r"""
        Computes a basis for the solution's Q matrices, using the construction to determine forced
        zero eigenvectors. This method is used by ``change_problem_bases`` and
        ``change_solution_bases``, and would not usually be invoked directly.
        """
        self.state("compute_flag_bases", "yes")

        num_types = len(self._types)

        if use_blocks and self.state("compute_block_bases") != "yes":
            self.compute_block_bases()

        self._flag_bases = []

        sys.stdout.write("Creating bases")
        sys.stdout.flush()

        for ti in range(num_types):

            if use_blocks:
                num_blocks, block_sizes, block_offsets = block_structure(self._block_bases[ti])
            else:
                num_blocks, block_sizes, block_offsets = 1, [len(self._flags[ti])], [0]

            BS = []

            for bi in range(num_blocks):

                Z = self._zero_eigenvectors[ti]

                if use_blocks:
                    B = (self._block_bases[ti].subdivision(bi, 0) * Z.T).T
                else:
                    B = Z

                B = B.echelon_form()

                nzev = B.rank()
                B = B[:nzev, :]

                if nzev == 0:
                    B = identity_matrix(QQ, block_sizes[bi], sparse=True)

                elif nzev == block_sizes[bi]:
                    pass

                else:
                    B = B.stack(B.right_kernel().basis_matrix())

                if use_blocks:
                    B = B * self._block_bases[ti].subdivision(bi, 0)

                if not keep_rows:
                    B = B[nzev:, :]  # delete rows corresponding to zero eigenvectors

                if B.nrows() > 0:
                    BS.append(B)

            M = block_matrix([[B] for B in BS], subdivide=True)

            if M.nrows() == 0:
                M = matrix(self._field, 0, len(self._flags[ti]), sparse=True)
            else:
                M = safe_gram_schmidt(M)

            M.set_immutable()
            self._flag_bases.append(M)

            sys.stdout.write(".")
            sys.stdout.flush()

        sys.stdout.write("\n")

        self._inverse_flag_bases = []

        for ti in range(num_types):

            M = copy(self._flag_bases[ti])

            if use_smaller:
                MT = M.T
                MT.set_immutable()
                self._inverse_flag_bases.append(MT)
                for j in range(M.nrows()):
                    M[j, :] /= sum([x ** 2 for x in M.row(j)])
                M.set_immutable()
                self._flag_bases[ti] = M

            else:
                for j in range(M.nrows()):
                    M[j, :] /= sum([x ** 2 for x in M.row(j)])
                MT = M.T
                MT.set_immutable()
                self._inverse_flag_bases.append(MT)

    def compute_products(self):
        r"""
        Computes the products of the flags. This method is by default called from
        ``generate_flags``, and so would normally not need to be invoked directly.
        """
        self.state("compute_products", "yes")

        num_types = len(self._types)
        graph_block = make_graph_block(self._graphs, self._n)
        self._product_densities_arrays = []

        #sys.stdout.write("Computing products")
        print("Computing products...")

        from tqdm import tqdm
        import multiprocessing as mp
        from copy import deepcopy
        import math
        
        # print("Applying pool to "+str(num_types)+" types in parallel")
        
        arguments = []
        for ti in range(num_types):
            arguments.append( (self._types[ti], self._flags[ti], self._n, self._flag_cls, self._graphs) )
        
        # print("Using "+str(mp.cpu_count())+" cores")
        
        p = mp.Pool()
        for rarray in p.starmap(process_products_mp, tqdm(arguments)):
            self._product_densities_arrays.append(rarray)
        p.close()
        
        
        # for ti in tqdm(range(num_types)):
        #     
        #     tg = self._types[ti]
        #     s = tg.n
        #     m = (self._n + s) / 2
        # 
        #     flags_block = make_graph_block(self._flags[ti], m)
        #     
        #     rarray = self._flag_cls.flag_products(graph_block, tg, flags_block, None)
        #     self._product_densities_arrays.append(rarr
        
        # from tqdm import tqdm
        # 
        # for ti in tqdm(range(num_types)):
        # 
        #     tg = self._types[ti]
        #     s = tg.n
        #     m = (self._n + s) / 2
        # 
        #     flags_block = make_graph_block(self._flags[ti], m)
        #     rarray = self._flag_cls.flag_products(graph_block, tg, flags_block, None)
        #     self._product_densities_arrays.append(rarray)
        # 
        #     #sys.stdout.write(".")
        #     #sys.stdout.flush()
        # 
        # #sys.stdout.write("\n")

    def _set_block_matrix_structure(self):

        self.state("set_block_matrix_structure", "yes")

        self._block_matrix_structure = []

        for ti in self._active_types:

            num_blocks, block_sizes, block_offsets = 1, [len(self._flags[ti])], [0]

            # Remove zero-sized blocks
            bi = 0
            while bi < num_blocks:
                if block_sizes[bi] == 0:
                    num_blocks -= 1
                    del block_sizes[bi]
                    del block_offsets[bi]
                else:
                    bi += 1

            for bi in range(num_blocks):
                self._block_matrix_structure.append((ti, block_sizes[bi], block_offsets[bi]))

    def _get_block_matrix_structure(self, ti):

        num_blocks = 0
        block_indices = []
        block_offsets = []
        block_sizes = []

        for bi in range(len(self._block_matrix_structure)):
            b = self._block_matrix_structure[bi]
            if b[0] == ti:
                num_blocks += 1
                block_indices.append(bi)
                block_sizes.append(b[1])
                block_offsets.append(b[2])

        return num_blocks, block_sizes, block_offsets, block_indices
        
        
    def solve_sdp(self, show_output=False, solver="csdp",
        force_sharp_graphs=False, force_zero_eigenvectors=False,
        check_solution=True, tolerance=1e-5, show_sorted=False, show_all=False,
        use_initial_point=False, import_solution_file=None, csdp_settings=None):
        r"""
        Solves a semi-definite program to get a bound on the problem.

        INPUT:

         - ``show_output`` - Boolean (default: False). Whether to display output from the SDP
           solver.

         - ``solver`` - String (default: "csdp"). The SDP solver command to use. This can be one
           of the following:

            - "csdp" (Default) : the CSDP solver.
            - "sdpa" : the SDPA solver.
            - "sdpa_dd" : the double precision variant of SDPA
            - "sdpa_qd" : the quadruple precision variant of SDPA
            - "dsdp" : the DSDP solver.

            Note that the SDP solver must be present on the system; and it should be in a
            directory listed in PATH. The name of the solver should be "csdp", "sdpa", "sdpa_dd",
            "sdpa_qd" or "dsdp".

         - ``force_sharp_graphs`` - Boolean (default: False). If True, then the SDP is set up so
           that graphs that are supposed to be sharp are not given any "slack". Generally, this
           option is not particularly useful. It can sometimes improve the "quality" of a solution.

         - ``check_solution`` - Boolean (default: True). Whether to run ``check_solution`` to see
           if the bound appears to be tight; i.e. whether it is sufficiently close to the density
           given by the extremal construction.

         - ``tolerance`` - Number (default: 0.00001). This argument is passed to
           ``check_solution``. If a graph has a coefficient whose absolute difference with the
           bound is less than ``tolerance``, then it is considered to be sharp.

         - ``show_sorted`` - Boolean (default: False). This argument is passed to
           ``check_solution``. Whether to sort the sharp graphs according to their coefficients.
           If False, then they will be displayed in order of the graph indices.

          - ``show_all`` - Boolean (default: False). This argument is passed to
            ``check_solution``. If True, then the coefficients of all the graphs will be displayed
            instead of just the sharp graphs. In this case, the graphs that appear to be sharp are
            annotated with an "S" and the graphs that are forced to be sharp by the construction
            are annotated with a "C". (If these sets are not identical, then there is a problem.)

          - ``use_initial_point`` - Boolean (default: False). Whether to write an initial point
            file for the SDP solver. The initial point file is not used unless the solver is CSDP.
            Using an initial point can speed up the computation, but occasionally causes problems;
            so this option is False by default.

          - ``import_solution_file`` - Filename or None (default: None). If not None, then the SDP
            solver will not be run; instead the output file from a previous run of an SDP solver
            will be read. Care should be taken to ensure that the file being imported is for
            exactly the same problem, as minimal sanity-checking is done.
        """

        if import_solution_file is None:

            if solver=='csdp':
                self.write_csdp_settings_file(csdp_settings)
            
            if self.state("write_sdp_input_file") != "yes":
                self.write_sdp_input_file(force_sharp_graphs=force_sharp_graphs,
                                          force_zero_eigenvectors=force_zero_eigenvectors)
            if use_initial_point and self.state("write_sdp_initial_point_file") != "yes":
                self.write_sdp_initial_point_file()
            self._run_sdp_solver(show_output=show_output, solver=solver,
                                 use_initial_point=use_initial_point)

        else:

            self._sdp_output_filename = import_solution_file
            # pretend we have run the solver!
            self.state("run_sdp_solver", "force_yes")

        self._read_sdp_output_file()

        if check_solution:
            self.check_solution(tolerance=tolerance, show_sorted=show_sorted, show_all=show_all)

            
    def write_csdp_settings_file(self, settings):
        
        if settings is None:
            settings = {
                'axtol': 1.0e-8,
                'atytol': 1.0e-8,
                'objtol': 1.0e-8,
                'pinftol': 1.0e8,
                'dinftol': 1.0e8,
                'maxiter': 10000,
                'minstepfrac': 0.90,
                'maxstepfrac': 0.97,
                'minstepp': 1.0e-8,
                'minstepd': 1.0e-8,
                'usexzgap': 1,
                'tweakgap': 0,
                'affine': 0,
                'printlevel': 1,
                'perturbobj': 1,
                'fastmode': 0
            }
            
        self._csdp_settings_filename = os.path.join(str(SAGE_TMP), "param.csdp")
        
        if os.path.exists(self._csdp_settings_filename):
            sys.stdout.write("Deleting existing CSDP settings file...\n")
            os.remove(self._csdp_settings_filename)
        
        sys.stdout.write("Writing CSDP settings file...\n")
        
        with open(self._csdp_settings_filename, "w") as f:
            for k,v in settings.items():
                assert k in [
                    'axtol', 'atytol', 'objtol', 'pinftol', 'dinftol', 'maxiter', 'minstepfrac',
                    'maxstepfrac', 'minstepp', 'minstepd', 'usexzgap',  'tweakgap',  'affine',    
                    'printlevel', 'perturbobj', 'fastmode', 
                ]
                f.write(f"{k}={v}\n")
            
    # TODO: add option for forcing sharps

    def write_sdp_input_file(self, force_sharp_graphs=False, force_zero_eigenvectors=False):
        r"""
        Writes an input file for the SDP solver, specifying the SDP to be solved. This method is
        by default called by ``solve_sdp``.

        INPUT:

         - ``force_sharp_graphs`` - Boolean (default: False). If True, then the SDP is set up so
           that graphs that are supposed to be sharp are not given any "slack". Generally, this
           option is not particularly useful. It can sometimes improve the "quality" of a solution.
        """
        num_graphs = len(self._graphs)
        num_types = len(self._types)
        num_active_densities = len(self._active_densities)
        num_density_coeff_blocks = len(self._density_coeff_blocks)

        if num_active_densities < 1:
            raise NotImplementedError("there must be at least one active density.")

        if num_density_coeff_blocks < 1:
            raise NotImplementedError("there must be at least one density coefficient block.")

        if self.state("set_block_matrix_structure") != "yes":
            self._set_block_matrix_structure()
        total_num_blocks = len(self._block_matrix_structure)

        if force_zero_eigenvectors:
            num_extra_matrices = sum(self._zero_eigenvectors[ti].nrows() for ti in self._active_types)
        else:
            num_extra_matrices = 0

        self.state("write_sdp_input_file", "yes")

        self._sdp_input_filename = os.path.join(str(SAGE_TMP), "sdp.dat-s")

        sys.stdout.write("Writing SDP input file...\n")

        with open(self._sdp_input_filename, "w") as f:

            # num constraints
            f.write("%d\n" % (num_graphs + num_density_coeff_blocks + num_extra_matrices,))

            # num blocks in each constraint
            f.write("%d\n" % (total_num_blocks + 3 + (1 if force_zero_eigenvectors else 0),))

            # block sizes
            f.write("1 ")
            for b in self._block_matrix_structure:
                f.write("%d " % b[1])

            f.write("-%d -%d" % (num_graphs, num_active_densities))
            if force_zero_eigenvectors:
                f.write(" -%d" % num_extra_matrices)
            f.write("\n")

            # RHS of the SDP problem
            f.write("0.0 " * num_graphs)
            f.write("1.0 " * num_density_coeff_blocks)
            f.write("0.0 " * num_extra_matrices)
            f.write("\n")

            # objective function (\delta)
            if not self._minimize:
                f.write("0 1 1 1 -1.0\n")
            else:
                f.write("0 1 1 1 1.0\n")

            if force_zero_eigenvectors:
                for mi in range(num_extra_matrices):
                    f.write("0 %d %d %d %s\n" % (total_num_blocks + 4, mi + 1, mi + 1, "1.0" if self._minimize else "-1.0"))

            # slack vars and bound c for each constraint
            for i in range(num_graphs):
                if not self._minimize:
                    f.write("%d 1 1 1 -1.0\n" % (i + 1,))
                else:
                    f.write("%d 1 1 1 1.0\n" % (i + 1,))
                # if not graph sharp, add buffer var to make constraint equality
                if not (force_sharp_graphs and i in self._sharp_graphs):
                    f.write("%d %d %d %d 1.0\n" % (i + 1, total_num_blocks + 2, i + 1, i + 1))

            # add objective function to the SDP
            for i in range(num_graphs):
                for j in range(num_active_densities):
                    d = self._densities[self._active_densities[j]][i]
                    if d != 0:
                        if self._minimize:
                            d *= -1
                        f.write("%d %d %d %d %s\n" % (i + 1, total_num_blocks + 3, j + 1, j + 1, d.n(digits=64)))

            # set constant equal to 1
            for i in range(num_density_coeff_blocks):
                for di in self._density_coeff_blocks[i]:
                    if di in self._active_densities:
                        j = self._active_densities.index(di)
                        f.write("%d %d %d %d 1.0\n" % (num_graphs + i + 1, total_num_blocks + 3, j + 1, j + 1))

            # fill block_matrix with entries stored in product_densities_arrays
            for ti in self._active_types:

                num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

                for row in self._product_densities_arrays[ti]:
                    gi = row[0]
                    j = row[1]
                    k = row[2]
                    bi = num_blocks - 1
                    if bi > 0:
                        while block_offsets[bi] > j:
                            bi -= 1
                        j -= block_offsets[bi]
                        k -= block_offsets[bi]
                    value = Integer(row[3]) / Integer(row[4])
                    f.write("%d %d %d %d %s\n" %
                            (gi + 1, block_indices[bi] + 2, j + 1, k + 1, value.n(digits=64)))

            # TODO: get working with blocks, inactive types
            if force_zero_eigenvectors:
                mi = 0
                for ti in self._active_types:
                    nf = len(self._flags[ti])
                    for zi in range(self._zero_eigenvectors[ti].nrows()):
                        for j in range(nf):
                            for k in range(j, nf):
                                value = self._zero_eigenvectors[ti][zi, j] * self._zero_eigenvectors[ti][zi, k]
                                if value != 0:
                                    f.write("%d %d %d %d %s\n" %
                                            (num_graphs + num_density_coeff_blocks + mi + 1, ti + 2, j + 1, k + 1, value.n(digits=64)))
                        f.write("%d %d %d %d -1.0\n" % (num_graphs + num_density_coeff_blocks + mi + 1, total_num_blocks + 4, mi + 1, mi + 1))
                        mi += 1
        
    def get_sdp(self):
        # Some old code
        num_graphs = len(self._graphs)
        num_types = len(self._types)
        num_active_densities = len(self._active_densities)
        num_density_coeff_blocks = len(self._density_coeff_blocks)

        if num_active_densities < 1:
            raise NotImplementedError("there must be at least one active density.")

        if num_density_coeff_blocks < 1:
            raise NotImplementedError("there must be at least one density coefficient block.")

        if self.state("set_block_matrix_structure") != "yes":
            self._set_block_matrix_structure()
        total_num_blocks = len(self._block_matrix_structure)

        num_extra_matrices = 0   

        self.state("write_sdp_input_file", "yes")
        
        # J is number of graphs
        J = num_graphs

        # I is number of types
        I = len(self._active_types)

        # M gives the matrix sizes, i.e., the number of flags per type
        M = []
        for ti in self._active_types:
            num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)
            assert num_blocks == len(block_sizes) == 1
            M.append(block_sizes[0])

        # c gives the target density
        assert num_active_densities == 1
        c = [self._densities[self._active_densities[0]][i] for i in range(num_graphs)]

        # C gives the pair flag denisities
        C = {i: {j: [] for j in range(J)} for i in range(I)}

        for ti in self._active_types:
            num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

            for row in self._product_densities_arrays[ti]:
                gi = row[0]
                j = row[1]
                k = row[2]
                bi = num_blocks - 1
                if bi > 0:
                    while block_offsets[bi] > j:
                        bi -= 1
                    j -= block_offsets[bi]
                    k -= block_offsets[bi]
                value = Integer(row[3]) / Integer(row[4])
                C[ti][gi].append((j, k, Integer(row[3]) / Integer(row[4])))
        
        return I, J, M, C, c

    # TODO: handle no sharp graphs

    def write_sdp_initial_point_file(self, small_change=1/Integer(10)):
        r"""
        Writes an initial point file for the SDP solver. The zero matrix gives a feasible primal
        point. If a construction has been set, then this will be used to generate a feasible dual
        point, otherwise a zero matrix will be used for this as well. Using this method can reduce
        the time it takes the SDP solver to find a solution - typically by a third.

        INPUT:

         - ``small_change`` - Number (default: 1/10). Both primal and dual points must be
           perturbed by adding a small positive amount to the leading diagonals, in order that the
           matrices have no zero eigenvalues. Smaller values are not necessarily better here. The
           optimal value seems to depend on the problem and solver used.
        """

        num_graphs = len(self._graphs)
        num_types = len(self._types)
        num_active_densities = len(self._active_densities)

        self._sdp_initial_point_filename = os.path.join(str(SAGE_TMP), "sdp.ini-s")

        if self.state("set_block_matrix_structure") != "yes":
            self._set_block_matrix_structure()
        total_num_blocks = len(self._block_matrix_structure)

        self.state("write_sdp_initial_point_file", "yes")

        sys.stdout.write("Writing SDP initial point file...\n")

        with open(self._sdp_initial_point_filename, "w") as f:

            if self.state("set_construction") == "yes":

                for gi in range(num_graphs):
                    if gi in self._sharp_graphs:
                        si = self._sharp_graphs.index(gi)
                        f.write("%s " % self._sharp_graph_densities[si].n(digits=64))
                    else:
                        f.write("0.0 ")

                if not self._minimize:
                    f.write("%s\n" % (-self._target_bound).n(digits=64))
                else:
                    f.write("%s\n" % self._target_bound.n(digits=64))

                f.write("1 1 1 1 %s\n" % small_change.n(digits=64))

                for ti in range(num_types):

                    nf = len(self._flags[ti])
                    z_matrix = matrix(self._field, nf, nf)

                    for row in self._product_densities_arrays[ti]:
                        gi = row[0]
                        if not gi in self._sharp_graphs:
                            continue
                        si = self._sharp_graphs.index(gi)
                        j = row[1]
                        k = row[2]
                        value = Integer(row[3]) / Integer(row[4])
                        z_matrix[j, k] += value * self._sharp_graph_densities[si]

                    for j in range(nf):
                        z_matrix[j, j] += small_change

                    num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

                    for bi in range(num_blocks):
                        for j in range(block_sizes[bi]):
                            for k in range(j, block_sizes[bi]):
                                value = z_matrix[block_offsets[bi] + j, block_offsets[bi] + k]
                                if value != 0:
                                    f.write("1 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, k + 1, value.n(digits=64)))

                for gi in range(num_graphs):
                    if gi in self._sharp_graphs:
                        si = self._sharp_graphs.index(gi)
                        value = self._sharp_graph_densities[si]
                    else:
                        value = 0
                    if value <= 0:
                        value = small_change
                    f.write("1 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, value.n(digits=64)))

                for j in range(num_active_densities):
                    f.write("1 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, small_change.n(digits=64)))

            else:

                for gi in range(num_graphs + 1):
                    f.write("0.0 ")
                f.write("\n")
                f.write("1 1 1 1 %s\n" % small_change.n(digits=64))
                for ti in range(num_types):
                    num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)
                    for bi in range(num_blocks):
                        for j in range(block_sizes[bi]):
                            f.write("1 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, j + 1, small_change.n(digits=64)))
                for gi in range(num_graphs):
                    f.write("1 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, small_change.n(digits=64)))
                for j in range(num_active_densities):
                    f.write("1 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, small_change.n(digits=64)))

            # TODO: make this an exact Q check.
            if self.state("check_exact") == "yes":

                f.write("2 1 1 1 %s\n" % self._bound.n(digits=64))

                for ti in range(num_types):

                    num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

                    for bi in range(num_blocks):
                        for j in range(block_sizes[bi]):
                            for k in range(j, block_sizes[bi]):
                                value = self._exact_Q_matrices[ti][block_offsets[bi] + j, block_offsets[bi] + k]
                                if j == k:
                                    value += small_change
                                if value != 0:
                                    f.write("2 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, k + 1, value.n(digits=64)))

                for gi in range(num_graphs):
                    value = self._bound - self._bounds[gi]
                    if value <= 0:
                        value = small_change
                    f.write("2 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, value.n(digits=64)))

                for j in range(num_active_densities):
                    value = self._exact_density_coeffs[self._active_densities[j]]
                    if value <= 0:
                        value = small_change
                    f.write("2 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, value.n(digits=64)))

            else:

                densities = [sum([self._densities[j][gi] for j in self._active_densities])
                             / num_active_densities for gi in range(num_graphs)]

                bound = min(densities) if self._minimize else max(densities)

                value = bound
                if value <= 0:
                    value = small_change
                f.write("2 1 1 1 %s\n" % value.n(digits=64))

                num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

                for ti in range(num_types):

                    num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

                    for bi in range(num_blocks):
                        for j in range(block_sizes[bi]):
                            f.write("2 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, j + 1, small_change.n(digits=64)))

                for gi in range(num_graphs):
                    value = (bound - densities[gi]) * (-1 if self._minimize else 1)
                    if value <= 0:
                        value = small_change
                    f.write("2 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, value.n(digits=64)))

                value = Integer(1) / num_active_densities
                for j in range(num_active_densities):
                    f.write("2 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, value.n(digits=64)))

    # TODO: report error if problem infeasible

    def _run_sdp_solver(self, show_output=False, solver="csdp", use_initial_point=False):

        self.state("run_sdp_solver", "yes")

        previous_directory = os.getcwd()
        os.chdir(str(SAGE_TMP))

        print("Now in directory "+str( os.getcwd() ))
        
        if solver == "csdp":
            cmd = "%s %s sdp.out" % (cdsp_cmd, self._sdp_input_filename)

            if use_initial_point and self.state("write_sdp_initial_point_file") == "yes":
                cmd += " %s" % self._sdp_initial_point_filename

        elif solver == "dsdp":
            cmd = "%s %s -gaptol 1e-18 -print 1 -save sdp.out" % (dsdp_cmd, self._sdp_input_filename)

        elif solver == "sdpa":
            cmd = "%s -ds %s -o sdpa.out" % (sdpa_cmd, self._sdp_input_filename)

        elif solver == "sdpa_dd":
            cmd = "%s -ds %s -o sdpa.out" % (sdpa_dd_cmd, self._sdp_input_filename)

        elif solver == "sdpa_qd":
            cmd = "%s -ds %s -o sdpa.out" % (sdpa_qd_cmd, self._sdp_input_filename)

        else:
            raise ValueError("unknown solver.")

        sys.stdout.write(f"Running SDP solver with command '{cmd}'...\n")

        # For maximization problems, the objective value returned by the SDP solver
        # must be negated.
        obj_value_factor = 1.0 if self._minimize else -1.0

        child = pexpect.spawn(cmd, timeout=2**40)
        obj_val = None
        self._sdp_solver_output = ""
        sys.stdout.write("Reading output file...\n")

        # import time
        # time.sleep(1)

        while True:
            if child.eof():
                break
            try:
                child.expect("\r\n")
                line = child.before.decode("utf-8").strip()  + "\n"
                self._sdp_solver_output += line

                if show_output:
                    sys.stdout.write(line)

                if "Primal objective value:" in line:  # CSDP
                    print(f"Updating objective value to {line.split()[-1]}")
                    obj_val = self._approximate_field(line.split()[-1]) * obj_value_factor
                elif "objValPrimal" in line:  # SDPA
                    obj_val = self._approximate_field(line.split()[-1]) * obj_value_factor
                elif "DSDP Solution" in line:  # DSDP: seems to print absolute value
                    obj_val = self._approximate_field(line.split()[-1])

            except OverflowError:
                continue

            except pexpect.EOF:
                break

        child.close()
        self._sdp_solver_returncode = child.exitstatus

        sys.stdout.write("Returncode is %d. Objective value is %s.\n" % (
            self._sdp_solver_returncode, obj_val))

        # TODO: if program is infeasible, a returncode of 1 is given,
        # and output contains "infeasible"

        if "sdpa" in solver:

            with open("sdpa.out", "r") as inf:
                with open("sdp.out", "w") as f:

                    found, diagonal = False, False
                    t, row, col = 0, 1, 1

                    for line in inf:
                        if line[:6] == "yMat =":
                            break
                    else:
                        raise ValueError

                    for line in inf:

                        if line == "}":
                            break
                        elif line[:3] == "{ {":
                            t += 1
                            row = 1
                            diagonal = False
                        elif line[:2] == "{+" or line[:2] == "{-":
                            t += 1
                            row = 1
                            diagonal = True

                        line = line.replace("{", "")
                        line = line.replace("}", "")
                        col = 1
                        for a in line.split(","):
                            try:
                                v = a.strip()
                                vf = float(v)  # only done to see if we get ValueError
                                if diagonal:
                                    f.write("2 %d %d %d %s\n" % (t, row, col, v))
                                    row += 1
                                elif row <= col:
                                    f.write("2 %d %d %d %s\n" % (t, row, col, v))
                                col += 1
                            except ValueError:
                                pass

                        if col > 1:  # at least one number found...
                            row += 1

        self._sdp_output_filename = os.path.join(str(SAGE_TMP), "sdp.out")
        os.chdir(previous_directory)

    # TODO: read in dual solution

    def _read_sdp_output_file(self):

        self.state("read_solution", "yes")

        num_graphs = len(self._graphs)
        num_types = len(self._types)
        num_densities = len(self._densities)

        if self.state("set_block_matrix_structure") != "yes":
            self._set_block_matrix_structure()
        num_blocks = len(self._block_matrix_structure)

        with open(self._sdp_output_filename, "r") as f:

            self._sdp_Q_matrices = [matrix(self._approximate_field,
                len(self._flags[ti]), len(self._flags[ti])) for ti in range(num_types)]

            self._sdp_density_coeffs = [self._approximate_field(0) for i in range(num_densities)]

            for line in f:
                numbers = line.split()
                if numbers[0] != "2":
                    continue
                bi = int(numbers[1]) - 2
                if bi == num_blocks + 1:
                    j = int(numbers[2]) - 1
                    di = self._active_densities[j]
                    self._sdp_density_coeffs[di] = self._approximate_field(numbers[4])
                    continue
                if bi < 0 or bi >= num_blocks:
                    continue
                j = int(numbers[2]) - 1
                k = int(numbers[3]) - 1
                ti, size, offset = self._block_matrix_structure[bi]
                j += offset
                k += offset

                self._sdp_Q_matrices[ti][j, k] = self._approximate_field(numbers[4])
                self._sdp_Q_matrices[ti][k, j] = self._sdp_Q_matrices[ti][j, k]

        for ti in range(num_types):
            self._sdp_Q_matrices[ti].set_immutable()

    def check_solution(self, tolerance=1e-5, show_sorted=False, show_all=False):
        r"""
        Checks the approximate floating point bound given by the SDP solver, and determines which
        graphs appear to be sharp. The apparently sharp graphs are compared against the graphs
        that the construction forces to be sharp.

        INPUT:

         - ``tolerance`` - Number (default: 0.00001) If a graph has a coefficient whose absolute
           difference with the bound is less than ``tolerance``, then it is considered to be sharp.

         - ``show_sorted`` - Boolean (default: False). Whether to sort the sharp graphs according
           to their coefficients. If False, then they will be displayed in order of the graph
           indices.

          - ``show_all`` - Boolean (default: False). If True, then the coefficients of all the
            graphs will be displayed instead of just the sharp graphs. In this case, the graphs
            that appear to be sharp are annotated with an "S" and the graphs that are forced to be
            sharp by the construction are annotated with a "C". (If these sets are not identical,
            then there is a problem.)
        """
        self.state("check_solution", "yes")

        num_types = len(self._types)
        num_graphs = len(self._graphs)
        num_densities = len(self._densities)

        sys.stdout.write("Checking numerical bound...\n")

        fbounds = [sum([self._densities[j][i] * self._sdp_density_coeffs[j] for j in range(num_densities)]) for i in range(num_graphs)]

        for ti in self._active_types:
            for row in self._product_densities_arrays[ti]:
                gi, j, k, numer, denom = row
                d = Integer(numer) / Integer(denom)
                value = self._sdp_Q_matrices[ti][j, k]
                if j != k:
                    d *= 2
                if not self._minimize:
                    fbounds[gi] += d * value
                else:
                    fbounds[gi] -= d * value

        if not self._minimize:
            bound = max(fbounds)
        else:
            bound = min(fbounds)

        self._sdp_bounds = fbounds

        if self.state("set_construction") == "yes":

            if abs(bound - self._approximate_field(self._target_bound)) < tolerance:
                sys.stdout.write("Bound of %s appears to have been met.\n" % self._target_bound)
            else:
                sys.stdout.write("Warning: bound of %s appears to have not been met.\n" % self._target_bound)
                return
            sharp_graphs = self._sharp_graphs
        else:
            sharp_graphs = []  # set dummy sharp_graphs

        apparently_sharp_graphs = [gi for gi in range(num_graphs) if abs(fbounds[gi] - bound) < tolerance]

        if show_sorted or show_all:

            if not self._minimize:
                sorted_indices = sorted(range(num_graphs), key=lambda i: -fbounds[i])
            else:
                sorted_indices = sorted(range(num_graphs), key=lambda i: fbounds[i])

            for gi in sorted_indices:
                if gi in apparently_sharp_graphs:
                    sys.stdout.write("S")
                elif not show_all:
                    break
                else:
                    sys.stdout.write(" ")
                if gi in self._sharp_graphs:
                    sys.stdout.write("C")
                else:
                    sys.stdout.write(" ")
                sys.stdout.write(" %s : graph %d (%s) " % (fbounds[gi], gi, self._graphs[gi]))
                sys.stdout.write("\n")

        if self.state("set_construction") != "yes":
            return

        if not (show_sorted or show_all):
            sys.stdout.write("The following %d graphs appear to be sharp:\n" % len(apparently_sharp_graphs))
            for gi in apparently_sharp_graphs:
                sys.stdout.write("%.12f : graph %d (%s)\n" % (fbounds[gi], gi, self._graphs[gi]))

        extra_sharp_graphs = [gi for gi in apparently_sharp_graphs if not gi in self._sharp_graphs]
        missing_sharp_graphs = [gi for gi in self._sharp_graphs if not gi in apparently_sharp_graphs]

        if len(extra_sharp_graphs) > 0:
            sys.stdout.write("Warning: additional sharp graphs: %s\n" % (extra_sharp_graphs,))

        for gi in missing_sharp_graphs:
            sys.stdout.write("Warning: graph %d (%s) does not appear to be sharp.\n" % (gi, self._graphs[gi]))

    def import_solution(self, directory, complement=False):
        r"""
        Imports a solution found by Flagmatic 1.0 or 1.5.

        INPUT:

         - ``directory`` - the Flagmatic output directory, which must contain a flags.py file.

         - ``complement`` - Boolean (default: False). If True, then the solution will be assumed
           to be for the complementary problem. For example, if we are trying to minimize the
           density of k-cliques, the complementary problem is to minimize the density of
           independent sets of size k.
        """
        self.state("write_sdp_input_file", "yes")
        self.state("run_sdp_solver", "yes")

        num_graphs = len(self._graphs)
        num_types = len(self._types)
        num_densities = len(self._densities)

        if num_densities > 1:
            raise NotImplementedError

        sys.path.insert(0, directory)
        dont_write_bytecode = sys.dont_write_bytecode
        sys.dont_write_bytecode = True

        try:
            import flags
        except ImportError:
            sys.stdout.write("Cannot find flags.py in directory provided.\n")
            return

        # TODO: check admissible graphs, target bound, (sharps?).

        if flags.n != self._n:
            raise ValueError

        if self._flag_cls().oriented:
            if not ("oriented %d-graph" % self._flag_cls().r) in flags.description:
                raise ValueError
        else:
            if not ("%d-graph" % self._flag_cls().r) in flags.description:
                raise ValueError

        if flags.num_H != num_graphs:
            raise ValueError

        if flags.num_types != num_types:
            raise ValueError

        type_translations = []
        flag_translations = []

        for ti in range(num_types):
            tg = self._flag_cls(flags.types[ti])
            if complement:
                tg = tg.complement(True)
            for tj in range(num_types):
                if tg.is_labelled_isomorphic(self._types[tj]):
                    type_translations.append(tj)
                    break
            else:
                raise ValueError

            num_flags = len(self._flags[tj])
            ftr = []

            for fi in range(num_flags):
                fg = self._flag_cls(flags.flags[ti][fi])
                fg.t = tg.n
                if complement:
                    fg = fg.complement(True)
                for fj in range(num_flags):
                    if fg.is_labelled_isomorphic(self._flags[tj][fj]):
                        ftr.append(fj)
                        break
                else:
                    raise ValueError("solution has a flag that is not present.")

            flag_translations.append(ftr)

        self._sdp_Q_matrices = [matrix(self._approximate_field, len(self._flags[ti]),
                                len(self._flags[ti])) for ti in range(num_types)]

        try:
            f = open(directory + "/" + flags.out_filename, "r")
        except IOError:
            try:
                f = gzip.open(directory + "/" + flags.out_filename + ".gz", "rb")
            except IOError:
                print("Could not open %s or %s.gz" % (flags.out_filename, flags.out_filename))
                return

        for line in f:
            numbers = line.split()
            if numbers[0] != "2":
                continue
            ti = int(numbers[1]) - 2
            if ti >= 0 and ti < num_types:
                tj = type_translations[ti]
                if tj in self._active_types:
                    j = flag_translations[ti][int(numbers[2]) - 1]
                    k = flag_translations[ti][int(numbers[3]) - 1]
                    self._sdp_Q_matrices[tj][j, k] = numbers[4]
                    self._sdp_Q_matrices[tj][k, j] = self._sdp_Q_matrices[tj][j, k]

        f.close()

        self._sdp_density_coeffs = [1.0]

        for ti in range(num_types):
            self._sdp_Q_matrices[ti].set_immutable()

        sys.path.remove(directory)
        sys.dont_write_bytecode = dont_write_bytecode

    def make_exact(self, denominator=1024, meet_target_bound=True,
                   protect=None, use_densities=True, use_blocks=True, rank=None, show_changes=False,
                   check_exact_bound=True, diagonalize=True):
        r"""
        Makes an exact bound for the problem using the approximate floating point bound
        found by the SDP solver.

        INPUT:

         - ``denominator`` - Integer (default: 1024). The denominator to use when
            rounding. Higher numbers will cause the solution to be perturbed less, but can
            cause the resulting certificates to be larger.

         - ``meet_target_bound`` - Boolean (default: True). Determines whether the
           solution should be coerced to meet the target bound. If this is False, then a
           simpler method of rounding will be used. If there is no target bound, then this
           option will be disregarded, and the simpler method of rounding will be used.

         - ``protect`` - Array of integers or None (default: None). If an array of
            integers is given, then the entries of the matrices for the types with these
            indices will not be adjusted. (Using this option is not recommended.)

          - ``use_densities`` - Boolean (default: True). Whether to adjust the density
            coefficients to meet the target bound. This option will be disregarded if
            there is only one density.

          - ``use_blocks`` - Boolean (default: True). When computing the new basis for
             the solution's Q matrices, determines whether to give them a block structure
             with two blocks, using the invariant anti-invariant idea of Razborov.

          - ``rank`` - Integer or None (default: None). When computing the DR matrix,
             stop after ``rank`` columns have been found. This can save time in the
             case that the rank of the DR matrix is known (for example, from a previous
             run).

          - ``show_changes`` - Boolean (default: False). When meeting the target bound,
             display the changes being made to the matrix entries and the density
             coefficients.

          - ``check_exact_bound`` - Boolean (default: True). Whether to check the bound
             afterwards.

          - ``diagonalize`` - Boolean (default: True). Whether to diagonalize the Q
             matrices afterwards. If ``meet_target_bound`` is False, the Q matrices are
             always diagonalized.
        """

        if meet_target_bound and self.state("set_construction") != "yes":
            meet_target_bound = False
            sys.stdout.write("No target bound to meet.\n")

        if meet_target_bound:
            self.change_solution_bases(use_blocks=use_blocks)
            num_sharps = len(self._sharp_graphs)
        else:
            self._sdp_Qdash_matrices = self._sdp_Q_matrices
            # if non-tight, we don't have a separate diagonalization step
            self._exact_diagonal_matrices = []
            self._exact_r_matrices = []
            num_sharps = 0

        self.state("make_exact", "yes")

        num_types = len(self._types)
        num_graphs = len(self._graphs)
        num_densities = len(self._densities)

        if protect is None:
            protect = []

        q_sizes = [self._sdp_Qdash_matrices[ti].nrows() for ti in range(num_types)]

        def rationalize(f):
            return Integer(round(f * denominator)) / denominator

        sys.stdout.write("Rounding matrices")

        self._exact_Qdash_matrices = []

        for ti in range(num_types):

            if meet_target_bound:

                M = matrix(QQ, q_sizes[ti], q_sizes[ti], sparse=True)
                for j in range(q_sizes[ti]):
                    for k in range(j, q_sizes[ti]):
                        value = rationalize(self._sdp_Qdash_matrices[ti][j, k])
                        if value != 0:
                            M[j, k] = value
                            M[k, j] = value

            else:

                try:
                    LF = numpy.linalg.cholesky(self._sdp_Qdash_matrices[ti])
                    # TODO: Consider using this:
                    # LF = self._sdp_Qdash_matrices[ti].cholesky_decomposition()
                except numpy.linalg.linalg.LinAlgError:
                # except ValueError:
                    sys.stdout.write("Could not compute Cholesky decomposition for type %d.\n" % ti)
                    return
                L = matrix(QQ, q_sizes[ti], q_sizes[ti], sparse=True)
                for j in range(q_sizes[ti]):
                    for k in range(j + 1):  # only lower triangle
                        L[j, k] = rationalize(LF[j, k])
                L.set_immutable()
                M = L * L.T
                if not meet_target_bound:
                    D = identity_matrix(QQ, q_sizes[ti], sparse=True)
                    D.set_immutable()
                    self._exact_diagonal_matrices.append(D)
                    self._exact_r_matrices.append(L)
            
            row_div = self._sdp_Qdash_matrices[ti].subdivisions()[0]
            M.subdivide(row_div, row_div)
            self._exact_Qdash_matrices.append(matrix(self._field, M))
            sys.stdout.write(".")
            sys.stdout.flush()

        sys.stdout.write("\n")

        self._exact_density_coeffs = [rationalize(self._sdp_density_coeffs[di]) for di in range(num_densities)]

        # TODO: Make density coeffs in each block sum to 1. Ever necessary when >1 density?
        for db in self._density_coeff_blocks:
            if len(db) == 1:
                self._exact_density_coeffs[db[0]] = Integer(1)

        # Force all inactive densities to be 0 (they should be anyway).
        for j in range(num_densities):
            if not j in self._active_densities:
                self._exact_density_coeffs[j] = Integer(0)

        if meet_target_bound:

            self.state("meet_target_bound", "yes")

            triples = [(ti, j, k) for ti in self._active_types for j in range(q_sizes[ti])
                       for k in range(j, q_sizes[ti])]

            num_triples = len(triples)
            triples.sort()
            triple_to_index = dict((triples[i], i) for i in range(num_triples))

            R = matrix(self._field, num_sharps, num_triples, sparse=True)

            sys.stdout.write("Constructing R matrix")

            # TODO: only use triples that correspond to middle blocks.

            for ti in self._active_types:

                Ds = [matrix(QQ, len(self._flags[ti]), len(self._flags[ti]))
                      for si in range(num_sharps)]

                for row in self._product_densities_arrays[ti]:
                    gi = row[0]
                    if not gi in self._sharp_graphs:
                        continue
                    si = self._sharp_graphs.index(gi)
                    j = row[1]
                    k = row[2]
                    value = Integer(row[3]) / Integer(row[4])
                    Ds[si][j, k] = value
                    Ds[si][k, j] = value

                if self.state("transform_solution") == "yes":
                    B = self._inverse_flag_bases[ti]
                    for si in range(num_sharps):
                        Ds[si] = B.T * Ds[si] * B

                for si in range(num_sharps):
                    for j in range(q_sizes[ti]):
                        for k in range(j, q_sizes[ti]):
                            trip = (ti, j, k)
                            value = Ds[si][j, k]
                            if j != k:
                                value *= 2
                            if self._minimize:
                                value *= -1
                            R[si, triple_to_index[trip]] = value

                sys.stdout.write(".")
                sys.stdout.flush()
            sys.stdout.write("\n")

            density_cols_to_use = []
            DR = matrix(self._field, num_sharps, 0)  # sparsity harms performance too much here
            EDR = DR.T

            sys.stdout.write("Constructing DR matrix")

            # Only if there is more than one density
            if num_densities > 1 and use_densities:

                for j in self._active_densities:

                    if not rank is None and DR.ncols() == rank:
                        break

                    new_col = matrix(QQ, [[self._densities[j][gi]] for gi in self._sharp_graphs])
                    if new_col.is_zero():
                        continue
                    EDR = EDR.stack(new_col.T)
                    EDR.echelonize()
                    if EDR[-1, :].is_zero():
                        EDR = EDR[:-1, :]
                        sys.stdout.write("~")
                        sys.stdout.flush()
                        continue

                    DR = DR.augment(new_col)
                    density_cols_to_use.append(j)
                    sys.stdout.write(".")
                    sys.stdout.flush()

                sys.stdout.write("\n")
                sys.stdout.write("DR matrix (density part) has rank %d.\n" % DR.ncols())

            col_norms = {}
            for i in range(num_triples):
                n = sum(x**2 for x in R.column(i))
                if n != 0:
                    col_norms[i] = n

            # Use columns with greatest non-zero norm - change minus to plus to
            # use smallest columns (not sure which is best, or maybe middle?)

            cols_in_order = sorted(col_norms.keys(), key = lambda i : -col_norms[i])
            cols_to_use = []

            for i in cols_in_order:

                if not rank is None and DR.ncols() == rank:
                    break

                ti, j, k = triples[i]
                if ti in protect:  # don't use protected types
                    continue
                new_col = R[:, i: i + 1]
                if new_col.is_zero():
                    continue
                EDR = EDR.stack(new_col.T)
                EDR.echelonize()
                if EDR[-1, :].is_zero():
                    EDR = EDR[:-1, :]
                    sys.stdout.write("~")
                    sys.stdout.flush()
                    continue

                DR = DR.augment(new_col)
                cols_to_use.append(i)
                sys.stdout.write(".")
                sys.stdout.flush()

            sys.stdout.write("\n")
            sys.stdout.write("DR matrix has rank %d.\n" % DR.ncols())

            T = matrix(self._field, num_sharps, 1)

            for si in range(num_sharps):

                gi = self._sharp_graphs[si]
                T[si, 0] = self._target_bound

                for j in range(num_densities):
                    if not j in density_cols_to_use:
                        T[si, 0] -= self._exact_density_coeffs[j] * self._densities[j][gi]

                for i in range(num_triples):
                    if not i in cols_to_use:
                        ti, j, k = triples[i]
                        T[si, 0] -= self._exact_Qdash_matrices[ti][j, k] * R[si, i]

            FDR = matrix(self._field, DR)
            try:
                X = FDR.solve_right(T)
            except ValueError:
                if rank is None:
                    raise ValueError("could not meet bound.")
                else:
                    raise ValueError("could not meet bound (try increasing the value of ``rank``).")

            RX = matrix(self._approximate_field, X.nrows(), 1)

            for i in range(len(density_cols_to_use)):
                di = density_cols_to_use[i]
                RX[i, 0] = self._exact_density_coeffs[di]
                self._exact_density_coeffs[di] = X[i, 0]

            for i in range(len(density_cols_to_use), X.nrows()):
                ti, j, k = triples[cols_to_use[i - len(density_cols_to_use)]]
                RX[i, 0] = self._sdp_Qdash_matrices[ti][j, k]
                self._exact_Qdash_matrices[ti][j, k] = X[i, 0]
                self._exact_Qdash_matrices[ti][k, j] = X[i, 0]

            if show_changes:
                for i in range(X.nrows()):
                    sys.stdout.write("%.11s -> %.11s " % (RX[i,0], RDF(X[i,0])))
                    if i < len(density_cols_to_use):
                        sys.stdout.write("(density %d)\n" % density_cols_to_use[i])
                    else:
                        sys.stdout.write("(matrix %d, entry [%d, %d])\n" % triples[cols_to_use[i - len(density_cols_to_use)]])

        for ti in range(num_types):
            self._exact_Qdash_matrices[ti].set_immutable()

        if check_exact_bound:
            self.check_exact_bound()

    def check_exact_bound(self, diagonalize=True):
        r"""
        Usually called by ``make_exact``. If the solution was transformed, then computes
        the Q matrices from the Q' matrices. If the solution was adjusted to meet the
        target bound, the Q matrices are checked (numerically) for negative eigenvalues.
        In all cases the bound is checked.

        If ``diagonalize`` is set to True, then ``diagonalize`` will be called at the
        end.
        """
        num_types = len(self._types)
        num_graphs = len(self._graphs)
        num_densities = len(self._densities)

        if num_densities > 1:
            negative_densities = [j for j in range(num_densities) if self._exact_density_coeffs[j] < 0]
            if len(negative_densities) > 0:
                sys.stdout.write("Warning! Densities %s have negative coefficients, so the bound is not valid.\n" % negative_densities)
                return
            sys.stdout.write("All density coefficients are non-negative.\n")

        # If we didn't try to meet the target bound, then the method of rounding is_exact
        # guaranteed to produce positive-semidefinite matrices.
        if self.state("meet_target_bound") == "yes":
            negative_types = []
            very_small_types = []
            for ti in self._active_types:
                if self._exact_Qdash_matrices[ti].nrows() > 0:
                    eigvals = sorted(numpy.linalg.eigvalsh(self._exact_Qdash_matrices[ti]))
                    if eigvals[0] < 0.0:
                        negative_types.append(ti)
                    elif eigvals[0] < 1e-6:
                        very_small_types.append(ti)

            if len(negative_types) > 0:
                sys.stdout.write("Warning! Types %s have negative eigenvalues, so the bound is not valid.\n" % negative_types)
                return
            sys.stdout.write("All eigenvalues appear to be positive.\n")
            if len(very_small_types) > 0:
                sys.stdout.write("Types %s have very small eigenvalues (but this is probably OK).\n" % very_small_types)

        self.state("check_exact", "yes")

        self._exact_Q_matrices = []

        if self.state("transform_solution") == "yes":
            for ti in range(num_types):
                B = self._inverse_flag_bases[ti]
                M = B * self._exact_Qdash_matrices[ti] * B.T
                self._exact_Q_matrices.append(M)
        else:
            self._exact_Q_matrices = self._exact_Qdash_matrices

        bounds = [sum([self._densities[j][i] * self._exact_density_coeffs[j]
                  for j in range(num_densities)]) for i in range(num_graphs)]

        for ti in self._active_types:
            for row in self._product_densities_arrays[ti]:
                gi, j, k, numer, denom = row
                d = Integer(numer) / Integer(denom)
                value = self._exact_Q_matrices[ti][j, k]
                if j != k:
                    value *= 2
                if not self._minimize:
                    bounds[gi] += d * value
                else:
                    bounds[gi] -= d * value

        if self._field == QQ:
            if not self._minimize:
                bound = max(bounds)
            else:
                bound = min(bounds)
        else:
            # Sorting doesn't currently work for number fields with embeddings, so use float approximation.
            # TODO: Check if Sage 5.0 fixes this.
            if not self._minimize:
                bound = max(bounds, key=lambda x: float(x))
            else:
                bound = min(bounds, key=lambda x: float(x))

        self._bounds = bounds
        self._bound = bound

        if self.state("meet_target_bound") != "yes":
            sys.stdout.write("Bound is %s (%s).\n" % (bound, bound.n(digits=10)))
            return

        if self._field == QQ:
            if not self._minimize:
                violators = [gi for gi in range(num_graphs) if bounds[gi] > self._target_bound]
            else:
                violators = [gi for gi in range(num_graphs) if bounds[gi] < self._target_bound]
        else:
            if not self._minimize:
                violators = [gi for gi in range(num_graphs) if float(bounds[gi]) > float(self._target_bound)]
            else:
                violators = [gi for gi in range(num_graphs) if float(bounds[gi]) < float(self._target_bound)]

        sys.stdout.write("Bound of %s attained by:\n" % self._target_bound)
        unexpectedly_sharp = []
        for gi in range(num_graphs):
            if bounds[gi] == self._target_bound:
                sys.stdout.write("%s : graph %d (%s)\n" % (bounds[gi], gi, self._graphs[gi]))
                if not gi in self._sharp_graphs:
                    unexpectedly_sharp.append(gi)

        if len(unexpectedly_sharp) > 0:
            sys.stdout.write("Warning: the following graphs unexpectedly attain the bound: %s\n"
                             % ", ".join([str(gi) for gi in unexpectedly_sharp]))

        if len(violators) > 0:
            sys.stdout.write("Bound violated by:")
            for gi in violators:
                sys.stdout.write("%s : graph %d (%s)\n" % (bounds[gi], gi, self._graphs[gi]))

        if diagonalize:
            self.diagonalize()

    def make_exact_plus(self, limit_denominator=100, meet_target_bound=True,
                   protect=None, use_densities=True, use_blocks=True, rank=None, show_changes=False,
                   check_exact_bound=True, diagonalize=True):

        if meet_target_bound and self.state("set_construction") != "yes":
            meet_target_bound = False
            sys.stdout.write("No target bound to meet.\n")

        if meet_target_bound:
            self.change_solution_bases(use_blocks=use_blocks)
            num_sharps = len(self._sharp_graphs)
        else:
            self._sdp_Qdash_matrices = self._sdp_Q_matrices
            # if non-tight, we don't have a separate diagonalization step
            self._exact_diagonal_matrices = []
            self._exact_r_matrices = []
            num_sharps = 0

        self.state("make_exact", "yes")

        num_types = len(self._types)
        num_graphs = len(self._graphs)
        num_densities = len(self._densities)

        if protect is None:
            protect = []

        q_sizes = [self._sdp_Qdash_matrices[ti].nrows() for ti in range(num_types)]
        
        from fractions import Fraction
        
        def rationalize(f):
            frac = Fraction.from_float(float(f)).limit_denominator(int(limit_denominator))
            return Integer(frac.numerator) / Integer(frac.denominator)
            #return Integer(round(f * denominator)) / denominator

        sys.stdout.write("Rounding matrices")

        self._exact_Qdash_matrices = []

        for ti in range(num_types):

            if meet_target_bound:

                M = matrix(QQ, q_sizes[ti], q_sizes[ti], sparse=True)
                for j in range(q_sizes[ti]):
                    for k in range(j, q_sizes[ti]):
                        value = rationalize(self._sdp_Qdash_matrices[ti][j, k])
                        if value != 0:
                            M[j, k] = value
                            M[k, j] = value

            else:

                try:
                    LF = numpy.linalg.cholesky(self._sdp_Qdash_matrices[ti])
                    # TODO: Consider using this:
                    # LF = self._sdp_Qdash_matrices[ti].cholesky_decomposition()
                except numpy.linalg.linalg.LinAlgError:
                # except ValueError:
                    sys.stdout.write("Could not compute Cholesky decomposition for type %d.\n" % ti)
                    return
                L = matrix(QQ, q_sizes[ti], q_sizes[ti], sparse=True)
                for j in range(q_sizes[ti]):
                    for k in range(j + 1):  # only lower triangle
                        L[j, k] = rationalize(LF[j, k])
                L.set_immutable()
                M = L * L.T
                if not meet_target_bound:
                    D = identity_matrix(QQ, q_sizes[ti], sparse=True)
                    D.set_immutable()
                    self._exact_diagonal_matrices.append(D)
                    self._exact_r_matrices.append(L)
            
            row_div = self._sdp_Qdash_matrices[ti].subdivisions()[0]
            M.subdivide(row_div, row_div)
            self._exact_Qdash_matrices.append(matrix(self._field, M))
            sys.stdout.write(".")
            sys.stdout.flush()

        sys.stdout.write("\n")

        self._exact_density_coeffs = [rationalize(self._sdp_density_coeffs[di]) for di in range(num_densities)]

        # TODO: Make density coeffs in each block sum to 1. Ever necessary when >1 density?
        for db in self._density_coeff_blocks:
            if len(db) == 1:
                self._exact_density_coeffs[db[0]] = Integer(1)

        # Force all inactive densities to be 0 (they should be anyway).
        for j in range(num_densities):
            if not j in self._active_densities:
                self._exact_density_coeffs[j] = Integer(0)

        if meet_target_bound:

            self.state("meet_target_bound", "yes")

            triples = [(ti, j, k) for ti in self._active_types for j in range(q_sizes[ti])
                       for k in range(j, q_sizes[ti])]

            num_triples = len(triples)
            triples.sort()
            triple_to_index = dict((triples[i], i) for i in range(num_triples))

            R = matrix(self._field, num_sharps, num_triples, sparse=True)

            sys.stdout.write("Constructing R matrix")

            # TODO: only use triples that correspond to middle blocks.

            for ti in self._active_types:

                Ds = [matrix(QQ, len(self._flags[ti]), len(self._flags[ti]))
                      for si in range(num_sharps)]

                for row in self._product_densities_arrays[ti]:
                    gi = row[0]
                    if not gi in self._sharp_graphs:
                        continue
                    si = self._sharp_graphs.index(gi)
                    j = row[1]
                    k = row[2]
                    value = Integer(row[3]) / Integer(row[4])
                    Ds[si][j, k] = value
                    Ds[si][k, j] = value

                if self.state("transform_solution") == "yes":
                    B = self._inverse_flag_bases[ti]
                    for si in range(num_sharps):
                        Ds[si] = B.T * Ds[si] * B

                for si in range(num_sharps):
                    for j in range(q_sizes[ti]):
                        for k in range(j, q_sizes[ti]):
                            trip = (ti, j, k)
                            value = Ds[si][j, k]
                            if j != k:
                                value *= 2
                            if self._minimize:
                                value *= -1
                            R[si, triple_to_index[trip]] = value

                sys.stdout.write(".")
                sys.stdout.flush()
            sys.stdout.write("\n")

            density_cols_to_use = []
            DR = matrix(self._field, num_sharps, 0)  # sparsity harms performance too much here
            EDR = DR.T

            sys.stdout.write("Constructing DR matrix")

            # Only if there is more than one density
            if num_densities > 1 and use_densities:

                for j in self._active_densities:

                    if not rank is None and DR.ncols() == rank:
                        break

                    new_col = matrix(QQ, [[self._densities[j][gi]] for gi in self._sharp_graphs])
                    if new_col.is_zero():
                        continue
                    EDR = EDR.stack(new_col.T)
                    EDR.echelonize()
                    if EDR[-1, :].is_zero():
                        EDR = EDR[:-1, :]
                        sys.stdout.write("~")
                        sys.stdout.flush()
                        continue

                    DR = DR.augment(new_col)
                    density_cols_to_use.append(j)
                    sys.stdout.write(".")
                    sys.stdout.flush()

                sys.stdout.write("\n")
                sys.stdout.write("DR matrix (density part) has rank %d.\n" % DR.ncols())

            col_norms = {}
            for i in range(num_triples):
                n = sum(x**2 for x in R.column(i))
                if n != 0:
                    col_norms[i] = n

            # Use columns with greatest non-zero norm - change minus to plus to
            # use smallest columns (not sure which is best, or maybe middle?)

            cols_in_order = sorted(col_norms.keys(), key = lambda i : -col_norms[i])
            cols_to_use = []

            for i in cols_in_order:

                if not rank is None and DR.ncols() == rank:
                    break

                ti, j, k = triples[i]
                if ti in protect:  # don't use protected types
                    continue
                new_col = R[:, i: i + 1]
                if new_col.is_zero():
                    continue
                EDR = EDR.stack(new_col.T)
                EDR.echelonize()
                if EDR[-1, :].is_zero():
                    EDR = EDR[:-1, :]
                    sys.stdout.write("~")
                    sys.stdout.flush()
                    continue

                DR = DR.augment(new_col)
                cols_to_use.append(i)
                sys.stdout.write(".")
                sys.stdout.flush()

            sys.stdout.write("\n")
            sys.stdout.write("DR matrix has rank %d.\n" % DR.ncols())

            T = matrix(self._field, num_sharps, 1)

            for si in range(num_sharps):

                gi = self._sharp_graphs[si]
                T[si, 0] = self._target_bound

                for j in range(num_densities):
                    if not j in density_cols_to_use:
                        T[si, 0] -= self._exact_density_coeffs[j] * self._densities[j][gi]

                for i in range(num_triples):
                    if not i in cols_to_use:
                        ti, j, k = triples[i]
                        T[si, 0] -= self._exact_Qdash_matrices[ti][j, k] * R[si, i]

            FDR = matrix(self._field, DR)
            try:
                X = FDR.solve_right(T)
            except ValueError:
                if rank is None:
                    raise ValueError("could not meet bound.")
                else:
                    raise ValueError("could not meet bound (try increasing the value of ``rank``).")

            RX = matrix(self._approximate_field, X.nrows(), 1)

            for i in range(len(density_cols_to_use)):
                di = density_cols_to_use[i]
                RX[i, 0] = self._exact_density_coeffs[di]
                self._exact_density_coeffs[di] = X[i, 0]

            for i in range(len(density_cols_to_use), X.nrows()):
                ti, j, k = triples[cols_to_use[i - len(density_cols_to_use)]]
                RX[i, 0] = self._sdp_Qdash_matrices[ti][j, k]
                self._exact_Qdash_matrices[ti][j, k] = X[i, 0]
                self._exact_Qdash_matrices[ti][k, j] = X[i, 0]

            if show_changes:
                for i in range(X.nrows()):
                    sys.stdout.write("%.11s -> %.11s " % (RX[i,0], RDF(X[i,0])))
                    if i < len(density_cols_to_use):
                        sys.stdout.write("(density %d)\n" % density_cols_to_use[i])
                    else:
                        sys.stdout.write("(matrix %d, entry [%d, %d])\n" % triples[cols_to_use[i - len(density_cols_to_use)]])

        for ti in range(num_types):
            self._exact_Qdash_matrices[ti].set_immutable()

        if check_exact_bound:
            self.check_exact_bound()
            
    def diagonalize(self):
        r"""
        For each matrix Q, produces a matrix R and a diagonal matrix M such that
        Q = R * M * R.T, where R.T denotes the transpose of R. Usually called from
        ``make_exact``. Note that if the solution has not been adjusted to meet a target
        bound, a simpler method of rounding is performed, and diagonalization is done
        at the same time.
        """

        self.state("diagonalize", "yes")

        self._exact_diagonal_matrices = []
        self._exact_r_matrices = []

        sys.stdout.write("Diagonalizing")

        for ti in range(len(self._types)):
            R, M = LDLdecomposition(self._exact_Qdash_matrices[ti])
            self._exact_diagonal_matrices.append(M)
            if self.state("transform_solution") == "yes":
                R = self._inverse_flag_bases[ti] * R
            self._exact_r_matrices.append(R)
            sys.stdout.write(".")
            sys.stdout.flush()
        sys.stdout.write("\n")

        # Q can now be computed as Q = R * M * R.T

        sys.stdout.write("Verifying")

        for ti in range(len(self._types)):
            Q = self._exact_r_matrices[ti] * self._exact_diagonal_matrices[ti] * self._exact_r_matrices[ti].T
            if Q != self._exact_Q_matrices[ti]:
                raise ValueError  # TODO: choose appropriate error
            sys.stdout.write(".")
            sys.stdout.flush()
        sys.stdout.write("\n")

    def describe(self):
        r"""
        Returns a human-readable description of the problem. This is used by
        ``make_certificate``.
        """

        description = self._flag_cls.description() + "; "
        description += "minimize " if self._minimize else "maximize "
        description += self._describe_density()

        forbidden = []
        for g in self._forbidden_graphs:
            forbidden.append(str(g))
        for g in self._forbidden_induced_graphs:
            forbidden.append("induced %s" % g)
        for pair in self._forbidden_edge_numbers:
            forbidden.append("induced %d.%d" % pair)
        if len(forbidden) > 0:
            description += "; forbid " + ", ".join(forbidden)

        return description

    def _describe_density(self):

        if len(self._density_graphs) == 0:
            return ""
        elif len(self._density_coeff_blocks) == 1 and len(self._density_coeff_blocks[0]) == 1:
            density = self._density_graphs[0]
            if len(density) == 1 and density[0][1] == 1:
                return "%s density" % density[0][0]
            else:
                return "density expression %s" % self._flag_cls.format_combination(density)
        else:
            return "combination of quantum graphs"


    def show_large_densities(self, larger_than=1e-4):

        self.state("run_sdp_solver", "ensure_yes")

        num_densities = len(self._densities)

        densities_to_use = []
        for j in range(num_densities):
            if self._sdp_density_coeffs[j] > larger_than:
                densities_to_use.append(j)

        sys.stdout.write("Densities: %s\n" % (densities_to_use,))

        sys.stdout.write("Coefficients: %s\n" % ([self._sdp_density_coeffs[j] for j in densities_to_use],))

        sys.stdout.write("Other densities: %s\n" % ([di for di in range(num_densities) if not di in densities_to_use],))

    # def show_independent_densities(self):

    #     self.state("run_sdp_solver", "ensure_yes")
    
    #     num_sharps = len(self._sharp_graphs)
    #     num_densities = len(self._densities)
    
    #     densities_to_use = []
        
    #     if len(self._sdp_density_coeffs) > 0:
    #         density_indices = sorted(range(num_densities), key=lambda i: -self._sdp_density_coeffs[i])
    #     else:
    #         density_indices = range(num_densities)
        
    #     DR = matrix(self._field, 0, num_sharps, sparse=True)
    #     EDR = matrix(self._field, 0, num_sharps, sparse=True)
                
    #     sys.stdout.write("Constructing DR matrix")
        
    #     for j in density_indices:
    #         new_row = matrix(QQ, [[self._densities[j][gi] for gi in self._sharp_graphs]], sparse=True)
    #         if new_row.is_zero():
    #             continue
    #         try:
    #             X = EDR.solve_left(new_row)
    #             continue
    #         except ValueError:
    #             DR = DR.stack(new_row)
    #             EDR = EDR.stack(new_row)
    #             EDR.echelonize()
    #             densities_to_use.append(j)
    #             sys.stdout.write(".")
    #             sys.stdout.flush()
            
    #     sys.stdout.write("\n")
    #     sys.stdout.write("Rank is %d.\n" % DR.nrows())

    #     sys.stdout.write("Densities: %s\n" % (densities_to_use,))

    #     sys.stdout.write("Coefficients: %s\n" % ([self._sdp_density_coeffs[j] for j in densities_to_use],))

        
    def _augment_certificate(self, data):
        
        if len(self._assumptions) == 0:
            return
        
#         assumption_strings = []
#         for assumption in self._assumptions:
#             axs = []
#             for g, coeff in assumption[1]:
#                 if coeff == 1:
#                     axs.append(str(g))
#                 else:
#                     cs = str(coeff)
#                     if " " in cs:
#                         cs = "(%s)" % cs
#                     axs.append("%s*%s" % (cs, g))
#             assumption_strings.append("[%s] %s >= 0" % (assumption[0], " + ".join(axs)))

        assumption_strings = ["[%s] %s >= 0" %
            (assumption[0], self._flag_cls.format_combination(assumption[1]))
            for assumption in self._assumptions]
        
        data["assumptions"] = assumption_strings
        data["assumption_flags"] = self._assumption_flags
        
        data["admissible_graph_densities"] = self._densities
        data["density_coefficients"] = self._exact_density_coeffs


    def add_exact_matrices(self, types, qmatrices, rmatrices):
        r"""
        Adds Qdash and R matrices [possibly diagonal] (the ones such that R*Qdash*R^T = Q).
        Not possible to do if Qdash matrices already exist.
        Preferably, Qdash is diagonal, but not necessary -- if not diagonal, R will be
        multiplied by a matrix P such that PQdashP^T is diagonal: (RP)D(RP)^T = Q

        INPUT:
        - types: list of types in respective order to matrices, i.e. ["2:", "2:12"]
        - qmatrices: list of Qdash matrices in respective order to types
                    matrices are assumed to have rational entries; each Q is a rational
                    matrix.
        - rmatrices: list of R matrices; form as for Qdash matrices
        """

        if len(qmatrices) == 0:
            return

        if not(len(qmatrices) == len(types) and len(qmatrices) == len(rmatrices)):
            sys.stdout.write("Number of types must equal number of matrices.\n")
            raise ValueError

        if self._exact_Qdash_matrices or self._exact_diagonal_matrices:
            sys.stdout.write("Qdash matrices already exist. Cannot add new matrices.\n")
            raise ValueError
        
        self._exact_Qdash_matrices = [matrix() for t in types]
        
        for t in range(types):
            tp = GraphFlag(types[t])
            if tp in self.types:
                indx = (self.types).index(tp)
                self._exact_Qdash_matrices[indx] = qmatrices[t]
                self._exact_r_matrices[indx] = rmatrices[t]
                
                
    def write_manual_certificate(self, filename, bound):
        r"""
        Writes a certificate in a (hopefully) clear format, in JSON, to the file specified
        by ``filename``. For more information about the contents of the certificates, see
        the User's Guide.

        Does NOT require exact bound to be computed. ONLY saves problem information.

        INPUT:
        
        - filename: name of the file to store certificate 
        - bound: exact bound (at this point exact bound is not known to Flagmatic)
        """


        # TODO: check whether add_exact_matrices was called

        if self.state("write_sdp_input_file") != "yes":
            sys.stdout.write("Cannot write pre-certificate. Not enough information about the Problem.\n")
            raise ValueError

        self._bound = Rational(bound) # set bound manually

        def upper_triangular_matrix_to_list(M):
            return [list(M.row(i))[i:] for i in range(M.nrows())]

        def matrix_to_list(M):
            return [list(M.row(i)) for i in range(M.nrows())]

        data = {
            "description": self.describe(),
            "bound": self._bound, # needs to be passed manually (not yet available)
            "order_of_admissible_graphs": self._n,
            "number_of_admissible_graphs": len(self._graphs),
            "admissible_graphs": self._graphs,
            "number_of_types": len(self._types),
            "types": self._types,
            "numbers_of_flags": [len(L) for L in self._flags],
            "flags": self._flags,
            "qdash_matrices": [upper_triangular_matrix_to_list(M) for M in qdash_matrices],
            "r_matrices": [matrix_to_list(M) for M in r_matrices],
        }

        if len(self._density_graphs) == 1:
            data["admissible_graph_densities"] = self._densities[0]

        # Allow subclasses to add more things (like assumptions)
        self._augment_certificate(data)

        def default_handler(O):
            # Only output an int if it is less than 2^53 (to be valid JSON).
            if O in ZZ and O < 9007199254740992:
                return int(Integer(O))
            return repr(O)

        try:
            with open(filename, "w") as f:
                json.dump(data, f, indent=4, default=default_handler)
            sys.stdout.write("Written pre-certificate.\n")

        except IOError:
            sys.stdout.write("Cannot open file for writing.\n")
        

        
    def write_certificate(self, filename):
        r"""
        Writes a certificate in a (hopefully) clear format, in JSON, to the file specified
        by ``filename``. For more information about the contents of the certificates, see
        the User's Guide.

        Does require exact bound to be computed beforehand.
        """

        self.state("write_certificate", "yes")

        def upper_triangular_matrix_to_list(M):
            return [list(M.row(i))[i:] for i in range(M.nrows())]

        def matrix_to_list(M):
            return [list(M.row(i)) for i in range(M.nrows())]

        if self.state("meet_target_bound") != "yes" or self.state("diagonalize") == "yes":
            qdash_matrices = self._exact_diagonal_matrices
            r_matrices = self._exact_r_matrices
        else:
            qdash_matrices = self._exact_Qdash_matrices
            r_matrices = self._inverse_flag_bases

        data = {
            "description": self.describe(),
            "bound": self._bound,
            "order_of_admissible_graphs": self._n,
            "number_of_admissible_graphs": len(self._graphs),
            "admissible_graphs": self._graphs,
            "number_of_types": len(self._types),
            "types": self._types,
            "numbers_of_flags": [len(L) for L in self._flags],
            "flags": self._flags,
            "qdash_matrices": [upper_triangular_matrix_to_list(M) for M in qdash_matrices],
            "r_matrices": [matrix_to_list(M) for M in r_matrices],
        }

        if len(self._density_graphs) == 1:
            data["admissible_graph_densities"] = self._densities[0]

        # Allow subclasses to add more things
        self._augment_certificate(data)

        def default_handler(O):
            # Only output an int if it is less than 2^53 (to be valid JSON).
            if O in ZZ and O < 9007199254740992:
                return int(Integer(O))
            return repr(O)

        try:
            with open(filename, "w") as f:
                json.dump(data, f, indent=4, default=default_handler)
            sys.stdout.write("Written certificate.\n")

        except IOError:
            sys.stdout.write("Cannot open file for writing.\n")


    def verify_stability(self, tgraph, fgraph):

        """
        INPUT:
        - "tgraph": graph string of the type that we want to use in the proof
        - "fgraph": usually the minimum graph that will be blown up (F in PDF file)

        CLAIM 1: forbidding tgraph gives strictly worse bound
        CLAIM 2: Q_tgraph has rank 1 less than dimension
        CLAIM 3: unique embeddability
        CLAIM 4: every sharp graph of order N admits a strong homomorphism into F
        """

        # check if construction set
        # check if density met

        if self._stable:
            print("The problem is already stable! Nothing to verify.\n")
            return
        

        # start off modestly
        claim1 = False # forbidding type gives strictly worse bound
        claim2 = False # rk(Q) = dim(Q)-1
        claim3 = False # uniquely embeddable && different vertices attach differently to type
        claim4 = False # every sharp graphs strongly embeds into F

        
        # first check if the bound obtained is exact and sharp??
        if self.state("check_exact") != "yes":
            sys.stdout.write("Exact bound was not verified. Call make_exact() first.\n")
            return

        self.write_certificate("cert1.js")
        
        thebound = self._bound
        print("The problem is exact with a sharp bound of", str(thebound)+".\n")

        Tgraph = None
        try:
            Tgraph = GraphFlag(tgraph)
            tindex = (self.types).index(Tgraph) # CAREFULL! Tgraph may NOT be in self.types!!
            print("You selected the following type:", tgraph, "\n")
        except:
            raise ValueError

        # ---------- CLAIM 2 -----------
        # check the rank/dim of Q
        theQ = self._exact_Q_matrices[tindex]
        if theQ.dimensions()[0] - theQ.rank() == 1:
            print("dim(Q)-rank(Q) = 1")
            claim2 = True


        # ---------- CLAIM 3 -----------
        # every sharp graph of order N admits a strong homomorphism into F
        Fgraph = None
        try:
            Fgraph = GraphFlag(fgraph)
        except:
            raise ValueError

        # self._sharp_graphs stores subgraphs of construction that have non-zero density in it
        c_subgraphs = self._construction.subgraphs(self._n)
        num_embeddable = 0
        num_sharps = 0
        for i in len(self.graphs):
            if self._bounds[i] == self._bound: # gi sharp
                num_sharps += 1
                if i in self._sharp_graphs: # gi embeddable
                    num_embeddable += 1
        sys.stdout.write("num_embeddable = %d" %num_embeddable)
        if num_embeddable == num_sharps:
            print("slfkjdsfoakjpsdfoajskdf Every sharp graph is embeddable into a blowup of", fgraph+".\n")
            claim4 = True

        # ---------- CLAIM 4 -----------
        # check if type is uniquely embeddable into F
        claim3a = False # unique embeddability
        claim3b = False # distinct attachment

        # Tgraph --> Fblowup is a unique embedding if:
        # p(Tgraph,Fgraph)*(Fgraph.n choose Tgraph.n)*|Aut(Tgraph)|/|Aut(Fgraph)|
        dens = Fgraph.subgraph_density(Tgraph)
        sageT = Graph(Tgraph.n)
        for e in Tgraph.edges:
            sageT.add_edge(e)
        sageF = Graph(Fgraph.n)
        for e in Fgraph.edges:
            sageF.add_edge(e)
        Taut_group_order = sageT.automorphism_group().order()
        Faut_group_order = sageF.automorphism_group().order()
        if dens*binomial(Fgraph.n,Tgraph.n)*Taut_group_order/Faut_group_order == 1:
            claim3a = True
            

        # work towards claim3b:
        
        Flabelled = copy(Fgraph)
        Flabelled.t = Fgraph.n # now the copy is fully labelled
        perms = Permutations(range(1,Fgraph.n+1))

        # first retrieve how to embed Tgraph into Fgraph
        if claim3a:
            for perm in perms:
                Ffresh = copy(Flabelled)
                Ffresh.relabel(perm)
                is_embeddable = True
                for i in range(1,Tgraph.n+1):
                    for j in range(1,Tgraph.n+1):
                        e = [i,j]
                        e.sort()
                        e = tuple(e)
                        if (e in Tgraph.edges and e in Ffresh.edges) or (not(e in Tgraph.edges) and not(e in Ffresh.edges)):
                            continue
                        else:
                            is_embeddable = False # could break after this, improve once correct
                        
                if is_embeddable:
                    Flabelled = copy(Ffresh)
                    break
                
        # have our Ffresh with Tgraph as an induced subgraph on first Tgraph.n vertices of Ffresh
        # now check if attachments are all distinct
        restricted_neighbourhoods = [[] for x in range(Flabelled.n)]
        nn = Tgraph.n
        for x,y in Flabelled.edges:
            # if (x,y) edge has x in Tgraph
            if x < nn+1: restricted_neighbourhoods[y-1].append(x)
            # if (x,y) edge has y in Tgraph
            if y < nn+1: restricted_neighbourhoods[x-1].append(y)

        if len(restricted_neighbourhoods) == len(set(map(tuple,restricted_neighbourhoods) )):
            claim3b = True

        claim3 = claim3a and claim3b
        if claim3:
            print(tgraph, "is uniquely embeddable into", fgraph, "and different vertices of", fgraph, "attach differently to", tgraph+".\n")

        # --------- CLAIM 1 ----------
        self.forbid(tgraph)
        # reset states
        for state_name, state_value in self._states.items():
            if state_name == 'specify':
                self._states[state_name] = "yes"
            else:
                self._states[state_name] = "no"
        self.generate_flags(order=self._n)
        self.write_sdp_input_file()
        self.solve_sdp(import_solution_file=None)
        self.make_exact()

        if not self._minimize:
            if self._bound < thebound:
                print("Forbidding", Tgraph, "yields a bound of", self._bound, "which is strictly less than", str(thebound)+".")
                claim1 = True

        if self._minimize:
            if self._bound > thebound:
                print("Forbidding", Tgraph, "yields a bound of", self._bound, "which is strictly more than", str(thebound)+".")
                claim1 = True
            

        if claim1: print("Claim 1: verified.")
        else: print("Claim 1: NOT verified.")
        if claim2: print("Claim 2: verified.")
        else: print("Claim 2: NOT verified.")
        if claim3: print("Claim 3: verified.")
        else: print("Claim 3: NOT verified.")
        if claim4: print("Claim 4: verified.")
        else: print("Claim 4: NOT verified.")
                    
        if claim1 and claim2 and claim3 and claim4:
            self._stable = True
            print("\nOooh la la - early Christmas! The problem is stable!\n")
            self.write_certificate("cert2.js")
            print("Certificates written into 'cert1.js' and 'cert2.js'.")
        return
        


    def verify_robust_stability(self, tgraph, fgraph=None):

        """
        INPUT:
        - "tgraph": graph string of the type that we want to use in the proof
        - "fgraph": usually the minimum graph that will be blown up (B in the paper)

        ASSUMPTION 1.2: Forbidden graphs are not subgraphs of extremal construction.

        CLAIM 0:  Bound is sharp, exact, and construction is a blow-up (balanced, weights are equal in Flagmatic)
        CLAIM 1:  |tgraph| <= N-2, and forbidding tgraph gives strictly worse bound
        CLAIM 2a: There exists unique (up to automorph of tgraph in F) strong homomorphism f from tgraph to F
        CLAIM 2b: Every distinct x,y in V(F), neigh_F(x) \cap f(V(tgraph)) \neq neigh_F(y) \cap f(V(tgraph))
        CLAIM 3:  Every sharp graph of order N admits a strong homomorphism into F

        Thm 4.1:
        1. CLAIM 0
        2. CLAIM 1 & 2
        3. CLAIM 3
        Then the problem is robustly F-stable.

        NOTE:
        Fgraph, or F, is denoted by B in the paper.
        """

        print("Verifying robust stability...")

        
        # check if work needed
        if self._robustly_stable:
            print("The problem is already robustly stable! Nothing to verify.\n")
            return
        
        # start off modestly
        assumption12 = False # blowups of F are admissible
        claim0 = False # previous work OK, i.e. flag algebras problem OK and bound sharp
        claim1 = False # forbidding type gives strictly worse bound
        claim2 = False # uniquely embeddable && different vertices attach differently to type
        claim3 = False # every sharp graphs strongly embeds into F

        # preliminaries
        
        Fblow = self._construction
        if Fblow is None:
            raise ValueError("You first need to set extremal construction and it must be a Blowup Construction.")
        Fgraph = Fblow.graph
        forbidden = self._forbidden_graphs
        forbidden_induced = self._forbidden_induced_graphs

        Tgraph = None
        try:
            if self._flag_cls().r == 2:
                Tgraph = GraphFlag(tgraph)
                self.tgraph = Tgraph
                if not (fgraph is None):
                    Fgraph = GraphFlag(fgraph)
                    print("You selected the following type:", tgraph, "and the following F graph:", fgraph, "\n")
                else:
                    print("You selected the following type:", tgraph, "and F graph was taken from extremal construction:", Fgraph, "\n")
            elif self._flag_cls().r == 3:
                Tgraph = ThreeGraphFlag(tgraph)
                if not (fgraph is None):
                    Fgraph = ThreeGraphFlag(fgraph)
                    print("You selected the following type:", tgraph, "and the following F graph:", fgraph, "\n")
                else:
                    print("You selected the following type:", tgraph, "and F graph was taken from extremal construction:", Fgraph, "\n")
        except:
            raise ValueError


        
        # --------- ASSUMPTION 1.2 ----------
        # blowups of F are admissible
        # -----------------------------------

        if forbidden: # if any forbidden graphs
          for g in forbidden:
            if g in Fblow.subgraphs(g.n):
                break
            else:
                assumption12 = True
        else: # vacuously true
            assumption12 = True
                
        # ---------- CLAIM 0 -----------
        # claim0a: previous work OK, i.e. flag algebras problem OK and bound sharp
        # claim0b: construction is a blowup
        # claim0a and claim0b
        # ------------------------------
        claim0a = False
        claim0b = False
        
        # first check if the bound obtained is exact and sharp?
        # claim0a
        if not self.state("check_exact") == "yes":
            sys.stdout.write("Exact bound was not verified. Call make_exact() first.\n")
            return
        claim0a = True
        
        thebound = self._bound
        print("The problem is exact with a sharp bound of", str(thebound)+".\n")
        
        # check if the construction is a blowup (then we implicitly have the part ratios)
        # claim0b
        if not isinstance(self._construction, BlowupConstruction):
            print("The construction is not a blow-up construction.")
        else:
            claim0b = True

        claim0 = claim0a and claim0b
        


        # ---------- CLAIM 3 -----------
        # every sharp graph of order N admits a strong homomorphism into F
        # ------------------------------
        

        c_subgraphs = self._construction.subgraphs(self.n) # these are subgraphs of construction
        num_embeddable = 0
        num_sharp = 0
        for gi in range(len(self.graphs)):
            if self._bounds[gi] == self._bound:
                num_sharp += 1
                if self.graphs[gi] in c_subgraphs:
                    num_embeddable += 1
                else:
                    print("Not every sharp graph is embeddable in the blowup of", Fgraph, ". For example, ", self._graphs[gi], " isn't.\n")
                    print("Claim 3 is FALSE.")
                    break
        if num_embeddable == num_sharp:
            print("Every sharp graph is embeddable into a blowup of", Fgraph, ".\n")
            claim3 = True



        # ---------- CLAIM 2 -----------
        # claim2a: tgraph is uniquely embeddable
        # claim2b: different vertices of F attach differently to tgraph (when embedded in F)
        
        # ------------------------------
        claim2a = False # unique embeddability
        claim2b = False # distinct attachment

        otuples = Tuples(range(1,Fgraph.n+1), Tgraph.n)        
        coTgraph = Tgraph.complement()
        coFgraph = Fgraph.complement()
        sageF = Graph(Fgraph.n)

        if Fgraph.n > 0:

            sys.stdout.write("Verifying Claim 2...\n")
            
            for e in Fgraph.edges:
                sageF.add_edge(e)
            Faut_group_order = sageF.automorphism_group().order()
        
            strong_hom_count = 0
            for tpl in otuples:
                # for each map into F, check if it induces T
                edge_missing = False
                for edge in Tgraph.edges:
                    if edge_missing == True:
                        break
                    imedge1 = (tpl[edge[0]-1],tpl[edge[1]-1])
                    imedge2 = (tpl[edge[1]-1],tpl[edge[0]-1])
                    if imedge1 in Fgraph.edges or imedge2 in Fgraph.edges:
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
                    if imcoedge1 in coFgraph.edges or imcoedge2 in coFgraph.edges:
                        continue
                    else:
                        coedge_missing = True
                        break
                    
                if coedge_missing or edge_missing:
                    continue # this wasn't a strong hom embedding of T into F
                else:
                    strong_hom_count += 1
                    
                    
            claim2a = strong_hom_count == Faut_group_order            
            
            
            # work towards claim2b:
            # find how to strongly embed Tgraph into Fgraph == (find induced copy in Blowup(Fgraph) )
            
            claim2b_message = ""
        

            # this is quite big for brute force, just skip.
            
            NF = Fgraph.n
            NT = Tgraph.n
            Fgraph_vertex_set = range(1,NF+1)
            
            combs = Combinations(range(1,NF+1), NT)
            neighbourhoods_in_TinF = list()
            
            for comb in combs:
                
                TinF = Fgraph.degenerate_induced_subgraph(comb)
                
                if TinF == Tgraph:
                    comb_c = copy(Fgraph_vertex_set)
                    for x in comb:
                        comb_c.remove(x)
                        
                    for v in comb_c:
                        v_neighbourhood_in_TinF = list()
                        for u in comb:
                            if (u,v) in Fgraph.edges or (v,u) in Fgraph.edges:
                                v_neighbourhood_in_TinF.append(u)
                                neighbourhoods_in_TinF.append(v_neighbourhood_in_TinF)
                break

            num_neighbourhoods_in_TinF = len(neighbourhoods_in_TinF)
        
            for nbhd in neighbourhoods_in_TinF:
                nbhd.sort()
            
            
            all_pairs_are_ok = True
            for nbhd1_i in range(num_neighbourhoods_in_TinF-1):
                this_pair_is_ok = True # to begin with
                for nbhd2_i in range(nbhd1_i+1, num_neighbourhoods_in_TinF):
                    if neighbourhoods_in_TinF[nbhd1_i] == neighbourhoods_in_TinF[nbhd2_i]:
                        this_pair_is_ok = False
                        break
                if this_pair_is_ok == False:
                    all_pairs_ar_ok = False
                    break

            if all_pairs_are_ok == True:
                claim2b = True
                sys.stdout.write("Claim 2b is True. All vertices of F attach differently to T.\n")
                sys.stdout.flush()
            else:
                sys.stdout.write("Claim 2b is False.\n")
                sys.stdout.flush()

        
            claim2 = claim2a and claim2b
            if claim2:
                print(Tgraph, "is uniquely embeddable into", Fgraph, "and different vertices of", Fgraph, "attach differently to", str(Tgraph)+".\n")
            else:
                if not claim2a:
                    print("Claim 2a is False.\n")
                if not claim2b:
                    print("Claim 2b is False.\n")

        else: # if Fgraph is too large
            
            claim2b_message = "Fgraph is too big. Not verifying claim2b.\n"
            
        # ---------- CLAIM 1 -----------
        # forbidding Tgraph will decrease objective function
        # ------------------------------

        newproblem = None
        
        if Tgraph.n > 1:
            newproblem = GraphProblem(order=self._n,
                                       forbid_induced=self._forbidden_induced_graphs+[tgraph],
                                       forbid=self._forbidden_graphs,
                                       forbid_homomorphic_images=self._forbid_homomorphic_images,
                                       density=self._density_graphs[0],
                                       minimize=self._minimize,
                                       type_orders=self._type_orders,
                                       types=self._types_from_input,
                                       max_flags=self._max_flags,
                                       compute_products=self._compute_products,
                                       mode=self._mode)


            newproblem.solve_sdp(import_solution_file=None)
            newproblem.make_exact()
            

            if not self._minimize:
                if newproblem._bound < thebound:
                    print("Forbidding", Tgraph, "yields a bound of", newproblem._bound, "which is strictly less than", str(thebound)+".")
                    claim1 = True
                else:
                    print("Forbidding", Tgraph, "yields a bound of", newproblem._bound, "which is at least", str(thebound)+".")

            else:
                if newproblem._bound > thebound:
                    print("Forbidding", Tgraph, "yields a bound of", newproblem._bound, "which is strictly more than", str(thebound)+".")
                    claim1 = True
                else:
                    print("Forbidding", Tgraph, "yields a bound of", newproblem._bound, "which is at most", str(thebound)+".")

        else: #case: Tgraph.n <= 1
            if not self._minimize:
                if 0 < thebound:
                    print("Forbidding", Tgraph, "yields a bound of", 0, "which is strictly less than", str(thebound)+".")
                    claim1 = True
                else:
                    print("Forbidding", Tgraph, "yields a bound of", 0, "which is at least", str(thebound)+".")

            else: 
                if 1 > thebound:
                    print("Forbidding", Tgraph, "yields a bound of", 1, "which is strictly more than", str(thebound)+".")
                    claim1 = True
                else:
                    print("Forbidding", Tgraph, "yields a bound of", 0, "which is at most", str(thebound)+".")
                    


        

        print() # newline

        # ASSUMPTION 1.2
        if assumption12:
            print("Assumption 1.2 verified.")
        else:
            print("Assumption 1.2 NOT verified.")

        # CONDITION 1, THM 4.1
        if claim0:
            print("Condition 1 of Thm 4.1 verified.")
        else:
            print("Condition 1 of Thm 4.1 NOT verified.")

        # CONDITION 2, THM 4.1
        if claim1 and claim2:
            print("Condition 2 of Thm 4.1 verified.")
        elif claim1 and not claim2:
            print("Condition 2 of Thm 4.1 NOT verified.")
            print("Forbidding your type graph", Tgraph, "does yield a strictly worse bound, but the rest of Condition 2 in Thm 4.1 does NOT hold.")
        elif claim2 and not claim1:
            print("Condition 2 of Thm 4.1 NOT verified.")
            print("Forbidding your type graph", Tgraph, "does NOT yield a strictly worse bound, but the rest of Condition 2 in Thm 4.1 holds.")
        else:
            print("Condition 2 of Thm 4.1 NOT verified.")
            print("Forbidding your type graph", Tgraph, "does NOT yield a strictly worse bound, and the rest of Condition 2 in Thm 4.1 does NOT hold either.")

        # CONDITION 3, THM 4.1
        if claim3:
            print("Condition 3 of Thm 4.1 verified.")
        else:
            print("Condition 3 of Thm 4.1 NOT verified.")
                    
        if assumption12 and claim0 and claim1 and claim2 and claim3:
            self._robustly_stable = True
            if newproblem:
                newproblem.write_certificate("cert_robust_stab.js") # only writes certificate of the FA problem with forbidden Tgraph
            else: # newproblem is None
                print("Certificate is trivial, since everything is forbidden in a graph that avoids a 1-vertex subgraph.\n")
            self.write_certificate("cert_flag_alg.js")
            print("\n[OK] ROBUST STABILITY\n")
        return
    

    def verify_perfect_stability(self, fgraph=None):
        """Verify conditions of Thm 7.1 for Perfect Stability.

        Conditions:
        1. Assumption 2.1
        2. Robust fgraph-stability
        3. Assumption 5.1.1 (each forbidden graph is twin-free)

        and one of the following:

        4. rk(Q_tau) = dim(Q_tau)-1
        5. \lambda(Forb(\F \cup F)) < \lambda(Forb(\F)) (assuming \lambda is being maximized)

        
        1. and 2. are automatically met if problem is robustly stable (prerequisite to this method)

        
        CLAIM 1: all forbidden graphs are twin-free
        CLAIM 2: rk(Q_tau) = dim(Q_tau)-1
        CLAIM 3: forbidding B resulting in a smaller extremal value (assuming \lambda is being maximized)
        

        INPUT
        
        - fgraph: graph string of F graph that is blown up into the extremal construction. 
        If it is not provided as an argument, it will be taken from the construction.
        """
        
        print("Verifying perfect stability...")


        # nothing to do
        if self._perfectly_stable:
            print("The problem is already perfectly stable! Nothing to verify.\n")
            return

        # missing prerequisite
        if not self._robustly_stable:
            print("Please verify robust stability first (call verify_robust_stability() method).\n")
            return

        # start modestly
        claim1 = False
        claim2 = False
        claim3 = False
        
        # ----------- CLAIM 1 ------------
        # all forbidden graphs are twin-free
        # twin-free: no x,y have same neighbours
        # --------------------------------

        twins_exist = False
        twins_graph = None # forbidden graph F that has twins
        twins_vx = None # 2-tuple of x,y twins in F
        
        for g in self._forbidden_graphs:

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
                print("CLAIM 1 NOT verified.\n(At least one of the forbidden graphs is not twin-free.")
                print("Example:", str(twins_graph)+". Consider vertices", twins_vx[0], "and", str(twins_vx[1])+ ".)\n")
            else:
                print("CLAIM 1 verified. All forbidden graphs are twin-free.\n")
                claim1 = True
                


        # ----------- CLAIM 2 -----------
        # rk(Q_tau) = dim(Q_tau)-1
        # -------------------------------

        # note that tau = tgraph (just different notation)

        tau_index = None
        try:
            tau_index = self.types.index(self.tgraph)
            Q_tau = self._exact_Q_matrices[tau_index]
            Q_tau_rk = Q_tau.rank()
            if Q_tau.dimensions()[0]-1 == Q_tau_rk:
                claim2 = True
                print("CLAIM 2 verified. Q_tau has dimensions", Q_tau.dimensions(), "and rank", str(Q_tau_rk)+".\n")
            else:
                print("CLAIM 2 NOT verified. Matrix Q_tau of dimensions", Q_tau.dimensions(), "has rank", str(Q_tau_rk)+".\n")
                
        except Exception:
            print("Most probably, tgraph is not among the Flagmatic's types. Will attempt to verify Claim 3 now...\n")
            




        if not claim2: # only verify claim3 if claim2 failed.
            
            
            # ----------- CLAIM 3 ------------
            # Adding F to forbidden family reduces \lambda value
            # --------------------------------
            
            Fgraph = None
            if fgraph == None:
                Fgraph = self._construction.graph
            else:
                try:
                    Fgraph = GraphFlag(fgraph)
                except ValueError:
                    print("Ooops! F graph could not be initialized. Try providing a nicer one. More pink!")
                    
                    
            original_bound = self._bound
            perfstabproblem = GraphProblem(order=self._n,
                                           forbid_induced=self._forbidden_induced_graphs+[str(Fgraph)],
                                           forbid=self._forbidden_graphs,
                                           forbid_homomorphic_images=self._forbid_homomorphic_images,
                                           density=self._density_graphs[0],
                                           minimize=self._minimize,
                                           type_orders=self._type_orders,
                                           types=self._types_from_input,
                                           max_flags=self._max_flags,
                                           compute_products=self._compute_products,
                                           mode=self._mode)
            
            perfstabproblem.solve_sdp(solver="csdp")
            perfstabproblem.make_exact()
            
            new_bound = perfstabproblem._bound
            
            if self._minimize == False:
                if new_bound < original_bound:
                    print("The bound of", new_bound, "is strictly less than the original bound of", original_bound, ", OK.")
                    claim3 = True
                else:
                    print("The bound of", new_bound, "is NOT strictly less than the original bound of", original_bound, ".")
            else:
                if new_bound > original_bound:
                    print("The bound of", new_bound, "is strictly more than the original bound of", original_bound, ", OK.")
                    claim3 = True
                else:
                    print("The bound of", new_bound, "is NOT strictly more than the original bound of", original_bound,".")
                    
            if claim3 == True:
                print("CLAIM 3 verified: fgraph is \lambda-minimal.\n")
            else:
                print("CLAIM 3 NOT verified.\n(With "+str(Fgraph)+" forbidden, the bound is  not strictly smaller than the current bound.")
 

        if claim1 and (claim2 or claim3):
            self._perfectly_stable = True
            if not claim2:
                perfstabproblem.write_certificate("cert_perf_stab.js")
            print("[OK] PERFECT STABILITY")
         
                        
        """
        #==========================================================
        # construct matrix M
        #==========================================================
        
        C = self._construction
        B = C.graph
        B.make_minimal_isomorph()
        Bc = B.complement()
        tau = self.tgraph
        tau_index = self.types.index(tau)
        tau_flags = self.flags[tau_index]

        maps_tau_B = itertools.product(range(1,B.n+1), repeat=tau.n) # generates all (tau.n)-tuples with entries in [1,...,B.n]
        
        # for each map h in maps_tau_B
        # check if "(i,j) edge/nonedge in tau" implies "(h(i),h(j)) edge/nonedge in B"
        # if yes, add it to X
        X = list()
        for h in maps_tau_B:
            is_hom = True
            
            for e in tau.edges:
                he = [h[e[0]-1],h[e[1]-1]]
                he.sort()
                he = tuple(he)
                if not he in B.edges:
                    # h not homomorphism
                    is_hom = False
                    break
            if is_hom:
                tauc = tau.complement()
                for ec in tauc.edges:
                    hec = [h[ec[0]-1], h[ec[1]-1]]
                    hec.sort()
                    hec = tuple(hec)
                    if not hec in Bc.edges:
                        # h not homomorphism
                        is_hom = False
                        break
            if is_hom:
                X.append(h)
        
        # for each f in X, construct matrix Mf
        numflags = len(self.flags[tau_index])
        D = matrix(QQ, ones_matrix(1,B.n))
        Q_tau = self._sdp_Q_matrices[tau_index]

        count = 0
        for f_tuple in X:
            f = list(f_tuple)

            Mf = [[0 for x in range(B.n)] for y in range(numflags)]

            for j in range(numflags): # rows of Mf
                Fj = self.flags[tau_index][j]
                
                for w in range(1,B.n+1): # columns of Mf
                    fj = f+[w]
                    edges_match_so_far = True
                    for i in range(tau.n):
                        if not edges_match_so_far:
                            break
                        pairF = [i+1,tau.n+1]
                        pairB = [fj[i],w]
                        pairF.sort()
                        pairB.sort()
                        pairF = tuple(pairF)
                        pairB = tuple(pairB)

                        if pairB in B.edges and pairF in Fj.edges:
                            continue
                        elif not (pairB in B.edges) and not (pairF in Fj.edges):
                            continue
                        else:
                            if f_tuple == (4,3,1):
                                print("true C", pairF, pairB, w, Fj, B)
                            edges_match_so_far = False

                    if edges_match_so_far:
                        Mf[j][w-1] = 1

            Mf = matrix(Mf)

            Mf_dash = Mf.transpose()*self._exact_Q_matrices[tau_index]*Mf
            D = D.transpose().augment(Mf_dash.transpose()).transpose()


        b = [0 for x in range(D.dimensions()[0])]
        b[0] = 1
        b = vector(b)
        a = D.solve_right(b)
        print(a)
        """
                                
                        
                        
            

            
    
def ThreeGraphProblem(order=None, **kwargs):
    r"""
    Returns a Problem object, that will represent a Turán-type 3-graph problem. For help
    with Problem objects, enter

    sage: help(Problem)
    """
    return Problem(ThreeGraphFlag, order, **kwargs)


def GraphProblem(order=None, **kwargs):
    r"""
    Returns a Problem object, that will represent a Turán-type graph problem. For help
    with Problem objects, enter

    sage: help(Problem)
    """
    return Problem(GraphFlag, order, **kwargs)


def OrientedGraphProblem(order=None, **kwargs):
    r"""
    Returns a Problem object, that will represent a Turán-type oriented graph problem. For
    help with Problem objects, enter

    sage: help(Problem)
    """
    return Problem(OrientedGraphFlag, order, **kwargs)


def TwoMultigraphProblem(order=None, **kwargs):
    r"""
    Returns a Problem object, that will represent a Turán-type 2-multigraph problem. For
    help with Problem objects, enter

    sage: help(Problem)
    """
    return Problem(TwoMultigraphFlag, order, **kwargs)


def ThreeMultigraphProblem(order=None, **kwargs):
    r"""
    Returns a Problem object, that will represent a Turán-type 3-multigraph problem. For
    help with Problem objects, enter

    sage: help(Problem)
    """
    return Problem(ThreeMultigraphFlag, order, **kwargs)


def fpd(tp, flg1, flg2, grph):
    """
    Return the value of p(flg1,flg2;grph) if the given type is tp.
    Function name stands for 'flag product density'.
    
    INPUT:
    
    Use strings. Make sure order of graph is |flg1|+|flg2|-|tp|
    - tp # type on which both flags are based
    
    - flg1, flg2 # flags to be multiplied

    - grph # graph with respect to which density will be taken

    EXAMPLE:
    
    sage: fpd("2:", "3:13(2)", "3:23(2)", "4:1234")
    sage: 1/3
    """
    
    try:
        # read imput and convert it to flagmatic's default representation
        t = GraphFlag(tp)
        t.make_minimal_isomorph()

        f1 = GraphFlag(flg1)
        f1.make_minimal_isomorph()

        f2 = GraphFlag(flg2)
        f2.make_minimal_isomorph()

        gr = GraphFlag(grph)
        gr.make_minimal_isomorph()

    except ValueError:

        print("You are feeding unhealthy things to the function!")
    

    #not a 100% way to avoid name collision, but good enough...    
    the_most_ridiculous_name = GraphProblem() 
    
    gblock = make_graph_block([gr], gr.n)
    fblock1 = make_graph_block([f1], f1.n)
    fblock2 = make_graph_block([f2], f2.n)
    
    try:

        prod = the_most_ridiculous_name._flag_cls.flag_products(gblock, t, fblock1, fblock2)

    except ValueError:

        print("Something went wrong when multiplying flags. Make sure your flag is on the right type etc.")
        sys.exit(1)

    num_val = Rational((prod[0][3], prod[0][4]))

    return num_val


def fpds(tp, flg1, flg2, nn):
    """
    Return a linear combination of graphs on nn vertices that such
    that p(flg1,flg2; graph) > 0 for each term in the combination
    (here graph is on nn vertices).
    Function name stands for 'flag product densities'
    
    INPUT:
    
    Use strings. Make sure nn is at least |flg1|+|flg2|-|tp|
    
    - tp # type on which both flags are based
    
    - flg1, flg2 # flags to be multiplied

    - nn # order of the graph class

    EXAMPLE:
    
    sage: fpds("2:", "3:13(2)", "3:23(2)", 4)
    sage: [(4:1234, 1/3), (4:121324, 1/12)]
    """
    
    try:
        t = GraphFlag(tp)
        f1 = GraphFlag(flg1)
        f2 = GraphFlag(flg2)
        n_gr = int(nn)

        t.make_minimal_isomorph() # use the right flag representation
        f1.make_minimal_isomorph()    
        f2.make_minimal_isomorph()

    except ValueError:
        print("You are feeding unhealthy things to the function!")


    # check for correct value input...
    if t.n != f1.t or t.n != f2.t:
        raise TypeError("Your input is inconsistent (flags must be on the type that you specify).")
    
    if nn < f1.n + f2.n - t.n:
        raise TypeError("Your input violates |flg1|+|flg2|-|tp| <= nn.")

    the_most_ridiculous_name = GraphProblem()
    grphs = the_most_ridiculous_name._flag_cls.generate_graphs(n_gr)
    
    gblock = make_graph_block(grphs, n_gr)
    fblock1 = make_graph_block([f1], f1.n)
    fblock2 = make_graph_block([f2], f2.n)
    
    try:
        prod = the_most_ridiculous_name._flag_cls.flag_products(gblock, t, fblock1, fblock2)
    except ValueError:
        print("You are feeding unhealthy things to the function!")
        sys.exit(1)
    

    # write product as list of 2-tuples, each being (graph, coeff)
    lin_comb = list()
    for item in prod:
        lin_comb.append((grphs[item[0]], Rational((item[3],item[4]))))

    return lin_comb

def dens(graph, family_dimension):
    """
    Return the graph represented as a linear combination of graphs on
    n = family_dimension vertices.

    INPUT:
    
    - graph: graph to be represented as a lin combination. Any string
      representation will do as input.

    - family_dimension: order of graphs that form quantum graph
      representation of the input graph.

    EXAMPLE:
    
    sage: dens("3:121323", 4)
    sage: [(4:121323, 1/4), (4:12131423, 1/4), (4:1213142324, 1/2), (4:121314232434, 1)]
    """

    try:

        g = GraphFlag(graph)
        g.make_minimal_isomorph()
        
        n_gr = Integer(family_dimension)
        
        the_most_ridiculous_name2 = GraphProblem();
        grphs = the_most_ridiculous_name2._flag_cls.generate_graphs(n_gr)
        
    except ValueError:
        print("You are feeding unhealthy things to the function!")
        sys.exit(1)

    quantum_graph = list()

    for h in grphs:
        dgh = h.subgraph_density(g)
        if dgh > 0:
            quantum_graph.append((h, dgh))
            
    return quantum_graph


