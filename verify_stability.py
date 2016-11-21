"""
By Theorem 7.1, the following needs to be verified.

ROBUST STABILITY:
-----------------
CLAIM 0:
- Assumption 2.1:
   (i)   admissible graphs avoid all forbidden graphs (can't check that all such graphs are included in admissible)
   (ii)  Fgraph is on [m], and blowups of Fgraph are admissible

CLAIM 1:
- Assumption 5.1.1: Each forbidden graph is twin-free.

CLAIM 2:
- Theorem 7.1.1: The FA bound is tight (with extremal constr graph = balanced blowup of Fgraph)

CLAIM 3:
- Theorem 7.1.2:
   (i)   |Tgraph| <= N-2
   (ii)  lambda(admissible) > lambda(admissible minus those that contain Tgraph)
   (iii) There is a strong hom f: Tgraph --> Fgraph
   (iv)  Neighbourhoods of distinct vertices x1,x2 of Fgraph intersect image of Tgraph in Fgraph differently.

CLAIM 4:
- Theorem 7.1.3: Every sharp graph of order N admits a strong hom into Fgraph


PERFECT STABILITY:
------------------

CLAIM 5: Problem is robustly stable.

CLAIM 6: One of the following:
- Theorem 7.1.i:  rk(Q_tau) = dim(Q_tau)-1
- Theorem 7.1.ii: lambda(admissible not containing Tgraph) < lambda(admissible)

"""

# Call inspect_certificate.py certFA.py --stability
#
# - Claim 0
# - Claim 1
# - Claim 2


