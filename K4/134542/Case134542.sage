load ../../genericSRG.sage
from itertools import product

global cann
cann = {}

# The graph we recive needs to have the vertices in [0..17]
# adjacent to the correct number of vertices in {x_0',x_0''} = [21,22]
def fixGraph(G):
    global cann
    ret = []

    toFix = set()

    for v in [0..17]:
        deg = G.subgraph([0..17]).degree(v)
        nbrnum = len( set(G[v]).intersection(set([18,19])) )

        if deg == 1 and nbrnum == 0:
            toFix.add(v)
        if deg == 0 and nbrnum == 1:
            toFix.add(v)
        if deg == 2 or (deg == 1 and nbrnum == 1)  or (deg == 0 and nbrnum == 2):
            G.add_edge(v, 21)
            G.add_edge(v, 22)

    X = G.subgraph(set(G) - toFix)
    X.relabel()
    if not isInterlacedFast(X):
        return []

    toFix = list(toFix)

    # At this point it remains to fix the vertices in toFix,
    # these are the vertices that can chose their neighbor in [21,22]
    for nbr in product( * [[21,22]]*len(toFix)):
        H = G.copy()
        H.add_edges( (nbr[i], toFix[i]) for i in xrange(len(toFix)) )
        s = H.canonical_label(partition=[[0..17],[18,19],[20],[21,22]]).graph6_string() #odstrani to
        if s not in cann:
            cann[s] = True
            if isInterlacedFast(H):


                # A final check is offfered by Lemma 9 saying that x_0',x_0''
                # must have 16-t neighbors in X_2^{\emptyset}, where t is the
                # number of its neighbors in \{x_1,x_2\}
                t = len( set(G[21]).intersection([18,19]) )
                if len(set(G[21]).intersection([0..17])) != 19-t:
                    continue

                t = len( set(G[22]).intersection([18,19]) )
                if len(set(G[22]).intersection([0..17])) != 19-t:
                    continue
                ret += [H]

    return ret


def extendVertex(G, v, toFix):
    global cann

    ret = []

    if v == 18 or v == 19:
        t = 3
    elif v == 21 or v == 22:
        t = 1
    else:
        t = 2

    for nbr in Combinations([23..26],t):
        H = G.copy()
        H.add_edges( (v, el) for el in nbr)
        s = H.canonical_label(partition=[[0..17],[18,19],[20],[21,22],[23..26]]).graph6_string()
        if s not in cann:
            cann[s] = True
            X = H.subgraph( set(H) - set(toFix) )
            X.relabel()
            if isInterlacedFast(X):
                ret += [H]
    return ret

# We obtain a graph with 23 vertices we add vertices [23,24,25,26] representing
# K_4. The vertices [18,19] then need 3 vertices in K_4, 21,22 need a singe vertex in K_4
# and the vertices 0..17 need 2 vertices in K_4
def extend(G):
    G.add_edges( Combinations( [23..26], 2) )
    generated = [G]
    toFix = [0..17] + [18,19] + [21,22]
    while toFix:
        v = toFix.pop()
        generated_tmp = []
        for G in generated:
            generated_tmp += extendVertex(G,v, toFix)
        generated = generated_tmp
        print len(generated)
    return generated




L = []
for G in graphs.nauty_geng("18 -D2"):
    G.add_edge(20,18) # 18,19 = x_1,x_2
    G.add_edge(20,19) # 20 = x_0

    # We have 4 cases to cover
        # 1. x_1,x_2 have no neibgbohrs in X_2^{\emptyset}
        # 2. Precisely one of x_1,x_2 has a neighbor in X_2
        # 3. They both have distinct neighbors in X_2
        # 4. They both have the same neighbor in X_2.

    # Case 1

    if isInterlacedFast(G):
        L += [G.copy()]

    # Case 2 - we assume x_1 has a neighbor in X_2^{-0} and this neigbor has degree at most 1.
    candNbr = [v for v in [0..17] if G.subgraph([0..17]).degree(v) <= 1]

    for nbr in candNbr:
        H = G.copy()
        H.add_edge(18, nbr)
        s = H.canonical_label(partition=[[0..17],[18,19],[20]]).graph6_string()
        if s not in cann:
            cann[s] = True
            if isInterlacedFast(H):
                L += [H]
    # Case 3
    for nbr1,nbr2 in Combinations(candNbr, 2):
        H = G.copy()
        H.add_edge(18, nbr1)
        H.add_edge(19, nbr2)
        s = H.canonical_label(partition=[[0..17],[18,19],[20]]).graph6_string()
        if s not in cann:
            cann[s] = True
            if isInterlacedFast(H):
                L += [H]

    candNbr = [v for v in [0..17] if G.subgraph([0..17]).degree(v) == 0]
    # Case 4
    for nbr in candNbr:
        H = G.copy()
        H.add_edge(18, nbr)
        H.add_edge(19, nbr)
        s = H.canonical_label(partition=[[0..17],[18,19],[20]]).graph6_string()
        if s not in cann:
            cann[s] = True
            if isInterlacedFast(H):
                L += [H]

print "We have obtained", len(L), "candidate graphs. We will add x_0', x_0'' now."

L2 = []
for G in L:
    G.add_edge(20, 21) # x_0', x_0'' = 21,22
    G.add_edge(20, 22)

# By Lemma 7 it follows, that x_1,x_2 are each adjacent to at least one of the vertices in [21,22]

    for nbr1,nbr2 in CartesianProduct( [ (21,),(22,), (21,22) ], [ (21,),(22,), (21,22) ] ):
        H = G.copy()
        H.add_edges( (19,el) for el in nbr1)
        H.add_edges( (18,el) for el in nbr2)

        # s = H.canonical_label(partition=[[0..17],[18,19],[20],[21,22]]).graph6_string() #to particijo bi bio potreba odstranit
        # if s not in cann:
        #     cann[s] = True
        L2 += [H]
print 'We got', len(L2), 'candidate graphs. We add the remaining edges now.'

L = []
c = 0
for G in L2:
    L += fixGraph(G)
print 'We got', len(L), 'candidate graphs.'

L3 = []
for G in L:
    L3 += extend(G)

o = open('candidates.out', 'w')
for G in L3:
    o.write(G.graph6_string() + '\n')
o.close()
