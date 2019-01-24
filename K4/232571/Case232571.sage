load ../../genericSRG.sage
from itertools import product

global cann
cann = {}

# We have that
#
# 0..18 are the vertices in X_2^0, each needing 2 neighbors in [23..26] 
# 19,20 = {x_0,x_1] each not adjacent with any vertex in [23..26]
# 21 needing 3 neighbors in [23..26]
# 22 needing 1 neighbor in [23..26]    
def extendVertex(G, v, toFix):
    global cann

    ret = []

    if v == 19 or v == 20:
        # t = 0 
        return [G]

    elif v == 21:
        t = 3
    elif v == 22:
        t=1
    else:
        t = 2

    for nbr in Combinations([23..26],t):
        H = G.copy()
        H.add_edges( (v, el) for el in nbr)
        s = H.canonical_label(partition=[[0..18],[19,20],[21],[22],[23..26]]).graph6_string()
        if s not in cann:
            cann[s] = True
            X = H.subgraph (set(H) - set(toFix) )
            X.relabel()
            if isInterlacedFast(X):
                ret += [H]
    return ret     

def fixGraph(G):
    for v in [0..18]:
        deg = G.subgraph([0..18]).degree(v)
        if deg == 2:
            G.add_edge(v,22)
        if deg == 0 and not G.has_edge(v,21):
            return []
        if deg == 1 and G.has_edge(v, 21):
            G.add_edge(v, 22)

    if isInterlacedFast(G):
        return [G]

    return []

def extend(G):
    
    G.add_edges( Combinations( [23..26], 2) )

    generated = [G]
    toFix = [0..22] 
    while toFix:
        v = toFix.pop()
        generated_tmp = []
        for G in generated:
            generated_tmp += extendVertex(G,v, toFix)
        generated = generated_tmp
        # print len(generated)
    return generated


L = []
for G in graphs.nauty_geng("19 -D2"):

    for x in [0..18]:
        G.add_edge(x,19)
    # 19 = x_0, 20 = x_1

    G.add_vertex(20)

    if not isInterlacedFast(G):
        continue

    G.add_edge(21,20) 
    G.add_edge(21,19) 
    # 21 = x_3
    # 22 = x_1'
    G.add_edge(22,20)


    # we cover two options based on whether  x_3 ~ x_1' or not
    
    H = G.copy()
    I = H.copy()
    # if it doesn' exists
    L += [I]

    # if it exists
    H.add_edge(22,21)
        
    candNbr = [v for v in [0..18] if G.subgraph([0..18]).degree(v) <= 1]
    for nbr in candNbr:
        I=H.copy()
        I.add_edge(21,nbr)
        L += [I]


print "We have obtained", len(L), "candidate graphs. We will add edges between x_0', x_1 and X_0^1 now."


L2 = []
for G in L:
    L2 += fixGraph(G)

print 'We got', len(L2), 'candidate graphs. We add K_4 now'

L3=[]
for G in L2:
    L3+=extend(G)
print 'We got', len(L3), 'candidate graphs.'
