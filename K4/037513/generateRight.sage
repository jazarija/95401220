load genericSRG.sage
from multiprocessing import Pool 
from itertools import combinations

# This code generates all graphs K_4 \cup X_3 \cup X_2^0

def extendVertex(G, v, toFix, cann):
    
    k = G.subgraph([0..9]).degree(v)

    t = int(G.has_edge(v,11)) + int(G.has_edge(v,12)) 

    ret = []

    for x,y in combinations([13..16], 2):
        if x != 16 and y != 16:
            m = 2
        else:
            m = 1

        if m+t+k > 3:
            continue

        H = G.copy()
        H.add_edge(v,x)
        H.add_edge(v,y)

        s = H.canonical_label(partition=[[0..9],[10,11,12],[13..16]]).graph6_string()
    
        if s not in cann:
            cann.add(s) 
            X = H.subgraph(set(H)-toFix)               
            X.relabel()
            if isInterlacedFast(X):
                ret+= [H]
    return ret  

# Vertices 8,9 and 10 are x_0, x_1 and x_2 while 11,12,13,14 are the vertices of K_4

def extend(G):

    G.add_vertices([10,11,12])
    G.add_edges(combinations([13..16], 2))

    G.add_edge(10,13)
    G.add_edge(10,14)
    G.add_edge(10,15)

    for v in  [0..9]: 
        G.add_edge(10,v)

    edges_to_X3 = [ [(11,13), (11,14), (11,15), (12,13), (12,14), (12,15)],
                    [(11,13), (11,14), (11,15), (12,13), (12,14), (12,16)],
                    [(11,13), (11,14), (11,16), (12,13), (12,15), (12,16)] ]

    edges_to_X2 = [[],[(9,11)],[(9,12)],[(9,11),(9,12)],[(8,11),(9,12)]]

    generated = []

    for edges1 in edges_to_X3:
        for edges2 in edges_to_X2:
            H = G.copy()
            H.add_edges(edges1+edges2)
            generated += [H]

    toFix = set([0..9])

    while toFix:
        
        v = toFix.pop()
        generated_tmp = []
        cann = set()
        for H in generated:
            generated_tmp += extendVertex(H, v, toFix, cann)

        generated = generated_tmp
        
    print 'We got plenty of graphs', len(generated)

    return generated

L = []

global cann
cann = set() 


# In this first part we construct all possible graphs X_2^0 which
# are labeled with the integers 0,..,9. The vertices 8 and 9 are the
# ones that are potentially adjacent to x_1 or x_2.    

for G in graphs.nauty_geng("-t -D2 8"):

    if not isInterlacedFast(G):
        continue

    G.add_vertices( [8,9] )

    for nbr1, nbr2 in combinations( list(combinations([0..7],2))+list(combinations([0..7],1))+[[]], 2):

        H = G.copy()
        H.add_edges( (8, el) for el in nbr1)
        H.add_edges( (9, el) for el in nbr2)

        if max(H.degree()) > 2:
            continue

        if not H.subgraph([0..8]).is_triangle_free():
            continue

        s1 = H.canonical_label(partition = [ [8,9], [0..7] ]).graph6_string()
        s2 = H.canonical_label(partition = [ [9], [0..8] ]).graph6_string()

        if s1 not in cann or s2 not in cann:
            cann.add(s1)
            cann.add(s2)

            if isInterlacedFast(H):
                L += [H]

        H2 = H.copy()
        H2.add_edge(8,9)

        if max(H2.degree()) > 2:
            continue

        s1 = H2.canonical_label(partition = [ [8,9], [0..7] ]).graph6_string()
        s2 = H2.canonical_label(partition = [ [9], [0..8] ]).graph6_string()

        if s1 not in cann or s2 not in cann:
            cann.add(s1) 
            cann.add(s2)

            if isInterlacedFast(H2):
                L += [H2]

print 'got' , len(L), 'right graphs'

cann = set()
    
L2 = []
p = Pool(8)

for el in p.imap(extend, L):
    L2 += el
print 'We got a final list of length', len(L2) 

o = open('right.g6','w')
for G in L2:
    o.write(G.graph6_string() + '\n')
o.close() 
