load genericSRG.sage
from sys import argv
from itertools import combinations, product

def deg_right(G,v):

    k = G.subgraph([17..26]).degree(v)
    t = int(G.has_edge(v,11)) + int(G.has_edge(v,12)) 

    # Note that the graphs are generated so that x_0 is adjacent with 13,14,15 on K_4
    m = int(not G.has_edge(v,16)) 
    degree_right = 11-t-m+k
    
    return degree_right


def extendVertexBetween(G, v, toFix, cann, degs):

    ret = []
    degree_right = deg_right(G,v)
    
    for nbr in combinations([0..9], degree_right):
        H = G.copy()
        H.add_edges((v,i) for i in nbr)
        
        for u in degs:
            D = degs[u]
            d = H.subgraph([17..26]+[u]).degree(u)
            if len(toFix) + d < D:
                break
        else:
            s = H.canonical_label(partition=[[0..9],[10,11,12],[13,14,15,16],[17..26]]).graph6_string()

            if s not in cann:
                X = H.subgraph(set(H) - toFix)
                cann.add(s)
                X.relabel()
                if isInterlacedFast(X):
                    ret += [H]

    return ret        

# This function will always assume that the vertices x_0,x_1,x_2 = {10,11,12} 
# of both G and H are joined to K_4 in the same way.
def extend(G, H):
    
    H.relabel({i:i+17 for i in [0..9]})
    G.add_edges(H.edges())
   
    # We calculate the expected number of neighbors for the vertices of X_2^0 
    degs = {}
    for v in [0..9]:
        k = G.subgraph([0..9]).degree(v)
        m = len(set(G[v]).intersection(G[10]).intersection([13..16]))
        t = int(G.has_edge(v,11)) + int(G.has_edge(v,12))
        degs[v] = 7+k+m+t
    print degs        

    generated = [G]
    
    toFix = set([17..26])
    V = sorted([17..26], key=lambda x: deg_right(G,x))

    while toFix:
        v = V.pop()
        toFix.remove(v)
        cann = set()
        generated_tmp = []
        for G in generated:
            generated_tmp += extendVertexBetween(G, v, toFix, cann, degs)
        generated = generated_tmp
    return generated

L = []
for line in open(argv[1]):
    g, h = line.split(' ')
    # note to self G je right graph
    L += extend(Graph(g), Graph(h))

print 'Got ', len(L), 'graphs'
o = open(argv[1]+'.out','w')
for G in L:
    o.write(G.graph6_string() + '\n')
o.close()    
