from sys import argv
load ../../genericSRG.sage

global cann
cann = {}

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
for line in open(argv[1]):
    G = Graph(line)
    L += extend(G)
o = open(argv[1]+'.out','w')
for G in L:
    o.write(G.graph6_string() + '\n')
