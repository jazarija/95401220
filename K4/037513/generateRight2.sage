from sys import argv
from itertools import combinations, product

def create_new(G):

    ret = []
    case = G.subgraph([16,10,11,12]).degree(16)

    if case == 0:
        permutat = Permutations([13..15])
        for perm in permutat:
            H = G.copy()
            H.relabel({[13..15][i]:perm[i] for i in [0..2]})
            ret += [H]
            H1 = G.copy()
            H1.relabel({11:12,12:11})
            ret += [H1]

    if case == 1:
        permutat = Permutations([13,14])
        for perm in permutat:
            H = G.copy()
            H.relabel({[13,14][i]:perm[i] for i in [0,1]})
            ret+=[H]


    if case == 2:
        permutat = Permutations([14,15])
        for perm in permutat:
            H = G.copy()
            H.relabel({[14,15][i]:perm[i] for i in [0,1]})
            H.relabel({[11,12][i]:(perm[i]-3) for i in [0,1]})
            ret+=[H]

    # print case, len(ret)
    return ret




L = []
L3=[]
for line in open("right.g6"):
    L3+=[Graph(line).graph6_string()]
    for g in create_new(Graph(line)):
        L += [g.graph6_string()]
L2=[]
print len(L)
L=list(set(L))
print len(L)
for g in L:
    if g not in L3:
        L2+=[g]
print 'Got ', len(L2), 'graphs'
o = open('right2.out','w')
for G in L2:
    o.write(G + '\n')
o.close()    
