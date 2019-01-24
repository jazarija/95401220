from itertools import product

right = [Graph(line) for line in open('right2.g6')]
upper = [Graph(line) for line in open('upper.g6')]
o = open('candidatePairs.g6', 'w')
for G,H in product(right, upper):
    if G.subgraph([10..16]).is_isomorphic(H.subgraph([10..16])):
        o.write(G.graph6_string() + ' ' + H.graph6_string() + '\n')
o.close()

