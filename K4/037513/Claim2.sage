load ../../genericSRG.sage

# The vertices of the clique are 0,1,2,3
# X_3 = {4,5,6}
G1 = Graph ( [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 3), (2, 4), (2, 5), (2, 6)])
G2 = Graph( [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (1, 2), (1, 3), (1, 4), (1, 5),(1, 6),(2, 3),(2, 4),(2, 5), (3, 6)])
G3 = Graph( [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 6), (3, 5), (3, 6)])

for el in subsets( [ (4,5), (4,6), (5,6)]):
    if len(el) == 0:
        continue
    for G in [G1,G2,G3]:
        H = G.copy()
        H.add_edges(e for e in el)
        if isInterlaced(H):
            print 'There are interlacing configurations.'
            exit(1)

print 'There is no way to introduce edges in X_3.' 
