load ../../genericSRG.sage

found = false

# 0..3 = K_4
# 4,5 = x_1,x_2
# 6 = x_0

graphs1 = [ Graph ( [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (4, 5)]),
            Graph( [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (3, 5), (4, 5)])]
for G in graphs1:
    for nbr in subsets([4,5]):
        H = G.copy()
        H.add_vertex(6)
        H.add_edges ( (6,v) for v in nbr)
        if isInterlaced(H):
            found = True
            break
if not found:
    print 'No interlacing candidate was found. The claim holds.'            
