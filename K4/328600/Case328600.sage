load ../../genericSRG.sage
from itertools import combinations

global cann
cann = {}


def extendVertex(G, v, toFix):
    ret = []
    for nbr in combinations([23..26], 2):
        H = G.copy()
        H.add_edges((v, el) for el in nbr)
        s = H.canonical_label(partition=[[0..2],[3..22],[23..26]]).graph6_string()
        if s not in cann:
            cann[s] = True
            X = H.subgraph(set(H) - set(toFix))
            X.relabel()
            if isInterlacedFast(X):
                ret += [H]
    return ret


# def final_add(G):
#     ret = []
#     G.add_vertex(23)
#     G.add_edges([(23, 1), (23, 2)])
#     for nbr1 in combinations([19..22], 2):
#         H = G.copy()
#         H.add_edges((23, el) for el in nbr1)
#         for nbr2 in Combinations([3..18]):
#             H2 = H.copy()
#             H2.add_edges((23, el) for el in nbr2)
#             s = H2.canonical_label(partition=[[0..2],[3..18],[19..22],[23]]).graph6_string()
#             if s not in cann:
#                 cann[s] = True
#                 H2.relabel()
#                 if isInterlacedFast(H2):
#                     ret += [H2]
#     return ret


def extend(G):
    G.add_edges(combinations([23..26], 2))
    generated = [G]
    toFix = [3..22]
    while toFix:
        v = toFix.pop()
        generated_tmp = []
        for G in generated:
            generated_tmp += extendVertex(G, v, toFix)
        generated = generated_tmp
        print len(generated)
    return generated
    # final_generated = []
    # for H in generated:
    #     final_generated += final_add(H)
    # print len(final_generated)
    # return final_generated





G = Graph()
G.add_vertices([0, 1, 2])

w = open('candidates.txt', 'w')
for H in graphs.nauty_geng("20 -d2 -D2"):
    H.relabel({i: i + 3 for i in range(20)})
    H = H.union(G)
    H.add_edges([(x, y) for x in range(3, 23) for y in range(2)])
    print isInterlaced(H), H.order(), H.size()
    ret = extend(H)
    if len(ret) > 0:
        print min([I.order() - I.spectrum().count(2) for I in ret])
    for J in ret:
        w.write(J.canonical_label().graph6_string() + '\n')
w.close()
