load ../genericSRG.sage
global cann
cann = {}

def fixTrianglePair(G, t1, t2):
    
    global cann
    ret = []

    u,v = t1[0], t2[0]
    
    # First we determine how many neighbors do t1 and t2 share
    # on K_5. If 0 then each vertex of t1 has precisely two
    # neighbors in t2 and one neighbor otherwise.
    inter = len (set(G[u]).intersection(G[v]).intersection([0..4]))

    if inter == 0:
        for nbr in Permutations( t2 ):
            H = G.copy()
            for i in xrange(3): 
                H.add_edges( (t1[i],el) for el in t2 if el != nbr[i])
            s = H.canonical_label(partition=[ [0..4],[v for v in H if v>4]]).graph6_string()
            if s not in cann:
                cann[s] = True
                ret += [H]                

    elif inter == 1:
        for nbr in Permutations( t2 ):
            H = G.copy()
            for i in xrange(3):
                H.add_edge( t1[i], nbr[i])
            s = H.canonical_label(partition=[ [0..4],[v for v in H if v>4]]).graph6_string()
            if s not in cann:
                cann[s] = True
                ret += [H]       
    elif inter == 2:
        H = G.copy()
        ret=[H]   
    return ret       


def findConfigurations(edges_to_fix):

    L = [graphs.CompleteGraph(5)]
    triangles = []
    global cann
    cann={}

    for x,y in edges_to_fix:
        for G in L:
            vc = max(G)+1
            cur_t =  (vc,vc+1,vc+2)
            G.add_edges(Combinations(cur_t,2))
            for v in cur_t:
                G.add_edge(v,x)
                G.add_edge(v,y)
            L2 = [G]
            for t in triangles:
                L3 = []
                for G in L2:
                    L3 += fixTrianglePair(G, cur_t, t)
                L2 = L3
        
        triangles += [cur_t]            
        L = [X for X in L2 if isInterlacedFast(X)]  
        print 'Currently dealing with', len(L), 'graphs', set( G.spectrum().count(2) for G in L)
    return L

e1 = [ (0,1), (0,2), (0,3),(0,4),(1,2),(1,3),(1,4),(2,3) ]
e2 = [ (0,1), (0,2), (0,3),(0,4),(2,4),(1,3),(1,4),(2,3) ]

o=open("triangles","w")
for G in findConfigurations(e2) + findConfigurations(e1):
    o.write(G.graph6_string()+"\n")


