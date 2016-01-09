spectrum = [40] + [2]*75+[-10]*19

# FIXME fix it.
def does_interlace(Y, X):
    for i in xrange(1,len(Y)+1):
        x1 = i-1
        y = i-1
        x2 = (len(X)-len(Y)+i)-1
        if not (X[x1] >= Y[y] and Y[y] >= X[x2]):
            return False
    return True            


def does_interlace_float(Y, X) :
    for i in xrange(1,len(Y)+1):
        x1 = i-1
        y = i-1
        x2 = (len(X)-len(Y)+i)-1
        if not ((X[x1]-Y[y] >= 0.00000001 or abs(X[x1]-Y[y]) < 1/10**11) and ( abs(Y[y]-X[x2]) < 1/10**11 or Y[y]-X[x2] >= 0.000001)):
            return False
    return True           


def spec(G):
    return sorted(G.eigenvalues(), reverse=True)

def partitionedAM(G, field=QQ):
    v = 95
    k = 40
    A = G.am()
    n = G.order()+1
    M = matrix(field,n,n)
    inbe = 0        
    for i in xrange(n-1):
        for j in xrange(i+1, n-1):
            M[i,j] = M[j,i] = A[i,j]
        M[i,n-1] = k-G.degree(i)
        inbe += M[i,n-1]
        M[n-1,i] = (k-G.degree(i))/(v-n+1)
    M[n-1,n-1] = 2*(k*v/2.0 - G.size()-inbe)/(v-n+1)
    return M

def isInterlaced(G):
    return does_interlace(spec(partitionedAM(G)), spectrum)

def isInterlacedFast(G):
    return does_interlace_float(spec(partitionedAM(G, field=RDF)), spectrum)

