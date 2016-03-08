#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <math.h>

#include <assert.h>
#include <errno.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

/* Size of input graphs - note this code will only work for size smaller than 62 */
#define N 19

/* Size/valency of the original graph */
#define NTOT 95
#define VAL 40
#define LAMBDA 12 
#define MU  20

#define EIG_MAX 2 
#define EIG_MIN -10

#define DIAG (EIG_MAX + (double) (VAL-EIG_MAX)/NTOT)
#define ADJ (-1 + (double) (VAL-EIG_MAX)/NTOT)
#define NONADJ ((double) (VAL-EIG_MAX)/NTOT)
/* Forbidden eigenvalue of the star complement */
#define SC_EIG 2

/* What is the order of the star complement ? */
#define STAR_COMPLEMENT_ORDER 20

/* Precision for floating point arithmetic things. */
#define EPS 0.000001

/* graph6 related things */
#define G6LEN(n)  (((n)*((n)-1)/2+5)/6+1) 
#define BIAS6 63
#define TOPBIT6 32


/* 
   At all times this holds the current graph given as a (N+1)x(N+1) adjacency matrix. 
   All the functions isInterlaced, partitionedAM,... always assume this holds the
   current version of our graph - being read from a file or after we have added a new
   vertex to it.
*/
static gsl_matrix *adj;

/* FIXME get rid of this at some point. */
static gsl_matrix *det_adj;

#if N+1 == STAR_COMPLEMENT_ORDER
static gsl_matrix *par_adj; 
static gsl_permutation *par_perm; 
#endif    


static gsl_matrix *E2;

/* Stuff for computing eigenvalues of par_adj */
static gsl_vector *eval_adj;
static gsl_eigen_symm_workspace *w_adj;


/* 
   We need a separate adjacency matrix to test whether our graphs have SC_EIG as an 
   eigenvalue 
*/
static gsl_permutation *p;

/* 
   For every vertex of G we store a bitset indicating its neighbors.
*/
static unsigned nbrs[N+1];

/* 
   Let i,j be two (not necesarily distinct) vertices. Then succWays[i][j][0] is

        BAD_EXTENSION if there exist a valid extension of N(i) \cap N(j) giving us
        a graph with 2 as an eigenvalue 

        OR

        0 <= k <= 2^{N-1} where k is the number of ways to extend N(i) \cap N(j). Note
        if k = 0 then there is no way to extend N(i) \cap N(j) and in particular
        the input subgraph is invalid unless:
            
            - i = j and the degree of i is VAL
            - i is adjacent to j and they have LAMBDA common neighbors.
            - i is not adjacent to j and they have MU common neighbors.
    
    Finally succWays[i][j][1,...,k+1] holds the bit vectors representing the neighbors
    of the newly added vertex.
*/
static unsigned *succWays[N][N];

/* Note. Make sure whatever you use here that its bigger than N */
#define BAD_EXTENSION (1U<<24) 

/* All extended graphs go to this file. */
static FILE *outfile;

/*   
     We often need to determine the number of bits set and it turns out
     this is one of the bottlenecks of this program. Hence we pre-compute
     this. 
*/
static unsigned bitCount[1<<(N+1)];

/*
    The following function accepts a graph6 string for a graph G of order N.

    It fills the adjacency matrix, edges and neighbor information for G.
*/
static void g6toadj(char *s) {

	char *p;
	int i,j,k,x = 0;

    gsl_matrix_set_zero(adj);
    memset(nbrs, 0, sizeof(nbrs));

    for (i = 0; i < N+1; i++) {
        gsl_matrix_set(E2, i, i, (EIG_MAX + (double) (VAL-EIG_MAX)/NTOT));
    }        

    p = s + 1;

    k = 1;
    for (j = 1; j < N; ++j) {
        for (i = 0; i < j; ++i) {
            if (--k == 0) {
		        k = 6;
		        x = *(p++) - BIAS6;
            }
	    
            if (x & TOPBIT6) {

                gsl_matrix_set(adj, i, j, 1);
                gsl_matrix_set(adj, j, i, 1);

                nbrs[i] |= (1U << j);
                nbrs[j] |= (1U << i);

                gsl_matrix_set(E2, i, j, ADJ);
                gsl_matrix_set(E2, j, i, ADJ);
    
            } else {
                gsl_matrix_set(E2, i, j, NONADJ);
                gsl_matrix_set(E2, j, i, NONADJ); 
            }
            x <<= 1;
        }
    }
}

/* 
    This function accepts a string capable of holding the graph6 representation of a graph
    of order N+1 (stored in adj) and returns the respective graph6 string for G 

    NOTE: The code was adapted from B.D.Mckay's nauty implementation.
*/   
static char *adjtog6(char *gcode) {

    unsigned i,j;
    int  k;
    
    char *p,x;

    p = gcode;
    *p++ = BIAS6+N+1;

    k = 6;
    x = 0;
   
    for (j = 1; j < N+1; ++j) {
        for (i = 0; i < j; ++i) {
            x <<= 1;
            if (nbrs[i] & (1U << j)) {
                x |= 1;
            }                

            if (--k == 0) {
                *p++ = BIAS6 + x;
                k = 6;
                x = 0;
            }
        }
    }

    if (k != 6) {
        *p++ = BIAS6 + (x << k);
    }      

    *p++ = '\n';
    *p = '\0';

    return gcode;
}


/*
    This function checks if the graph presented in adj satisfies the conditions for strong-regularity.

    In particular, it returns 
   
    0 if there are two adjacent vertices u and v having more than LAMBDA common neighbors

    1 otherwise.

    TODO. Add condition for non-adjacent vertices as well.

*/  
static unsigned srgCondition(void) {

    unsigned i,j;

    for (i = 0; i < N+1; i++) {
        /* It makes sense to do so because our graphs have small avg. deg */
        if (bitCount[ nbrs[i] ] <= (MU > LAMBDA ? LAMBDA : MU)) {
            continue;
        }

        for (j = i+1; j < N+1; j++) {
            unsigned val = bitCount[ nbrs[i] & nbrs[j] ];
            
            /* ADJACENT */
            if ( (nbrs[i] & (1U << j)) ) {
                if (val > LAMBDA) {
                    return 0;
                }
            } else if (val > MU) {
                return 0;
            }
        }         
    }
    return 1;
}

static unsigned posDefCondition(void) {

    gsl_matrix_memcpy(det_adj, E2);
    gsl_eigen_symm(det_adj, eval_adj, w_adj);

    unsigned i;

    for (i = 0; i < N+1; i++) {
        double eig =  gsl_vector_get (eval_adj, i);

        /* FIXME Argue about the best value here */
        if (eig < -0.00000001) {
            return 0;
        }
    }   
    return 1; 
}

static unsigned weakPosDefCondition(void) {

    gsl_matrix_memcpy(det_adj, E2);

    int signum;

    gsl_linalg_LU_decomp(det_adj, p , &signum);
    double det = gsl_linalg_LU_det(det_adj, signum);
    
    if (det < -0.00001) {
        return 0;
    }

    return 1; 
}

/*
    This function returns 1 if the graph stored in adj does not have SC_EIG as an eigenvalue and 0 otherwise.
*/    
static unsigned validSC(void) {

    unsigned i;

    /* FIXME is det_adj actually changed or can we work on adj???? */
    gsl_matrix_memcpy(det_adj, adj);
    
    for (i = 0; i < N+1; i++) {
        gsl_matrix_set(det_adj, i, i, -SC_EIG);   
    }
    int signum;

    gsl_linalg_LU_decomp(det_adj, p , &signum);
    double det = gsl_linalg_LU_det(det_adj, signum);

    if (fabs(det) > EPS) {
        return 1;
    }
    return 0;
}


#if N+1 == STAR_COMPLEMENT_ORDER
/* 
    This code returns 1 if and only if SC_EIG is an eigenvalue
    of the partitioned matrix of our current graph.
*/   
unsigned partAMcond(void) {

    unsigned i,j;
    unsigned total_edges = 0; 

    gsl_matrix_set_zero(par_adj);

    for (i = 0; i < N+1; i++) {
        gsl_matrix_set(par_adj, i, i, -SC_EIG);
        for (j = i+1; j < N+1; j++) {
            if (nbrs[i] & (1U << j))  {
                gsl_matrix_set(par_adj, i, j, 1);
                gsl_matrix_set(par_adj, j, i, 1);
            }   
        }
        unsigned deg = bitCount[ nbrs[i] ];

        gsl_matrix_set(par_adj, i, N+1, VAL-deg);
        gsl_matrix_set(par_adj, N+1, i, (double)(VAL-deg)/(NTOT-N-1));
        
        total_edges += deg;
    }

    gsl_matrix_set(par_adj, N+1, N+1, -SC_EIG + (double) 2*(NTOT*VAL/2 + total_edges/2 - (N+1)*VAL)/(NTOT-N-1));
    
    int signum;
    gsl_linalg_LU_decomp(par_adj, par_perm, &signum);
    double det = gsl_linalg_LU_det(par_adj, signum);

    /* The determinant is not zero. Hence SC_EIG is not an eigenvalue of our graph. */
    if (fabs(det) > 0.0000001) {
        return 0;
    }
    return 1;
}
#endif


/* 
    expand works as follows. Let G be the graph represented by adj and V the new
    vertex that we wish to introduce.

    For every subset S of V(G) we create the graph G' by joining V to the vertices in S.
    (the subset is called joinVerts in the code) 

    If G' satisfies the strong-regularity condition and is interlacing the eigenvalues of 
    our SRG we have found a new candidate. We store the subset S (as an unsigned int)
    into the variable codes and let graphCodes[S] represent its graph6 string.

    Now for every vertex v S and every pair u,v in S we have just found a valid extension
    of N(v) and N(v) \cap N(u) respectively. If SC_EIG is not an eigenvalue of G' we increase
    the number of good extensions of N(v) and N(v) \cap N(u). We call such a graph G' a good 
    extension.

    If however G' is not a good extension then we mark the fact that not all extensions of N(v)
    and N(v) \cap N(u) give good extensions.

    If in the end we find a vertex v such that all extensions of N(v)  are good or a pair u,v \in V(G)
    such that all extensions N(u) \cap N(v) are good we output the extension of a vertex (resp pair)
    giving us the least number of graphs.

    If this is not the case, we output all extensions stored in graphCodes. 
*/   
static unsigned goodCands = 0;
static unsigned nproc = 0;

char graphCodes[1<<N][G6LEN(N+1)+2];


static void expand(void) {

    unsigned i,j;
    unsigned joinVerts;

    
    /* each valid graph is stored into this string */
    unsigned codes[1<<N];
    unsigned totalGraphs = 0;    

    for (i = 0; i < N ; i++) {
        for (j = i ; j < N; j++) {
            succWays[i][j][0] = 0;
        }
    }

    for (joinVerts = 0; joinVerts < (1U << N); joinVerts++) {

        /* 
            Adding edges 
         
            NOTE: We need not clean E2 since we always set it appropriately. 
        */ 
        for (j = 0; j < N; j++) {
            if (joinVerts & (1U << j) ) {
                gsl_matrix_set(adj, N, j, 1);
                gsl_matrix_set(adj, j, N, 1);

                gsl_matrix_set(E2, N, j, ADJ);
                gsl_matrix_set(E2, j, N, ADJ);

                nbrs[j] |= (1U << N);

            } else {
                gsl_matrix_set(E2, N, j, NONADJ);
                gsl_matrix_set(E2, j, N, NONADJ);
            }
        }
        nbrs[N] = joinVerts;

        if (srgCondition() && /*weakPosDefCondition() &&*/ posDefCondition() 
            #if N+1 == STAR_COMPLEMENT_ORDER
            && partAMcond()                
            #endif                
                ) { 


            unsigned val = validSC();

            /* We only keep graphs not having 2 as an eigenvalue when we are doing our last step */            
            #if N+1 == STAR_COMPLEMENT_ORDER            
            if (val == 1) {
                adjtog6(graphCodes[joinVerts]);
                codes[totalGraphs++] = joinVerts;
            }
            #else
           /* Lets register that we got a good graph */
            adjtog6(graphCodes[joinVerts]);
            codes[totalGraphs++] = joinVerts;
            #endif            

            unsigned k,l;

            for (k = 0; k < N ; k++) {
                /* G' is not an extension for the vertex k */
                if (!(joinVerts & (1U << k))) {
                    continue;
                }
                    
                /* 
                    So we have added an edge to the vertex v, and this gives us a bad extension. 
                */
                if (val == 0) {
                    succWays[k][k][0] = BAD_EXTENSION;
                }
                /* 
                    We have a good extension and if all previous extensions were good as well
                    we can mark this fact by increasing the total number of good extensions found
                    so far 
                */
                else if (succWays[k][k][0] != BAD_EXTENSION) {
                    succWays[k][k][0]++;
                    /* Let us also store the bit string giving the successful code */
                    succWays[k][k][ succWays[k][k][0] ] = joinVerts;
                }
                                
                for (l = k+1; l < N ; l++) {
                    if (joinVerts & (1U << l)) {
                        /* We repeat essentialy the same thing as above yet considering 
                           two pairs of vertices 
                        */
                        if (val == 0) {
                            succWays[k][l][0] = BAD_EXTENSION;
                        } else if (succWays[k][l][0] != BAD_EXTENSION) {
                            succWays[k][l][0]++;
                            
                            succWays[k][l][ succWays[k][l][0] ] = joinVerts;
                        }
                    }               
                }
            }                  
        }
        /* Cleanup */            
        for (j = 0; j < N; j++) {
            if ( joinVerts & (1 << j) ) {
                gsl_matrix_set(adj, N, j, 0);
                gsl_matrix_set(adj, j, N, 0);
                nbrs[j] &= ~(1U << N);
            }
        }
    }
    
    /* 
       This is the last part of expand. We need to check if
       there is a choice of u,v so that all extensions of N(v) 
       or N(v) \cap N(u) are good.
    */
    int curMin = (1<<N);
    int min_i,min_j ;
    min_i = min_j = -1;

    for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {

            /* there is no way to extend i and j */
            if (succWays[i][j][0] == 0) {
                /* 
                   This implies we cannot extend a vertex. 

                   In our cases this means that G is invalid and
                   there is not much we can do.

                   NOTE. In dealing with other SRG's we'd have to add
                   a condition here that the degree of i < VAL !!! 
                */                   
                if (i == j) {
                    return;
                }

                unsigned val1 = bitCount[ nbrs[i] & nbrs[j] ];

                unsigned val2 = nbrs[i] & (1U<< j);
            
                if ( (val2 && val1 < LAMBDA) || (!val2 && val1 < MU)) {
                    return;
                }

                /* We just hit a vertex (pair) that is fully extended. */
                continue;
            }

            if (succWays[i][j][0] != BAD_EXTENSION && succWays[i][j][0] < curMin) {
                curMin = succWays[i][j][0];
                min_i = i;
                min_j = j;
            }
        }
    }

    /* There is no way to only obtain good extensions. Output all graphs */ 
    if (min_i == -1) {
        for (i = 0; i < totalGraphs; i++) {
            fputs(graphCodes[codes[i]], outfile);
        }
    } else {

        goodCands+=1;

        for (i = 1; i <= succWays[min_i][min_j][0]; i++) {
            fputs(graphCodes[ succWays[min_i][min_j][i] ], outfile);
        }
    }
}


void sig_handler(__attribute__((unused)) int signo) {
    fprintf(stderr, "Current progress %u\n", nproc);
}


int main(int argc, char *argv[]) {

    static FILE *infile;
    char line[G6LEN(N)+2];

    assert(argc > 1);

    infile = fopen(argv[1], "r");

    char buf[512];
    snprintf(buf, sizeof(buf), "%s.out", argv[1]);
    outfile = fopen(buf, "w");

    assert(infile && outfile);

    adj = gsl_matrix_alloc(N+1, N+1);
    E2 = gsl_matrix_alloc(N+1, N+1);

    det_adj = gsl_matrix_alloc(N+1, N+1);
    p = gsl_permutation_alloc(N+1);

#if N+1 == STAR_COMPLEMENT_ORDER
    par_adj = gsl_matrix_alloc(N+2, N+2);
    par_perm = gsl_permutation_alloc(N+2);

    assert(par_adj && par_perm);
#endif    

    assert(p && adj && E2 && det_adj);
    
    w_adj = gsl_eigen_symm_alloc(N+1);
    eval_adj = gsl_vector_alloc(N+1);

    assert(w_adj && eval_adj);

    signal(SIGUSR1, sig_handler);
    unsigned i,j;

    for (i = 0; i < N ; i++) {
        for (j = i; j < N; j++) {
            succWays[i][j] = malloc(sizeof(int) * ( (1U<< (N-1))+1));
            assert(succWays[i][j]);
        }
    }
    for (i = 0; i < (1U<<(N+1)); i++) {
        bitCount[i] = __builtin_popcount(i);
    }    

    while (1) {
        if (fgets(line, sizeof(line), infile) == NULL) 
            break;

        g6toadj(line);
        expand();
        nproc++;

	}
    printf("Successfuly extended %u graphs. Out of them %u were good candidates.\n", nproc, goodCands);
    return 0;

}
