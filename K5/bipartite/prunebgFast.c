/* 
    From nauty's directory execute:

         gcc -o genbg -O2 -DMAXN=32 -DPRUNE1=prune genbg.c gtools.o nauty1.o nautil1.o naugraph1.o schreier.o naurng.o <path to prunebg.c>  -I. -lgsl -lgslcblas -lm 
*/  

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include <assert.h>

#include "nauty.h"

/* Upper bound for the argument to geng */ 
#define N 32

/* ================================= SETTINGS FOR PRUNE =================================  */
#define NTOT 95 
#define VAL 40

#define LAMBDA 12
#define MU  20

#define SC_EIG 2
#define SC_ORDER 20

#define EIG_MAX 2 
#define EIG_MIN -10

#define DIAG (EIG_MAX + (double) (VAL-EIG_MAX)/NTOT)
#define ADJ (-1 + (double) (VAL-EIG_MAX)/NTOT)
#define NONADJ ((double) (VAL-EIG_MAX)/NTOT)

#define EPS 0.000001

const int eigenvalues[NTOT] = {40, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10};

/* ================================= SETTINGS FOR PRUNE =================================  */

static unsigned int firstCall = 1;

static gsl_matrix *adj[N+1];
static gsl_vector_complex  *eval[N+1];
static gsl_eigen_nonsymm_workspace *w[N+1];

static gsl_permutation *p; 
static gsl_eigen_symm_workspace *w_symm[N+1];
static gsl_vector *eval_symm[N+1];


static void init_memory(void) {

    unsigned i;

    p = gsl_permutation_alloc(SC_ORDER);

    assert(p);

    for (i = 1; i <= N; i++) {

        adj[i] = gsl_matrix_alloc(i+1,i+1);
        eval[i] = gsl_vector_complex_alloc(i+1);
        w[i] = gsl_eigen_nonsymm_alloc(i+1);

        eval_symm[i] = gsl_vector_alloc(i);
        w_symm[i] = gsl_eigen_symm_alloc(i);

        assert(w[i] && eval[i] && w_symm[i] && eval_symm[i] && adj[i]); 
    }
}

static void partitioned_am(unsigned *nbrs, unsigned n) {

    unsigned i,j;

    gsl_matrix_set_zero(adj[n]);

    unsigned total_edges = 0;

    for (i = 0; i < n; i++) {
        for (j = i+1; j < n; j++) {
            if (nbrs[i] & (1U << j)) {
                gsl_matrix_set(adj[n], i, j, 1);
                gsl_matrix_set(adj[n], j, i, 1);
            }    
        }
    }

    for (i = 0; i < n; i++) {

        unsigned deg = __builtin_popcount(nbrs[i]);

        gsl_matrix_set(adj[n], i, n, VAL-deg);
        gsl_matrix_set(adj[n], n, i, (double)(VAL-deg)/(NTOT-n));

        total_edges+=deg;
    }

    gsl_matrix_set(adj[n], n, n, (double) 2*(NTOT*VAL/2 + total_edges/2 - n*VAL)/(NTOT-n));
}

static unsigned pruneE2(gsl_matrix *E2, unsigned *nbrs, unsigned n) {

    unsigned i,j;

    gsl_matrix_set_all(E2, NONADJ);

    for (i = 0; i < n; i++) {
        gsl_matrix_set(E2, i,i, DIAG);
        for (j = i+1; j < n; j++) {
            if (nbrs[i] & (1U << j)) {
                gsl_matrix_set(E2, i, j, ADJ);
                gsl_matrix_set(E2, j, i, ADJ);
            }    
        }
    }

    gsl_eigen_symm(E2, eval_symm[n], w_symm[n]);

    for (i = 0; i < n; i++) {
        double eig =  gsl_vector_get(eval_symm[n], i);
        /* FIXME Argue about the best value here */
        if (eig < -0.00000001) {
            return 1;
        }
    }
    return 0; 
}



static int cmp(const void* elem1, const void* elem2) {

    const double *a = elem1, *b = elem2; 
    
    return *a > *b ? -1 : *a < *b ? 1 : 0;
}

static double *spectrum(unsigned n) {

    unsigned i;
    static double eigs[N+1];

    gsl_eigen_nonsymm (adj[n], eval[n], w[n]);    

    /* NOTE. Our matrix is of dimension (n+1) x (n+1). */
    for(i = 0; i < n+1; i++) {
        gsl_complex eval_i = gsl_vector_complex_get (eval[n], i);
        eigs[i] = GSL_REAL(eval_i);
    }

    qsort(eigs, n+1, sizeof(double), cmp); 

    return eigs;
}

/* 
   This function returns 0 if and only if 
   the sequence eigs_sc interlaces eigs 
*/
static unsigned does_interlace(double eigs_sc[N+1], unsigned n) {

    unsigned i;

    for (i = 0; i < n; i++) {

        double expr = eigenvalues[i]-eigs_sc[i];

        if (expr <= EPS && fabs(expr) >= EPS) { 
            return 1;    
        }

        expr = eigs_sc[i] - eigenvalues[NTOT-n+i];

        if (expr <= EPS && fabs(expr) >= EPS) {
            return 1;
        } 
    }
    return 0; 
}

static void makeNbrsGraph(graph *g, unsigned *nbrs, int n1, int n2) {
    
    unsigned i, j;
    set *gi;

    for (i = 0; i < n1; i++) {
        gi = GRAPHROW(g, i, 1);
        /* Fixme. Do we actually have n1+n2+1 here? */
        for (j = n1; j < n1+n2; j++) {
            if (ISELEMENT(gi, j)) {
                nbrs[i] |= (1U << j);
                nbrs[j] |= (1U << i);
            }    
        }
    }
    unsigned n = n1+n2;
    /* vertices n,n+1,n+2,n+3,n+4 will be the vertices of our clique */
    for (i = n ; i < n+5; i++) {
        for (j = i+1 ; j < n+5; j++) {
            nbrs[i] |= (1U << j);
            nbrs[j] |= (1U << i);
        }         
    }

    for (i = 0; i < n1; i++) {
        nbrs[i] |= (1U << n);
        nbrs[n] |= (1U << i);
        
        nbrs[i] |= (1U << (n+1));
        nbrs[n+1] |= (1U << i);
    }

    for (i = n1; i < n; i++) {
        nbrs[i] |= (1U << (n+1));
        nbrs[n+1] |= (1U << i);
        
        nbrs[i] |= (1U << (n+2));
        nbrs[n+2] |= (1U << i);
    }
}

unsigned is_validSC(gsl_matrix *adj, unsigned *nbrs, unsigned n) {
    
    unsigned i,j;
        
    gsl_matrix_set_zero(adj);
    
    for (i = 0; i < n; i++) {
        gsl_matrix_set(adj, i, i, -SC_EIG);
        for (j = i+1; j < n; j++) {
            if (nbrs[i] & (1U << j)) {
                gsl_matrix_set(adj, i, j, 1);
                gsl_matrix_set(adj, j, i, 1);
            }
        }
    }

    int signum;

    gsl_linalg_LU_decomp(adj, p , &signum);
    
    double det = gsl_linalg_LU_det(adj, signum);

    if (fabs(det) > EPS) {
        return 1;
    }

    /* The determinant is in fact 0 and hence SC_EIG is an eigenvalue of our graph 
     * and we need to continue 
    */
    return 0;
}

int prune (graph *g, int *deg, int n1, int n2, int maxn2) {

    if (firstCall) {
        init_memory();
        firstCall = 0;
    }

    unsigned nbrs[N];
    memset(nbrs, 0, sizeof(nbrs));
    makeNbrsGraph(g, nbrs, n1, n2);

    unsigned n = n1+n2;

    if (pruneE2(adj[n+5-1], nbrs, n+5)) { 
        return 1;
    }

    /* Comment me. */
    partitioned_am(nbrs, n+5);
    
    double *eigs = spectrum(n+5);

    if (does_interlace(eigs, n+6)) {
        return 1;
    }

    if (n+5 == SC_ORDER) {
        if (is_validSC(adj[SC_ORDER-1], nbrs, SC_ORDER)) {
            return 1;
        } 
    }
    return 0;
}    
