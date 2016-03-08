/* 
   Compile as

   gcc compGraph2graph6.c -O2 -lgsl -lgslcblas -o compGraph2graph6
*/   
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

/* Order of input graphs - note this code will only work for orders smaller than 62 */
#define N 20

/* Maximum order we expect for the the comparability graph */
#define MAX_ORDER 75000

/* Forbidden eigenvalue of the star complement */
#define SC_EIG 2 

/* Specified precision used while doing floating point operations we treat EPS as zero */
#define EPS 0.000001

#define EXPECTED_CLIQUE 75

/* graph6 related things */

#define SIZELEN(n) ((n)<=SMALLN?1:((n)<=SMALLISHN?4:8))
#define G6LEN(n) (SIZELEN(n) \
         + ((size_t)(n)/12)*((size_t)(n)-1) + (((size_t)(n)%12)*((size_t)(n)-1)+11)/12)

#define SMALLN 62
#define BIAS6 63
#define TOPBIT6 32
#define SMALLISHN 258047
#define MAXBYTE 126
#define C6MASK 63

static gsl_matrix *mat, *mat_inv;
static gsl_permutation *perm;

static gsl_vector *vecs[1<<N];
static gsl_vector *vecs_prod[1<<N];

static gsl_vector *vec_j;

/* storing the current graph6 string of our graph */
static char line[G6LEN(N)+2];

static char gcode[G6LEN(MAX_ORDER)+3]; 

/* number of currently processed graphs. */ 
static unsigned nproc = 0;
static unsigned skipped = 0;

static FILE *outFile;

static void stringtomat(char *s) {

	char *p;
	int i,j,k,x = 0;

    /* Clear the adjacency matrix */
    gsl_matrix_set_zero(mat);

    p = s + 1;

    k = 1;
    
    gsl_matrix_set(mat, 0, 0, SC_EIG);

    for (j = 1; j < N; ++j) {
        gsl_matrix_set(mat, j, j, SC_EIG);
        for (i = 0; i < j; ++i) {
            if (--k == 0) {
        		k = 6;
		        x = *(p++) - BIAS6;
            }
	    
            if (x & TOPBIT6) {
                gsl_matrix_set(mat, i, j, -1);
                gsl_matrix_set(mat, j, i, -1);
            }
            x <<= 1;
        }
    }

    int signum;
    gsl_linalg_LU_decomp(mat, perm , &signum);
    gsl_linalg_LU_invert(mat, perm, mat_inv);

}

/* Some graph6 string thingie */
static void encodegraphsize(const int n, char **pp) {
    char *p;

    p = *pp;
    if (n <= SMALLN) 
        *p++ = BIAS6 + n;
    else {
        *p++ = MAXBYTE;
        *p++ = BIAS6 + (n >> 12);
        *p++ = BIAS6 + ((n >> 6) & C6MASK);
        *p++ = BIAS6 + (n & C6MASK);
    }
    *pp = p;
}

/* Given that mat_inv is the M = SC_EIG*I - A we 
   compute all binary vectors x such that x M x^t == SC_EIG
   and x M j^t == -1. 

   Finally for any pair x,y of such vectors we add an edge to our 
   graph if x M y^t is either 0 or -1.
*/

/* This is a global thingie since declaring it localy as 
 * an array of size 1<<N breaks the stack limit */
gsl_vector *verts[1<<N];

static void constructGraph(void) {

    double res;
    unsigned i, j;

    /* After the first iteration this holds the number of 
       vertices of the obtained graph. The respective vertices
       are stored in vecs_prod[0]...vecs_prod[cache_size-1].
    */ 
    unsigned cache_size = 0; 

    for (i = 0; i < 1<<N; i++) {
        gsl_blas_dsymv(CblasUpper, 1, mat_inv, vecs[i], 0, vecs_prod[cache_size]);

        gsl_blas_ddot(vecs_prod[cache_size], vec_j, &res);
    
        if (fabs(res+1) < EPS) {
            gsl_blas_ddot(vecs_prod[cache_size], vecs[i], &res);

            if (fabs(res-SC_EIG) < EPS) {
                verts[cache_size] = vecs[i];
                cache_size+=1;
            }                
        }
    }
    
    if (cache_size < EXPECTED_CLIQUE) {
        skipped++;
        return;
    }

    assert(cache_size < MAX_ORDER);

    char *p = gcode;
    encodegraphsize(cache_size, &p);

    int k = 6, x = 0;

    for (i = 1; i < cache_size; i++) {
        for (j = 0; j < i; j++) {
            x <<= 1;
            gsl_blas_ddot(vecs_prod[i], verts[j], &res);

            /* We have an edge */
            if (fabs(res) <= EPS || fabs(res+1) <= EPS) {
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

    fputs(gcode, outFile);
}


static void init_vectors(void) {
    
    unsigned i,j;

    vec_j = gsl_vector_alloc(N);

    assert(vec_j);

    gsl_vector_set_all(vec_j, 1);

    for (i = 0; i < 1<<N; i++) {
        vecs[i] = gsl_vector_calloc(N);
        vecs_prod[i] = gsl_vector_alloc(N);

        assert(vecs[i] && vecs_prod[i]);            

        /* We fill the i'th vector of vecs */
        for (j = 0; j < N; j++) {
            if ( i & (1<<j) ) {
                gsl_vector_set(vecs[i], j, 1);
            }
        }
    }   
}


int main(int argc, char **argv) {
    
    static FILE *infile;

    assert (argc > 1);
    
    infile = fopen(argv[1], "r");
    mat = gsl_matrix_alloc(N, N);
    mat_inv = gsl_matrix_alloc(N, N);
    perm = gsl_permutation_calloc(N);

    char buf[512];

    snprintf(buf, sizeof(buf), "%s.out", argv[1]);
    outFile = fopen(buf, "w");

    assert (outFile && infile && mat && mat_inv && perm);
    
    init_vectors();

    while (1) {
        if (!fgets(line, sizeof(line), infile)) 
            break;
        stringtomat(line);
        constructGraph();
        nproc++;
	}
    
    printf("Successfuly processed: %u graphs. Skipped: %u\n" , nproc, skipped);

    return 0;
}
