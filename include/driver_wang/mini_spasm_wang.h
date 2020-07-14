#ifndef INCLUDE_MINI_SPASM_H_W
#define INCLUDE_MINI_SPASM_H_W

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>

// typedef uint64_t u64;
// typedef uint64_t u32;

/* --- primary SpaSM routines and data structures --- */


struct matrix_entry_t {
    int i;
    int j;
    double x;
};

typedef struct {                /* matrix in compressed-sparse row format */
	int nzmax;                    /* maximum number of entries */
	int n;                        /* number of rows */
	int m;                        /* number of columns */
	int *p;                       /* row pointers (size n+1) */
	int *j;                       /* column indices, size nzmax */
	double *x;                 /* numerical values, size nzmax (optional) */
} spasm;

typedef struct {                /* matrix in triplet form */
	int nzmax;                    /* maximum number of entries */
	int nz;                       /* # entries */
	int n;                        /* number of rows */
	int m;                        /* number of columns */
	int *i;                       /* row indices, size nzmax */
	int *j;                       /* column indices (size nzmax) */
	double *x;                 /* numerical values, size nzmax (optional) */
} spasm_triplet;


/* example (this is Matrix/t1)

		[ 4.5  0.0  3.2  0.0 ]
		[ 3.1  2.9  0.0  0.9 ]
A = [ 0.0  1.7  3.0  0.0 ]
		[ 3.5  0.4  0.0  1.0 ]

Triplet form (nz != -1) :

i = {   2,   1,   3,   0,   1,   3,   3,   1,   0,   2 }
j = {   2,   0,   3,   2,   1,   0,   1,   3,   0,   1 }
x = { 3.0, 3.1, 1.0, 3.2, 2.9, 3.5, 0.4, 0.9, 4.5, 1.7 }

the coefficients may appear in any order.

Compressed Row form :

p = {   0,             3,             6,        8,      10 }
i = {   0,   1,   3,   1,   2,   3,   0,   2,   1,   3 }
x = { 4.5, 3.1, 3.5, 2.9, 1.7, 0.4, 3.2, 3.0, 0.9, 1.0 }

In particular, the actual number of nnz is p[n]. Coefficients of a row need not be sorted by column index.

The numerical values are optional (useful for storing a sparse graph, or the pattern of a matrix). */

int spasm_nnz(const spasm * A);
void *spasm_malloc(size_t size);
void spasm_triplet_free(spasm_triplet * A);
spasm *spasm_csr_alloc(int n, int m, int nzmax);
void spasm_csr_free(spasm * A);
spasm_triplet *spasm_load_mm(FILE * f);
double spasm_wtime();
void spasm_csr_gemv(const spasm *A, const double *x, double *y);
void spasm_triplet_gemv(const spasm_triplet *T, const double *x, double *y);
void spasm_triplet_transpose(const spasm_triplet *T, spasm_triplet *R);

/* utilities */
static inline int spasm_max(int a, int b) {
	return (a > b) ? a : b;
}

static inline int spasm_min(int a, int b) {
	return (a < b) ? a : b;
}

// static inline void spasm_swap(int *a, int i, int j) {
// 	int x = a[i];
// 	a[i] = a[j];
// 	a[j] = x;
// }

// static inline int spasm_row_weight(const spasm * A, int i) {
// 	int *Ap = A->p;
// 	return Ap[i + 1] - Ap[i];
// }

// #endif /* INCLUDE_MINI_SPASM_H */
#include <assert.h>                                                              
#include <math.h>                                                                
#include <err.h>                                                                 
#include <sys/time.h>                                                            
                                                                                 
                                                         
#include "mmio_wang.h"                                                                
                                                                                 
double spasm_wtime() {                                                           
  struct timeval ts;                                                             
  gettimeofday(&ts, NULL);                                                       
  return (double)ts.tv_sec + ts.tv_usec / 1E6;                                   
}                                                                                
                                                                                 
                                                                                 
int spasm_nnz(const spasm * A) {                                                 
  return A->p[A->n];                                                             
}                                                                                
                                                                                 
                                                                                 
void *spasm_malloc(size_t size) {                                                
  void *x = aligned_alloc(64, size);                                             
  if (x == NULL)                                                                 
    err(1, "malloc failed");                                                     
  return x;                                                                      
}                                                                                
                                                                                 
/* return a string representing n in 4 bytes */                                  
static void spasm_human_format(int64_t n, char *target) {                        
  if (n < 1000) {                                                                
    sprintf(target, "%d", (int) n);                                              
    return;                                                                      
  }                                                                              
  if (n < 1000000) {                                                             
    sprintf(target, "%.1fk", n / 1e3);                                           
    return;                                                                      
  }                                                                              
  if (n < 1000000000) {                                                          
    sprintf(target, "%.1fM", n / 1e6);                                           
    return;                                                                      
  }                                                                              
  if (n < 1000000000000ll) {                                                     
    sprintf(target, "%.1fG", n / 1e9);                                           
    return;                                                                      
  }                                                                              
  if (n < 1000000000000000ll) {                                                  
    sprintf(target, "%.1fT", n / 1e12);                                          
    return;                                                                      
  }                                                                              
}                                                                                
                                                                                 
/* free a sparse matrix */                                                       
void spasm_csr_free(spasm * A) { 
  if (A == NULL)                                                                 
    return;                                                                      
  free(A->p);                                                                    
  free(A->j);                                                                    
  free(A->x);   /* trick : free does nothing on NULL pointer */                  
  free(A);                                                                       
}                                                                                
                                                                                 
void spasm_triplet_free(spasm_triplet * T)                                       
{                                                                                
  free(T->i);                                                                    
  free(T->j);                                                                    
  free(T->x);   /* trick : free does nothing on NULL pointer */                  
  free(T);                                                                       
}                                                                                
                                                                                 
/* allocate a sparse matrix (triplet form) */                                    
spasm_triplet *spasm_triplet_alloc(int nzmax)                                    
{                                                                                
  spasm_triplet *T;                                                              
                                                                                 
  T = (spasm_triplet*)spasm_malloc(sizeof(spasm_triplet));                                       
  T->m = 0;                                                                      
  T->n = 0;                                                                      
  T->nzmax = nzmax;                                                              
  T->nz = 0;                                                                     
  T->i = (int*)spasm_malloc(nzmax * sizeof(int));                                      
  T->j = (int*)spasm_malloc(nzmax * sizeof(int));                                      
  T->x = (double*)spasm_malloc(nzmax * sizeof(double));                                   
  return T;                                                                      
}                                                                                
                                                                                 
/* allocate a sparse matrix (compressed-row form) */                             
spasm *spasm_csr_alloc(int n, int m, int nzmax)                                  
{                                                                                
  spasm *A;                                                                      
                                                                                 
  A = (spasm *)spasm_malloc(sizeof(spasm));  /* allocate the cs struct */                 
  A->m = m;   /* define dimensions and nzmax */                                  
  A->n = n;                                                                      
  A->nzmax = nzmax;                                                              
  A->p = (int*)spasm_malloc((n + 1) * sizeof(A->p));                                   
  A->j = (int*)spasm_malloc(nzmax * sizeof(A->j));                                     
  A->x = (double*)spasm_malloc(nzmax * sizeof(A->x));                                     
  return A;                                                                      
}                                                                                
                                                                                 

/* add an entry to a triplet matrix; enlarge it if necessary */                  
static inline void spasm_add_entry(spasm_triplet * T, int i, int j, double x)    
{                                                                                
  assert((i >= 0) && (j >= 0));                                                  
  int nz = T->nz;                                                                
  assert(nz < T->nzmax);                                                         
                                                                                 
  T->i[nz] = i;                                                                  
  T->j[nz] = j;                                                                  
  T->x[nz] = x;                                                                  
                                                                                 
  T->nz = nz + 1;                                                                
  T->n = spasm_max(T->n, i + 1);                                                 
  T->m = spasm_max(T->m, j + 1);                                                 
}                                                                                
                                                                                 
/*                                                                               
 * Load a matrix in MatrixMarket sparse format.                                  
 * Heavily inspired by the example program:                                      
 *     http://math.nist.gov/MatrixMarket/mmio/c/example_read.c                   
 */                                                                              
spasm_triplet *spasm_load_mm(FILE * f)                                           
{                                                                                
  MM_typecode matcode;                                                           
  int n, m, nnz;                                                                 
                                                                                 
  double start = spasm_wtime();                                                  
  if (mm_read_banner(f, &matcode) != 0)                                          
    errx(1, "Could not process Matrix Market banner.\n");                        
  char *typecode = mm_typecode_to_str(matcode);                                  
                                                                                 
  if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode))                          
    errx(1, "Matrix Market type: [%s] not supported", typecode);                 
                                                                                 
  if (!mm_is_general(matcode))                                                   
    errx(1, "Matrix market type [%s] not supported",  typecode);                 
                                                                                 
  if (mm_is_symmetric(matcode) || mm_is_skew(matcode))                           
    errx(1, "Matrix market type [%s] not supported",  typecode);                 
                                                                                 
  int real = mm_is_real(matcode) || mm_is_integer(matcode);                      
  int pattern = mm_is_pattern(matcode);                                          
                                                                                 
  if (!real && !pattern)                                                         
    errx(1, "Matrix market type [%s] not supported",  typecode);                 
                                                                                 
  if (mm_read_mtx_crd_size(f, &n, &m, &nnz) != 0)                                
    errx(1, "Cannot read matrix size");                                          
                                                                                 
  char s_nnz[16];                                                                
  spasm_human_format(nnz, s_nnz);                                                
  fprintf(stderr, "[IO] loading %d x %d MTX [%s] %s NNZ ...", n, m, typecode, s_nnz);
  fflush(stderr);                                                                
  free(typecode); // mm_typecode_to_str return a malloc'd char*                  
                                                                                 
  spasm_triplet *T = spasm_triplet_alloc(nnz);                                   
                                                                                 
  for (int i = 0; i < nnz; i++) {                                                
    int u, v;                                                                    
    double x;                                                                    
    if (real) {                                                                  
      if (3 != fscanf(f, "%d %d %lg\n", &u, &v, &x))                             
        errx(1, "parse error entry %d\n", i);                                    
      spasm_add_entry(T, u - 1, v - 1, x);                                       
    } else {                                                                     
      if (2 != fscanf(f, "%d %d\n", &u, &v))                                     
        errx(1, "parse error entry %d\n", i);                                    
      spasm_add_entry(T, u - 1, v - 1, 1.0);                                     
    }                                                                            
  }                                                                              
  fprintf(stderr, " Actually %d x %d [%.1fs]\n", T->n, T->m, spasm_wtime() - start);
  return T;                                                                      
}
#endif
