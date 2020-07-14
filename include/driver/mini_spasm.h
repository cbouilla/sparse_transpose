#ifndef INCLUDE_MINI_SPASM_H
#define INCLUDE_MINI_SPASM_H

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

#endif /* INCLUDE_MINI_SPASM_H */
