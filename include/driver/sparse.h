///
/// \file mini_spasm.h
/// \author Charles Bouillaguet (Github: cbouilla)
/// Modified by Jérôme Bonacchi (Github: MarsParallax)
/// \brief This file contains structures and functions to manage sparse matrix
/// formats. Jérôme Bonacchi changed typedef, types, function names...
/// \date 2020
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_DRIVER_SPARSE_H
#define INCLUDE_DRIVER_SPARSE_H

#include <stdint.h>
#include <stdio.h>

///
/// \brief Type used for dimensions and entries.
///
typedef uint32_t u32;

///
/// \brief Triplet corresponding to a matrix entry.
///
typedef struct
{
  u32 i;    ///< the row index
  u32 j;    ///< the column index
  double x; ///< the numerical value
} mtx_entry;

///
/// \brief Matrix in COO format
///
typedef struct
{
  u32 nnz_max; ///< the maximum number of entries (only nonzero entries are
               /// usually stored, but explicit zero entries are allowed)
  u32 nnz;     ///< the actual number of entries
  u32 n;       ///< the number of rows
  u32 m;       ///< the number of columns
  u32 *i;      ///< the row indices (size == `nnz_max`)
  u32 *j;      ///< the column indices (size == `nnz_max`)
  double *x;   ///< the numerical values (size == `nnz_max`)
} mtx_COO;

///
/// \brief Matrix in compressed-sparse row (CSR) format
///
typedef struct
{
  u32 nnz_max; ///< the maximum number of entries (only nonzero entries are
               /// usually stored, but explicit zero entries are allowed)
  u32 n;       ///< the number of rows
  u32 m;       ///< the number of columns
  u32 *p;      ///< the row pointers (size == `n` + 1)
  u32 *j;      ///< the column indices (size == `nnz_max`)
  double *x;   ///< the numerical values (size == `nnz_max`)
} mtx_CSR;

/* example (this is Matrix/t1)

    [ 4.5  0.0  3.2  0.0 ]
    [ 3.1  2.9  0.0  0.9 ]
A = [ 0.0  1.7  3.0  0.0 ]
    [ 3.5  0.4  0.0  1.0 ]

COO form (nz != -1) :

i = {   2,   1,   3,   0,   1,   3,   3,   1,   0,   2 }
j = {   2,   0,   3,   2,   1,   0,   1,   3,   0,   1 }
x = { 3.0, 3.1, 1.0, 3.2, 2.9, 3.5, 0.4, 0.9, 4.5, 1.7 }

the coefficients may appear in any order.

Compressed Row (CSR) form :

p = {   0,             3,             6,        8,      10 }
i = {   0,   1,   3,   1,   2,   3,   0,   2,   1,   3 }
x = { 4.5, 3.1, 3.5, 2.9, 1.7, 0.4, 3.2, 3.0, 0.9, 1.0 }

In particular, the actual number of nnz is p[n]. Coefficients of a row need not
be sorted by column index.

The numerical values are optional (useful for storing a sparse graph, or the
pattern of a matrix). */

///
/// \brief Returns the number of entries in a matrix in CSR format.
/// Usually, only nonzero entries are stored, hence this is the number of
/// nonzero entries, but explicit zero entries are allowed.
///
/// \param[in] A the matrix in CSR format
/// \return u32 the number of entries in the matrix
///
u32 mtx_nnz(const mtx_CSR *A);

///
/// \brief Gets the time.
///
/// \return double the current time.
///
double spasm_wtime(void);

///
/// \brief Returns a string representing `n` in 4 bytes.
///
/// \param[in] n the integer to print
/// \param[out] target the string representing the integer in human readable
/// format
///
void spasm_human_format(const int64_t n, char *target);

///
/// \brief Allocates a vector.
///
/// \param[in] size the size of the vector to allocate
/// \return void* the pointer returned by `malloc`
///
void *spasm_malloc(const u32 size);

///
/// \brief Allocates a matrix in COO format.
///
/// \param[in] nnz_max the maximum number of nonzero entries
/// \return mtx_COO* the allocated matrix in COO format
///
mtx_COO *mtx_COO_alloc(const u32 nnz_max);

///
/// \brief Frees a matrix in COO format.
///
/// \param[in] A the matrix in COO format
///
void mtx_COO_free(mtx_COO *A);

///
/// \brief Allocates a matrix in CSR format.
///
/// \param[in] n the number of rows
/// \param[in] m the number of columns
/// \param[in] nnz_max the maximum number of nonzero entries
/// \return spasm* the allocated matrix in CSR format
///
mtx_CSR *mtx_CSR_alloc(const u32 n, const u32 m, const u32 nnz_max);

///
/// \brief Frees a matrix in CSR format.
///
/// \param[in] A the matrix in CSR format
///
void mtx_CSR_free(mtx_CSR *A);

///
/// \brief Loads the matrix from a file in Matrix Market format.
/// Heavily inspired by the example program:
/// http://math.nist.gov/MatrixMarket/mmio/c/example_read.c
///
/// \param[in] f the file descriptor of a file in Matrix Market format
/// \return mtx_COO* the allocated matrix in COO format
///
mtx_COO *mtx_load_mm(FILE *f);

///
/// \brief Mutiplies a matrix in CSR format with a vector.
///
/// \param[in] A the input matrix in CSR format
/// \param[in] x the input vector
/// \param[out] y the product vector
///
void mtx_CSR_gemv(const mtx_CSR *A, const double *x, double *y);

///
/// \brief Mutiplies a matrix in COO format with a vector.
///
/// \param[in] T the input matrix in COO format
/// \param[in] x the input vector
/// \param[out] y the product vector
///
void mtx_COO_gemv(const mtx_COO *T, const double *x, double *y);

///
/// \brief Transposes a matrix in COO format.
///
/// \param[in] T the input matrix in COO format
/// \param[out] R the output matrix in COO format
///
void mtx_COO_transpose(const mtx_COO *T, mtx_COO *R);

///
/// \brief Returns the maximum between two integers.
///
/// \return u32 `a` if it is greater than `b`, `b` otherwise
///
static inline u32 spasm_max(const u32 a, const u32 b)
{
  return (a > b) ? a : b;
}

///
/// \brief Returns the minimum between two integers
///
/// \return u32 `a` if it is lesser than `b`, `b` otherwise
///
static inline u32 spasm_min(const u32 a, const u32 b)
{
  return (a < b) ? a : b;
}

// static inline void spasm_swap(int *a, int i, int j) {
// 	int x = a[i];
// 	a[i] = a[j];
// 	a[j] = x;
// }

// static inline int spasm_row_weight(const mtx_CSR * A, int i) {
// 	int *Ap = A->p;
// 	return Ap[i + 1] - Ap[i];
// }

///
/// \brief Asserts that a matrix in COO format is equal to another matrix in
/// CSR format by doing a matrix-vector product
///
/// \param[in] T the matrix in COO format
/// \param[in] A the matrix in CSR format
///
void check(const mtx_COO *T, const mtx_CSR *A);

#endif /* INCLUDE_DRIVER_SPARSE_H */
