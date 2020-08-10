///
/// \file mini_spasm.h
/// \author // TODO
/// \brief This file implements structures and functions to manage sparse matrix
/// formats.
/// \date 2020-07-23
///
/// @copyright // TODO
///

#ifndef INCLUDE_DRIVER_MINI_SPASM_H
#define INCLUDE_DRIVER_MINI_SPASM_H

#include <stdio.h>

#include <stdint.h>

///
/// \brief Type used for dimensions and entries.
///
typedef uint32_t u32;

///
/// \brief A triplet corresponding to a matrix entry.
///
struct matrix_entry_t
{
  u32 i;    ///< the row index
  u32 j;    ///< the column index
  double x; ///< the numerical value
};

///
/// \brief Matrix in triplet (COO) format
///
typedef struct
{
  u32 nnz_max; ///< the maximum number of entries (only nonzero entries are
               ///< usually stored, but explicit zero entries are allowed)
  u32 nnz;     ///< the actual number of entries
  u32 n;       ///< the number of rows
  u32 m;       ///< the number of columns
  u32 *i;      ///< the row indices (size == `nnz_max`)
  u32 *j;      ///< the column indices (size == `nnz_max`)
  double *x;   ///< the numerical values (size == `nnz_max`)
} spasm_triplet;

///
/// \brief Matrix in compressed-sparse row (CSR) format
///
typedef struct
{
  u32 nnz_max; ///< the maximum number of entries (only nonzero entries are
               ///< usually stored, but explicit zero entries are allowed)
  u32 n;       ///< the number of rows
  u32 m;       ///< the number of columns
  u32 *p;      ///< the row pointers (size == `n` + 1)
  u32 *j;      ///< the column indices (size == `nnz_max`)
  double *x;   ///< the numerical values (size == `nnz_max`)
} spasm;

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
u32 spasm_nnz(const spasm *A);

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
/// \brief Allocates a matrix in triplet format.
///
/// \param[in] nnz_max the maximum number of nonzero entries
/// \return spasm_triplet* the allocated matrix in triplet format
///
spasm_triplet *spasm_triplet_alloc(const u32 nnz_max);

///
/// \brief Frees a matrix in triplet format.
///
/// \param[in] A the matrix in triplet format
///
void spasm_triplet_free(spasm_triplet *A);

///
/// \brief Allocates a matrix in CSR format.
///
/// \param[in] n the number of rows
/// \param[in] m the number of columns
/// \param[in] nnz_max the maximum number of nonzero entries
/// \return spasm* the allocated matrix in CSR format
///
spasm *spasm_csr_alloc(const u32 n, const u32 m, const u32 nnz_max);

///
/// \brief Frees a matrix in CSR format.
///
/// \param[in] A the matrix in CSR format
///
void spasm_csr_free(spasm *A);

///
/// \brief Loads the matrix from a file in Matrix Market format.
/// Heavily inspired by the example program:
/// http://math.nist.gov/MatrixMarket/mmio/c/example_read.c
///
/// \param[in] f the file descriptor of a file in Matrix Market format
/// \return spasm_triplet* the allocated matrix in triplet format
///
spasm_triplet *spasm_load_mm(FILE *f);

///
/// \brief Mutiplies a matrix in CSR format with a vector. //TODO assert dim
///
/// \param[in] A the input matrix in CSR format
/// \param[in] x the input vector
/// \param[out] y the product vector
///
void spasm_csr_gemv(const spasm *A, const double *x, double *y);

///
/// \brief Mutiplies a matrix in triplet format with a vector.//TODO assert dim
///
/// \param[in] T the input matrix in triplet format
/// \param[in] x the input vector
/// \param[out] y the product vector
///
void spasm_triplet_gemv(const spasm_triplet *T, const double *x, double *y);

///
/// \brief Transposes a matrix in triplet format.
///
/// \param[in] T the input matrix in triplet format
/// \param[out] R the output matrix in triplet format
///
void spasm_triplet_transpose(const spasm_triplet *T, spasm_triplet *R);

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

// static inline int spasm_row_weight(const spasm * A, int i) {
// 	int *Ap = A->p;
// 	return Ap[i + 1] - Ap[i];
// }

///
/// \brief Asserts that a matrix in triplet format is equal to another matrix in
/// CSR format by doing a matrix-vector product
///
/// \param[in] T the matrix in triplet format
/// \param[in] A the matrix in CSR format
///
void check(const spasm_triplet *T, const spasm *A);

#endif /* INCLUDE_DRIVER_MINI_SPASM_H */
