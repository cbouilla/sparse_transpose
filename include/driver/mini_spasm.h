#ifndef INCLUDE_DRIVER_MINI_SPASM_H
#define INCLUDE_DRIVER_MINI_SPASM_H

#include <stdint.h>
#include <stdio.h>

typedef uint32_t u32;

///
/// \brief A triplet corresponding to a matrix entry.
///
struct matrix_entry_t
{
  u32 i;    /// row index
  u32 j;    /// column index
  double x; /// numerical value
};

///
/// \brief Matrix in triplet format
///
typedef struct
{
  u32 nnz_max; /// maximum number of entries
  u32 nnz;     /// # entries
  u32 n;       /// number of rows
  u32 m;       /// number of columns
  u32 *i;      /// row indices, size nzmax
  u32 *j;      /// column indices (size nzmax)
  double *x;   /// numerical values, size nzmax (optional)
} spasm_triplet;

///
/// \brief Matrix in compressed-sparse row (CSR) format
///
typedef struct
{
  u32 nnz_max; /// maximum number of entries
  u32 n;       /// number of rows
  u32 m;       /// number of columns
  u32 *p;      /// row pointers (size n+1)
  u32 *j;      /// column indices, size nzmax
  double *x;   /// numerical values, size nzmax (optional)
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
/// \brief Returns the number of entries in a matrix in CSR format. This is
/// usually the number of nonzero entries, but explicit zero entries are
/// allowed.
///
/// \param A The matrix in CSR format
/// \return int The number of entries in `A`
///
int spasm_nnz(const spasm *A);

///
/// \brief Allocates a vector of size `size`.
///
/// \param size The size to allocate
/// \return void* The pointer returned by `malloc`
///
void *spasm_malloc(size_t size);

///
/// \brief Frees a matrix in triplet format
///
/// \param A The matrix in triplet format
///
void spasm_triplet_free(spasm_triplet *A);

///
/// \brief Allocates a matrix in CSR format
///
/// \param n The number of rows
/// \param m The number of columns
/// \param nzmax The maximum number of entries
/// \return spasm* The allocated matrix in CSR format
///
spasm *spasm_csr_alloc(int n, int m, int nzmax);

///
/// \brief Frees a matrix in CSR format
///
/// \param A The matrix in CSR format
///
void spasm_csr_free(spasm *A);

///
/// \brief Loads the matrix from a file in Matrix Market format
///
/// \param f The file in Matrix Market format
/// \return spasm_triplet* The allocated matrix in triplet format
///
spasm_triplet *spasm_load_mm(FILE *f);

///
/// \brief Gets the time.
///
/// \return double The current time.
///
double spasm_wtime();

///
/// \brief Mutiplies a matrix in CSR format with a vector
///
/// \param A The matrix in CSR format
/// \param x The input vector
/// \param y The resulting vector
///
void spasm_csr_gemv(const spasm *A, const double *x, double *y);

///
/// \brief Mutiplies a matrix in triplet format with a vector
///
/// \param T The matrix in triplet format
/// \param x The input vector
/// \param y The resulting vector
///
void spasm_triplet_gemv(const spasm_triplet *T, const double *x, double *y);

///
/// \brief Transposes a matrix in triplet format
///
/// \param T The input matrix in triplet format
/// \param R The output matrix in triplet format
///
void spasm_triplet_transpose(const spasm_triplet *T, spasm_triplet *R);

///
/// \brief Returns the maximum between two integers
///
/// \return int `a` if it is greater than `b`, `b` otherwise
///
static inline int spasm_max(int a, int b) { return (a > b) ? a : b; }

///
/// \brief Returns the minimum between two integers
///
/// \return int `a` if it is lesser than `b`, `b` otherwise
///
static inline int spasm_min(int a, int b) { return (a < b) ? a : b; }

// static inline void spasm_swap(int *a, int i, int j) {
// 	int x = a[i];
// 	a[i] = a[j];
// 	a[j] = x;
// }

// static inline int spasm_row_weight(const spasm * A, int i) {
// 	int *Ap = A->p;
// 	return Ap[i + 1] - Ap[i];
// }

#endif /* INCLUDE_DRIVER_MINI_SPASM_H */
