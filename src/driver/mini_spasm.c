///
/// \file mini_spasm.c
/// \author // TODO
/// \brief This file implements structures and functions to manage sparse matrix
/// formats.
/// \date 2020-07-23
///
/// @copyright // TODO
///

#include <assert.h>
#include <err.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <sys/time.h>

#include "mini_spasm.h"
#include "mmio.h"

double spasm_wtime()
{
  struct timeval ts;
  gettimeofday(&ts, NULL);
  return (double)ts.tv_sec + ts.tv_usec / 1E6;
}

u32 spasm_nnz(const spasm *A) { return A->p[A->n]; }

void *spasm_malloc(const u32 size)
{
  void *x = aligned_alloc(64, size);
  if (x == NULL)
    err(1, "malloc failed");
  return x;
}

void spasm_human_format(int64_t n, char *target)
{
  if (n < 1000) // 10^3
  {
    sprintf(target, "%d", (int)n);
    return;
  }
  if (n < 1000000) // 10^6
  {
    sprintf(target, "%.1fk", n / 1e3);
    return;
  }
  if (n < 1000000000) // 10^9
  {
    sprintf(target, "%.1fM", n / 1e6);
    return;
  }
  if (n < 1000000000000ll) // 10^12
  {
    sprintf(target, "%.1fG", n / 1e9);
    return;
  }
  if (n < 1000000000000000ll) // 10^15
  {
    sprintf(target, "%.1fT", n / 1e12);
    return;
  }
}

void spasm_csr_free(spasm *A)
{
  if (A == NULL)
    return;
  free(A->p);
  free(A->j);
  free(A->x); /* trick : free does nothing on NULL pointer */
  free(A);
}

void spasm_triplet_free(spasm_triplet *T)
{
  free(T->i);
  free(T->j);
  free(T->x); /* trick : free does nothing on NULL pointer */
  free(T);
}

spasm_triplet *spasm_triplet_alloc(const u32 nnz_max)
{
  spasm_triplet *T = spasm_malloc(sizeof(spasm_triplet));
  T->nnz_max = nnz_max;
  T->nnz = 0;
  T->n = 0;
  T->m = 0;
  T->i = spasm_malloc(nnz_max * sizeof(int));
  T->j = spasm_malloc(nnz_max * sizeof(int));
  T->x = spasm_malloc(nnz_max * sizeof(double));
  return T;
}

spasm *spasm_csr_alloc(const u32 n, const u32 m, const u32 nnz_max)
{
  spasm *A = spasm_malloc(sizeof(spasm));
  A->nnz_max = nnz_max;
  A->n = n;
  A->m = m;
  A->p = spasm_malloc((n + 1) * sizeof(A->p));
  A->j = spasm_malloc(nnz_max * sizeof(A->j));
  A->x = spasm_malloc(nnz_max * sizeof(A->x));
  return A;
}

///
/// \brief Adds an entry to a matrix in triplet format.
/// The dimensions of the matrix must be greater than both i and j. Dimensions
/// are not enlarge.
///
/// \param[in, out] T the input matrix in triplet format
/// \param[in] i the row index
/// \param[in] j the column index
/// \param[in] x the numerical value
///
inline void spasm_add_entry(spasm_triplet *T, const u32 i, const u32 j,
                            const double x)
{
  // Checks the matrix dimensions
  // assert((i >= 0) && (i < T->n) && (j >= 0) && (j < T->m));
  assert((i < T->n) && (j < T->m));
  const u32 nnz = T->nnz;
  assert(nnz < T->nnz_max);

  T->i[nnz] = i;
  T->j[nnz] = j;
  T->x[nnz] = x;
  // T->nnz++
  T->nnz = nnz + 1;

  // Enlarges the matrix dimensions if necessary
  // T->n = spasm_max(T->n, i + 1);
  // T->m = spasm_max(T->m, j + 1);
}

spasm_triplet *spasm_load_mm(FILE *f)
{
  MM_typecode matcode;
  u32 n, m, nnz;

  double start = spasm_wtime();
  if (mm_read_banner(f, &matcode) != 0)
    err(1, "Could not process Matrix Market banner.\n");
  char *typecode = mm_typecode_to_str(matcode);

  if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode))
    err(1, "Matrix Market type: [%s] not supported", typecode);

  if (!mm_is_general(matcode))
    err(1, "Matrix market type [%s] not supported", typecode);

  if (mm_is_symmetric(matcode) || mm_is_skew(matcode))
    err(1, "Matrix market type [%s] not supported", typecode);

  bool real = mm_is_real(matcode) || mm_is_integer(matcode);
  bool pattern = mm_is_pattern(matcode);

  if (!real && !pattern)
    err(1, "Matrix market type [%s] not supported", typecode);

  if (mm_read_mtx_crd_size(f, &n, &m, &nnz) != 0)
    err(1, "Cannot read matrix size");

  char s_nnz[16];
  spasm_human_format(nnz, s_nnz);
  fprintf(stderr, "[IO] loading %d x %d MTX [%s] %s NNZ ...", n, m, typecode,
          s_nnz);
  fflush(stderr);
  free(typecode); // mm_typecode_to_str return a malloc'd char*

  spasm_triplet *T = spasm_triplet_alloc(nnz);
  T->n = n;
  T->m = m;
  for (u32 i = 0; i < nnz; i++)
  {
    u32 u, v;
    double x;
    if (real)
    {
      if (3 != fscanf(f, "%d %d %lg\n", &u, &v, &x))
        err(1, "parse error entry %d\n", i);
      spasm_add_entry(T, u - 1, v - 1, x);
    }
    else
    {
      if (2 != fscanf(f, "%d %d\n", &u, &v))
        err(1, "parse error entry %d\n", i);
      spasm_add_entry(T, u - 1, v - 1, 1.0);
    }
  }
  fprintf(stderr, " Actually %d x %d [%.1fs]\n", T->n, T->m,
          spasm_wtime() - start);
  return T;
}

void spasm_triplet_gemv(const spasm_triplet *T, const double *x, double *y)
{
  const u32 *Ti = T->i;
  const u32 *Tj = T->j;
  const double *Tx = T->x;
  u32 m = T->m;
  u32 nnz = T->nnz;

  for (u32 j = 0; j < m; j++)
    y[j] = 0;
  for (u32 k = 0; k < nnz; k++)
  {
    u32 i = Ti[k];
    u32 j = Tj[k];
    y[j] += Tx[k] * x[i];
  }
}

void spasm_csr_gemv(const spasm *A, const double *x, double *y)
{
  const u32 *Ap = A->p;
  const u32 *Aj = A->j;
  const double *Ax = A->x;
  u32 n = A->n;
  u32 m = A->m;

  for (u32 j = 0; j < m; j++)
    y[j] = 0;
  for (u32 i = 0; i < n; i++)
    for (u32 k = Ap[i]; k < Ap[i + 1]; k++)
    {
      u32 j = Aj[k];
      y[j] += Ax[k] * x[i];
    }
}

void spasm_triplet_transpose(const spasm_triplet *T, spasm_triplet *R)
{
  R->i = T->j;
  R->j = T->i;
  R->x = T->x;
  R->n = T->m;
  R->m = T->n;
  R->nnz = T->nnz;
  R->nnz_max = T->nnz_max;
}

void check(const spasm_triplet *T, const spasm *A)
{
  const u32 n = T->n;
  const u32 m = T->m;
  const u32 nnz = T->nnz;

  // safety checks
  assert(A->n == n);
  assert(A->m == m);
  const u32 *Ap = A->p;
  const u32 *Aj = A->j;
  assert(Ap[0] == 0);
  for (u32 i = 0; i < n; i++)
    assert(Ap[i] <= Ap[i + 1]);
  assert(Ap[n] == nnz);
  for (u32 k = 0; k < nnz; k++)
  {
    // assert(0 <= Aj[k]); // useless with u32
    assert(Aj[k] < m);
  }

  // matrix-vector product
  double *X = (double *)malloc(n * sizeof(double));
  double *Ya = (double *)malloc(m * sizeof(double));
  double *Yb = (double *)malloc(m * sizeof(double));
  for (u32 i = 0; i < n; i++)
    X[i] = drand48();
  spasm_triplet_gemv(T, X, Ya);
  spasm_csr_gemv(A, X, Yb);

  double error = 0;
  for (u32 j = 0; j < m; j++)
  {
    double x = Ya[j] - Yb[j];
    error += x * x;
  }
  if (error > 1e-3)
    fprintf(stderr, "error = %f\n", error);
  assert(error < 1e-3);
  free(X);
  free(Ya);
  free(Yb);
}
