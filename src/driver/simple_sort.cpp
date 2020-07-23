///
/// \file simple_sort.cpp
/// \author Charles Bouillaguet and Jérôme Bonacchi
/// \brief This files uses std::sort and tbb::parallel_sort algorithms to
/// convert and to tranpose matrices.
/// \date 2020-07-09
///
/// @copyright Copyright (c) 2020
///

#include <algorithm>
#include <cstdint>
#include <omp.h>

extern "C"
{
#include "mini_spasm.h"
}

extern "C" void stdsort_compress(const spasm_triplet *T, spasm *A,
                                 struct matrix_entry_t *Te);
extern "C" void stdsort_transpose(const spasm *A, spasm *R,
                                  struct matrix_entry_t *Te);

#ifdef HAVE_TBB
#include <iostream>
#include <tbb/parallel_sort.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tbb_stddef.h>

extern "C" void tbbsort_compress(const spasm_triplet *T, spasm *A,
                                 struct matrix_entry_t *Te, int num_threads);
extern "C" void tbbsort_transpose(const spasm *A, spasm *R,
                                  struct matrix_entry_t *Te, int num_threads);
extern "C" void tbb_version();

#endif // HAVE_TBB

///
/// \brief Compares the row of the matrix_entry_t data type.
///
/// \return true if the row index of `a` is lesser than the row index of b
/// \return false if the row index of `a` is greater than, or equal to, the row
/// index of b
///
inline bool operator<(const struct matrix_entry_t &a,
                      const struct matrix_entry_t &b)
{
  return a.i < b.i;
}

void finalize(const u32 n, const u32 nnz, const matrix_entry_t *Te,
                     spasm *A)
{
  u32 *Ap = A->p;
  u32 *Aj = A->j;
  double *Ax = A->x;

  // scanning, copying and setting pointers
  // could be parallelized (but this is not completely trivial)
  int i = -1;
  for (u32 k = 0; k < nnz; k++)
  {
    u32 next_i = Te[k].i;
    u32 next_j = Te[k].j;
    double next_x = Te[k].x;
    if (i < next_i)
    {
      i++;
      while (i < next_i)
      {
        Ap[i] = k;
        i++;
      }
      Ap[i] = k;
    }
    // append entry to current line
    Aj[k] = next_j;
    Ax[k] = next_x;
  }
  // finalization
  i += 1;
  while (i <= n)
  {
    Ap[i] = nnz;
    i += 1;
  }
}

void stdsort_compress(const spasm_triplet *T, spasm *A,
                      struct matrix_entry_t *Te)
{
  u32 n = T->n;
  u32 nnz = T->nnz;
  u32 *Ti = T->i;
  u32 *Tj = T->j;
  double *Tx = T->x;

  for (u32 k = 0; k < nnz; k++)
  {
    Te[k].i = Ti[k];
    Te[k].j = Tj[k];
    Te[k].x = Tx[k];
  }

  std::sort(Te, Te + nnz);
  finalize(n, nnz, Te, A);
}

void stdsort_transpose(const spasm *A, spasm *R, struct matrix_entry_t *Te)
{
  u32 n = A->n;
  u32 m = A->m;
  u32 nnz = spasm_nnz(A);
  const u32 *Ap = A->p;
  const u32 *Aj = A->j;
  const double *Ax = A->x;

  for (u32 i = 0; i < n; i++)
    for (u32 k = Ap[i]; k < Ap[i + 1]; k++)
    {
      u32 j = Aj[k];
      double x = Ax[k];
      Te[k].i = j;
      Te[k].j = i;
      Te[k].x = x;
    }

  std::sort(Te, Te + nnz);
  finalize(m, nnz, Te, R);
}

#ifdef HAVE_TBB

void tbbsort_compress(const spasm_triplet *T, spasm *A,
                      struct matrix_entry_t *Te, const u32 num_threads)
{
  tbb::task_scheduler_init tsi(num_threads);

  u32 n = T->n;
  u32 nnz = T->nnz;
  u32 *Ti = T->i;
  u32 *Tj = T->j;
  double *Tx = T->x;

  // double a = spasm_wtime();
  omp_set_num_threads(num_threads);
#pragma omp parallel for
  for (u32 k = 0; k < nnz; k++)
  {
    Te[k].i = Ti[k];
    Te[k].j = Tj[k];
    Te[k].x = Tx[k];
  }
  // double b = spasm_wtime();
  tbb::parallel_sort(Te, Te + nnz);
  // double c = spasm_wtime();
  finalize(n, nnz, Te, A);

  // extra timing
  // double d = spasm_wtime();
  // printf("      Subtimes:\n");
  // printf("        Fill_Te: %.3f\n", b-a);
  // printf("        tbb::sort: %.3f\n", c-b);
  // printf("        finalize: %.3f\n", d-c);
}

void tbbsort_transpose(const spasm *A, spasm *R, struct matrix_entry_t *Te,
                       const u32 num_threads)
{
  tbb::task_scheduler_init tsi(num_threads);

  u32 n = A->n;
  u32 m = A->m;
  u32 nnz = spasm_nnz(A);
  const u32 *Ap = A->p;
  const u32 *Aj = A->j;
  const double *Ax = A->x;

  // double a = spasm_wtime();
  omp_set_num_threads(num_threads);
#pragma omp parallel for
  for (u32 i = 0; i < n; i++)
    for (u32 k = Ap[i]; k < Ap[i + 1]; k++)
    {
      u32 j = Aj[k];
      double x = Ax[k];
      Te[k].j = i;
      Te[k].i = j;
      Te[k].x = x;
    }

  // double b = spasm_wtime();
  tbb::parallel_sort(Te, Te + nnz);
  // double c = spasm_wtime();
  finalize(m, nnz, Te, R);

  // extra timing
  // double d = spasm_wtime();
  // printf("      Subtimes:\n");
  // printf("        Fill_Te: %.3f\n", b-a);
  // printf("        tbb::sort: %.3f\n", c-b);
  // printf("        finalize: %.3f\n", d-c);
}

void tbb_version()
{
  std::cout << "TBB version: compiled=" << TBB_INTERFACE_VERSION
            << ", runtime=" << tbb::TBB_runtime_interface_version() << '\n';
}

#endif // HAVE_TBB
