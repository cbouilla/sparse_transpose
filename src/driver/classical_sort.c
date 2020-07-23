///
/// \file classical_sort.c
/// \author Charles Bouillaguet and Jérôme Bonacchi
/// \brief This file implements the classical algorithms to convert and to
/// tranpose matrices.
/// \date 2020-07-09
///
/// @copyright Copyright (c) 2020
///

#include "classical_sort.h"
#include "mini_spasm.h"

void classical_compress(const spasm_triplet *T, spasm *A, u32 *W)
{
  const u32 *Ti = T->i;
  const u32 *Tj = T->j;
  const double *Tx = T->x;
  const u32 n = T->n;
  const u32 nnz = T->nnz;

  u32 *Ap = A->p;
  u32 *Aj = A->j;
  double *Ax = A->x;

  // initializing W
  for (u32 i = 0; i < n; i++) // gcc replaces this with AVX2-optimized memset
    W[i] = 0;

  // counting entries on each row
  for (u32 k = 0; k < nnz; k++)
  {
    // W[Ti[k]]++
    u32 i = Ti[k];
    u32 count = W[i];
    count += 1;
    W[i] = count;
  }

  // prefix-sum on W and initializing Ap
  u32 sum = 0;
  for (u32 i = 0; i < n; i++)
  {
    // W[i], sum = sum, W[i] + sum
    u32 count = W[i];
    Ap[i] = sum;
    W[i] = sum;
    sum += count;
  }
  Ap[n] = sum;

  // dispatching
  for (u32 k = 0; k < nnz; k++)
  {
    // Aj[l = W[Ti[k]]++] = Tj[k]
    // Ax[l] = Tx[k]
    u32 i = Ti[k];
    u32 j = Tj[k];
    double x = Tx[k];
    u32 l = W[i];
    Aj[l] = j;
    Ax[l] = x;
    l += 1;
    W[i] = l;
  }
}

void classical_transpose(const spasm *A, spasm *R, u32 *W)
{
  const u32 *Ap = A->p;
  const u32 *Aj = A->j;
  const double *Ax = A->x;
  const u32 n = A->n;
  const u32 m = A->m;

  u32 *Rp = R->p;
  u32 *Rj = R->j;
  double *Rx = R->x;

  // initializing W
  for (u32 i = 0; i < m; i++) /* gcc replaces this with AVX2-optimized memset */
    W[i] = 0;

  // counting entries on each column
  for (u32 i = 0; i < n; i++)
    for (u32 k = Ap[i]; k < Ap[i + 1]; k++)
    {
      u32 j = Aj[k];
      u32 count = W[j];
      count += 1;
      W[j] = count;
    }

  // prefix-sum on W and initializing Rp
  u32 sum = 0;
  for (u32 j = 0; j < m; j++)
  {
    u32 count = W[j];
    Rp[j] = sum;
    W[j] = sum;
    sum += count;
  }
  Rp[m] = sum;

  // dispatching
  for (u32 i = 0; i < n; i++)
    for (u32 k = Ap[i]; k < Ap[i + 1]; k++)
    {
      u32 j = Aj[k];
      double x = Ax[k];
      u32 l = W[j];
      Rj[l] = i;
      Rx[l] = x;
      l += 1;
      W[j] = l;
    }
}

/* Wang et. al variant using the extra array. Faster than the classical
variant...

void wang_transpose(const spasm_triplet * A, spasm * T, int *W, int *Z)
{
        (void) W;
        const int *Ai = A->i;
        const int *Aj = A->j;
        int *Tp = T->p;
        int *Tj = T->j;
        int m = A->m;
        int nnz = A->nnz;

        for (int j = 0; j < m; j++)
                Tp[j] = 0;
        for (int i = 0; i < n; i++)
                for (int k = Ap[i]; k < Ap[i + 1]; k++) {
                        int j = Aj[k];
                        int w = Tp[j];
                        Z[k] = w;
                        w += 1;
                        Tp[j] = w;
                }

        int s = 0;
        for (int j = 0; j < m; j++) {
                int w = Tp[j];
                Tp[j] = s;
                s += w;
        }
        Tp[m] = s;

        for (int k = 0; k < nnz; k++) {
                int j = Aj[k];
                int s = Z[k];
                int r = Tp[j];
                int i = Ai[k];
                int l = r + s;
                Tj[l] = i;
        }
}
*/
