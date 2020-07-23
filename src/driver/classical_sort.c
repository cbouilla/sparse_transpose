///
/// \file classical_sort.c
/// \author Charles Bouillaguet & Jérôme Bonacchi
/// \brief Implementation of the "classical" sort algorithm from Gustavson.
/// \date 2020-07-09
///
/// @copyright Copyright (c) 2020
///

#include "classical_sort.h"
#include "mini_spasm.h"

///
/// \brief Converts a matrix in triplet format into a matrix in CSR formatusing
/// Gustavon's algorithm.
///
/// \param T Input matrix in COO format
/// \param A Output matrix in CSR format
/// \param W Scratch space, size == #columns + 1
///
void classical_compress(const spasm_triplet *T, spasm *A, int *W)
{
  const int *Ti = T->i;
  const int *Tj = T->j;
  const double *Tx = T->x;
  const int n = T->n;
  const int nnz = T->nnz;

  int *Ap = A->p;
  int *Aj = A->j;
  double *Ax = A->x;

  // initializing W
  for (int i = 0; i < n; i++) // gcc replaces this with AVX2-optimized memset
    W[i] = 0;

  // counting entries on each row
  for (int k = 0; k < nnz; k++)
  {
    // W[Ti[k]]++
    int i = Ti[k];
    int count = W[i];
    count += 1;
    W[i] = count;
  }

  // prefix-sum on W and initializing Ap
  int sum = 0;
  for (int i = 0; i < n; i++)
  {
    // W[i], sum = sum, W[i] + sum
    int count = W[i];
    Ap[i] = sum;
    W[i] = sum;
    sum += count;
  }
  Ap[n] = sum;

  // dispatching
  for (int k = 0; k < nnz; k++)
  {
    // Aj[l = W[Ti[k]]++] = Tj[k]
    // Ax[l] = Tx[k]
    int i = Ti[k];
    int j = Tj[k];
    double x = Tx[k];
    int l = W[i];
    Aj[l] = j;
    Ax[l] = x;
    l += 1;
    W[i] = l;
  }
}

///
/// \brief Transposes a matrix in CSR format using Gustavon's algorithm.
///
/// \param A Input matrix in CSR format
/// \param R Output matrix in CSR format
/// \param W Scratch space, size == #columns + 1
///
void classical_transpose(const spasm *A, spasm *R, int *W)
{
  const int *Ap = A->p;
  const int *Aj = A->j;
  const double *Ax = A->x;
  const int n = A->n;
  const int m = A->m;

  int *Rp = R->p;
  int *Rj = R->j;
  double *Rx = R->x;

  // initializing W
  for (int i = 0; i < m; i++) /* gcc replaces this with AVX2-optimized memset */
    W[i] = 0;

  // counting entries on each column
  for (int i = 0; i < n; i++)
    for (int k = Ap[i]; k < Ap[i + 1]; k++)
    {
      int j = Aj[k];
      int count = W[j];
      count += 1;
      W[j] = count;
    }

  // prefix-sum on W and initializing Rp
  int sum = 0;
  for (int j = 0; j < m; j++)
  {
    int count = W[j];
    Rp[j] = sum;
    W[j] = sum;
    sum += count;
  }
  Rp[m] = sum;

  // dispatching
  for (int i = 0; i < n; i++)
    for (int k = Ap[i]; k < Ap[i + 1]; k++)
    {
      int j = Aj[k];
      double x = Ax[k];
      int l = W[j];
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
