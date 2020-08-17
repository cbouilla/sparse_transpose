/* sparse matrix "transposition"

Copyright 2019 Charles Bouillaguet.

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sparse.h"
#include "transpose7.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void planification(struct ctx_t *ctx, mtx_CSR *R, u32 *scratch,
                   double *scratch2)
{
  u32 *Rj = R->j;
  double *Rx = R->x;
  const u32 Rn = R->n;
  const u32 nnz = R->nnz_max;
  u8 bits = 0;
  u32 tmp = Rn;
  while (tmp > 0)
  {
    bits++;
    tmp >>= 1;
  }
  printf("Rn: %d, bits: %d, ", R->n, bits);
  ctx->bits = bits;
  // TODO depends on the number of threads
  u8 n = 0;
  if (bits <= L3_CACHE_SIZE)
  {
    ctx->radix[n] = bits;
    ctx->shift[n] = bits;
    bits -= ctx->radix[n];
    n++;
  }
  else
  {
    ctx->radix[n] = L3_CACHE_SIZE;
    bits -= ctx->radix[n];
    ctx->shift[n] = bits;
    n++;
    u8 s_shift = 0;
    while (bits > 0)
    {
      ctx->shift[n] = s_shift;
      if (bits <= L2_CACHE_SIZE)
      {
        ctx->radix[n] = bits;
        bits -= ctx->radix[n];
        s_shift += ctx->radix[n];
        n++;
      }
      else
      {
        u32 i = 1;
        while (i <= L2_CACHE_SIZE && i <= bits)
          i++;
        i--;
        ctx->radix[n] = i;
        bits -= ctx->radix[n];
        s_shift += ctx->radix[n];
        n++;
      }
    }
  }
  assert(bits == 0);
  ctx->n_passes = n;

  printf("n: %d, ", n);
  printf("radix: %d", ctx->radix[0]);
  for (int p = 1; p < n; p++)
    printf(", %d", ctx->radix[p]);
  printf("\n");

  // u8 n = ceil((double)bits / MAX_RADIX_BITS);
  // ctx->bits = bits;
  // ctx->n_passes = n;

  // TODO depends on the number of threads
  /* the first pass is done on the most significant bits */
  // ctx->radix[0] = ceil((double)bits / n);
  // bits -= ctx->radix[0];
  // ctx->shift[0] = bits;

  /* other passes are done on least significant bits first */
  // u8 s_shift = 0;
  // for (int p = 1; p < n; p++)
  // {
  //   ctx->shift[p] = s_shift;
  //   u8 r = ceil((double)bits / (n - p)); /* bits in p-th pass */
  //   ctx->radix[p] = r;
  //   printf(", %d", r);
  //   bits -= r;
  //   s_shift += r;
  // }
  // assert(bits == 0);

  /* common to all passes */
  u32 s_count = 0;
  for (int p = 0; p < n; p++)
  {
    if (p == 1)
      s_count = 0;
    u32 size = 1 << ctx->radix[p];
    ctx->n_buckets[p] = size;
    ctx->pCOUNT[p] = s_count;
    s_count += size;
    ctx->mask[p] = size - 1;
    if ((n - p) & 1)
    {
      ctx->OUTi[p] = Rj;
      ctx->OUTj[p] = scratch;
      ctx->OUTx[p] = Rx;
    }
    else
    {
      ctx->OUTi[p] = scratch + ((nnz | 63) + 1);
      ctx->OUTj[p] = scratch + 2 * ((nnz | 63) + 1);
      ctx->OUTx[p] = scratch2;
    }

    /* check alignment: the hardcoded value of 63 corresponds to the
       L1 cache size on modern cpus */
    unsigned long check = (unsigned long)ctx->OUTi[p];
    assert((check & 63) == 0);
    check = (unsigned long)ctx->OUTj[p];
    assert((check & 63) == 0);
    check = (unsigned long)ctx->OUTx[p];
    assert((check & 63) == 0);
  }
  ctx->par_count_size = ctx->n_buckets[0];
  ctx->seq_count_size = s_count;
}

u32 partitioning(struct ctx_t *ctx, const mtx_COO *A, cacheline *buffer,
                 u32 *tCOUNT, u32 *gCOUNT)
{
  const u32 *Ai = A->i;
  const u32 *Aj = A->j;
  const double *Ax = A->x;
  const u32 nnz = A->nnz;
  const u32 mask = ctx->mask[0];
  const u32 size = ctx->n_buckets[0];
  const u8 shift = ctx->shift[0];
  u32 *OUTi = ctx->OUTi[0];
  u32 *OUTj = ctx->OUTj[0];
  double *OUTx = ctx->OUTx[0];
  assert(size == mask + 1);

  /* parallel partitioning */
  const u32 t = omp_get_thread_num();
  const u32 T = omp_get_num_threads();

  u32 *COUNT = tCOUNT + t * size;
  memset(COUNT, 0, size * sizeof(*COUNT));

#pragma omp for schedule(static)
  // TODO vectorize ?
  for (u32 k = 0; k < nnz; k++)
  {
    const u32 j = Aj[k];
    const u32 b = (j >> shift) & mask;
    COUNT[b]++;
  }

  u32 last = 0;
#pragma omp master
  {
    u32 sum = 0;
    for (u32 i = 0; i < size; i++)
    {
      gCOUNT[i] = sum;
      for (u32 t = 0; t < T; t++)
      {
        const u32 w = tCOUNT[t * size + i];
        tCOUNT[t * size + i] = sum;
        sum += w;
      }
      if (sum > gCOUNT[i])
        last = i;
    }
    gCOUNT[size] = sum;
    assert(sum == nnz);
  }

#pragma omp barrier

  wc_prime(buffer, COUNT, size);
#pragma omp for schedule(static) nowait
  for (u32 k = 0; k < nnz; k++)
  {
    const u32 i = Ai[k];
    const u32 j = Aj[k];
    const double x = Ax[k];
    const u32 b = (j >> shift) & mask;
    wc_push(i, j, x, buffer, b, OUTi, OUTj, OUTx);
  }
  wc_purge(buffer, size, OUTi, OUTj, OUTx);

  // correctness check
  //   #pragma omp barrier
  //   for (u32 b = 0; b < size; b++)
  //       for (u32 k = gCOUNT[b]; k < gCOUNT[b + 1]; k++) {
  //           u32 i = OUTi[k];
  //           u32 j = OUTj[k];
  //           u32 c = (j >> shift) & mask;
  //           if (c != b) {
  //               printf("discrepancy in bucket %d for index %d: found (%d, %d)
  //   which should be in bucket %d\n", b, k, i, j, c); assert(0);
  //           }
  //       }

  return last + 1;
}

void histogram(const struct ctx_t *ctx, const u32 *Aj, const u32 lo,
               const u32 hi, const u8 n, u32 **W)
{
  const u8 *shift = ctx->shift;
  const u32 *mask = ctx->mask;

  switch (n)
  {
  case 2:
    // TODO vectorize ?
    for (u32 k = lo; k < hi; k++)
    {
      const u32 j = Aj[k];
      const u32 q1 = j & mask[1];
      W[1][q1]++;
    }
    break;
  case 3:
    // TODO vectorize ?
    for (u32 k = lo; k < hi; k++)
    {
      const u32 j = Aj[k];
      const u32 q1 = j & mask[1];
      const u32 q2 = (j >> shift[2]) & mask[2];
      W[1][q1]++;
      W[2][q2]++;
    }
    break;
  case 4:
    // TODO vectorize ?
    for (u32 k = lo; k < hi; k++)
    {
      const u32 j = Aj[k];
      const u32 q1 = j & mask[1];
      const u32 q2 = (j >> shift[2]) & mask[2];
      const u32 q3 = (j >> shift[3]) & mask[3];
      W[1][q1]++;
      W[2][q2]++;
      W[3][q3]++;
    }
    break;
  default:
    err(1,
        "Ask the programmer to hardcode more passes in (radix) transpose...\n");
  }
}

void transpose_bucket(struct ctx_t *ctx, cacheline *buffer, const u32 lo,
                      const u32 hi, mtx_CSR *R, const u32 bucket)
{
  const u8 n = ctx->n_passes;
  const u32 csize = ctx->seq_count_size;
  u32 COUNT[csize];
  u32 *W[n];

  // TODO vectorize ?
  for (u8 p = 1; p < n; p++)
    W[p] = COUNT + ctx->pCOUNT[p];

  u32 *INi = ctx->OUTi[0];
  u32 *INj = ctx->OUTj[0];
  double *INx = ctx->OUTx[0];

  // printf("histograming [%d:%d] with %d passes\n", lo, hi, n);
  memset(COUNT, 0, csize * sizeof(*COUNT));
  histogram(ctx, INj, lo, hi, n, W);

  for (u8 p = 1; p < n; p++)
  {
    const u8 shift = ctx->shift[p];
    const u32 mask = ctx->mask[p];
    u32 *OUTi = ctx->OUTi[p];
    u32 *OUTj = ctx->OUTj[p];
    double *OUTx = ctx->OUTx[p];

    /* prefix-sum */
    u32 sum = lo;
    u32 size = ctx->n_buckets[p];
    for (u32 i = 0; i < size; i++)
    {
      u32 w = W[p][i];
      W[p][i] = sum;
      sum += w;
    }

    wc_prime(buffer, W[p], size);
    for (u32 k = lo; k < hi; k++)
    {
      const u32 i = INi[k];
      const u32 j = INj[k];
      const double x = INx[k];
      const u32 b = (j >> shift) & mask;
      wc_push(i, j, x, buffer, b, OUTi, OUTj, OUTx);
    }
    wc_purge(buffer, size, OUTi, OUTj, OUTx);

    /* check
    for (u32 b = 0; b < size - 1; b++)
        for (u32 k = W[p][b]; k < W[p][b + 1]; k++) {
            u32 i = OUTi[k];
            u32 j = OUTj[k];
            u32 c = (j >> shift) & mask;
            if (c != b) {
                printf("pass %d, discrepancy in bucket %d for index %d: found
    (%d, %d) which should be in bucket %d\n", p, b, k, i, j, c); ASSERT(0);
            }
        }
    */
    INi = OUTi;
    INj = OUTj;
    INx = OUTx;
  }

  // Computing row pointers
  const u32 tmp = 1 << ctx->shift[0];
  const u32 lo_bucket = bucket * tmp;
  const u32 hi_bucket = spasm_min(R->n + 1, lo_bucket + tmp);
  u32 ptr = lo;
  for (u32 i = lo_bucket; i < hi_bucket; ++i)
  {
    R->p[i] = ptr;
    while (ptr < hi && ctx->OUTj[n - 1][ptr] == i)
    {
      ptr++;
    }
  }
}

void transpose(const mtx_COO *A, mtx_CSR *R, const u32 num_threads)
{
  const u32 nnz = A->nnz;
  u32 *Rp = R->p;

  struct ctx_t ctx;
  u32 *scratch = (u32 *)malloc_aligned(3 * ((nnz | 63) + 1) * sizeof(u32), 64);
  double *scratch2 =
      (double *)malloc_aligned(((nnz | 63) + 1) * sizeof(double), 64);
  planification(&ctx, R, scratch, scratch2);

#ifdef BIG_BROTHER
  printf("$$$     transpose:\n");
  printf("$$$        bits: %d\n", ctx.bits);
  printf("$$$        passes: \n");
  for (int p = 0; p < ctx.n_passes; p++)
    printf("$$$         - %d \n", ctx.radix[p]);
#endif

  const u32 size = ctx.par_count_size;
  u32 T = 1;
#ifdef _OPENMP
  omp_set_num_threads(num_threads);
  T = omp_get_max_threads();
#endif // _OPENMP
  u32 tCOUNT[T * size];
  u32 gCOUNT[size + 1];

  u32 non_empty;

#pragma omp parallel
  {
    cacheline *buffer = wc_alloc();
#ifdef BIG_BROTHER
    double start = wct_seconds();
#endif
    u32 tmp = partitioning(&ctx, A, buffer, tCOUNT, gCOUNT);

#pragma omp master
    {
      non_empty = tmp;

#ifdef BIG_BROTHER
      printf("$$$        buckets: %d\n", non_empty);
      printf("$$$        partitioning-wct: %.2f\n", wct_seconds() - start);
      u32 smallest = nnz;
      u32 biggest = 0;
      for (int i = 0; i < non_empty; i++)
      {
        u32 size = gCOUNT[i + 1] - gCOUNT[i];
        smallest = MIN(smallest, size);
        biggest = MAX(biggest, size);
      }
      printf("$$$        smallest-bucket: %d\n", smallest);
      printf("$$$        avg-bucket: %" PRId64 "\n", nnz / non_empty);
      printf("$$$        biggest-bucket: %d\n", biggest);
#endif
    }

#pragma omp barrier
#ifdef BIG_BROTHER
    double sub_start = wct_seconds();
#endif

    if (ctx.n_passes == 1)
    {
      // Computing row pointers
#pragma omp for simd schedule(static)
      for (u32 i = 0; i < R->n + 1; i++)
      {
        Rp[i] = gCOUNT[i];
      }
    }
    else
    {
#pragma omp for schedule(dynamic, 1)
      for (u32 i = 0; i < non_empty; i++)
      {
        // TODO if [gCOUNT[i], gCOUNT[i + 1]] is too small, call another sorting
        // algorithm
        transpose_bucket(&ctx, buffer, gCOUNT[i], gCOUNT[i + 1], R, i);
      }
    }

    free(buffer);

#pragma omp master
    {
#ifdef BIG_BROTHER
      printf("$$$        buckets-wct: %.2f\n", wct_seconds() - sub_start);
#endif
    }
  } // omp parallel
  free(scratch);
  free(scratch2);
}
