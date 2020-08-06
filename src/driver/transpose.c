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

// #include "cado.h"

/* the following should come after cado.h, which sets -Werror=all */
// #ifdef __GNUC__
// #pragma GCC diagnostic ignored "-Wunknown-pragmas"
// #endif
// #include <stdio.h>
// #include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// #include "portability.h"
#include "mini_spasm.h"
#include "transpose.h"
// #include "typedefs.h"
// #include "utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #define BIG_BROTHER

/* OK, so we can do this the easy way, or the hard way. */

#ifdef TRANSPOSE_EASY_WAY
/* simple, cheap and dirty.
   This is similar to algorithm 2 in "Parallel Transposition of Sparse Data
   Structures", by Hao Wang, Weifeng Liu, Kaixi Hou, and Wu-chun Feng.

   It is in fact a direct parallelization of "distribution sorting" using atomic
   memory accesses. This works decently well, because there are very few
   conflicts.

   The row pointers MUST point to the END of each row */
void transpose(uint64_t nnz, index_t *Ai, index_t *Aj, index_t Rn, index_t *Rp,
               index_t *Ri)
{
  (void)Rn;
/* dispatch entries */
#pragma omp parallel for schedule(static)
  for (uint64_t k = 0; k < nnz; k++)
  {
    index_t i = Ai[k];
    index_t j = Aj[k];
    index_t s;
#pragma omp atomic capture
    s = --Rp[j];
    Ri[s] = i;
  }
}
#else
/* The hard way.
   Uses a parallel radix sort with a software write-combining buffer.
   Relies on aligned_alloc (OK in C11) and OpenMP. */

/* Allocate the buffer */
static struct cacheline_t *wc_alloc()
{
  return (struct cacheline_t *)malloc_aligned(
      sizeof(struct cacheline_t) * (1 << MAX_RADIX_BITS), 64);
}

/* Setup the buffer for a new pass */
static inline void wc_prime(struct cacheline_t *buffer, const index_t *COUNT,
                            int n_buckets)
{
  for (int i = 0; i < n_buckets; i++)
  {
    buffer[i].row[CACHELINE_SIZE - 1] =
        COUNT[i]; // TODO mettre 0 poour push dans l'ordre de gauche à droite et
                  // éviter les deux branches dans flush
    buffer[i].col[CACHELINE_SIZE - 1] = COUNT[i] & (CACHELINE_SIZE - 1);
  }
}

/* Transfer data stored in this buffer to the output arrays.
   Assumption: this bucket is filled to the end */
static inline void wc_flush(struct cacheline_t *self, index_t count,
                            index_t start, index_t *OUTi, index_t *OUTj, double *OUTx)
{
  index_t target = count & ~(CACHELINE_SIZE - 1);
  if (start != 0)
  { /* incomplete flush */
    for (int i = start; i < CACHELINE_SIZE; i++)
    {
      OUTi[target + i] = self->row[i];
      OUTj[target + i] = self->col[i];
      OUTx[target + i] = self->value[i];
    }
  }
  else
  { /* complete cache line flush */
    store_nontemp_int(OUTi + target, self->row);
    store_nontemp_int(OUTj + target, self->col);
    store_nontemp_double(OUTx + target, self->value);
  }
  self->col[CACHELINE_SIZE - 1] = 0;
}

/* push an (i,j) pair into the buffer */
static inline void wc_push(index_t i, index_t j, double x, struct cacheline_t *buffer,
                           index_t bucket_idx, index_t *OUTi, index_t *OUTj, double *OUTx)
{
  struct cacheline_t *self = buffer + bucket_idx;
  index_t count = self->row[CACHELINE_SIZE - 1];
  index_t start = self->col[CACHELINE_SIZE - 1];
  index_t slot = count & (CACHELINE_SIZE - 1);
  self->row[slot] = i;
  self->col[slot] = j;
  self->value[slot] = x;
  if (slot == CACHELINE_SIZE - 1)
    wc_flush(self, count, start, OUTi, OUTj, OUTx);
  self->row[CACHELINE_SIZE - 1] = count + 1;
}

/* flush all buffer entries to the OUT arrays */
static inline void wc_purge(struct cacheline_t *buffer, int n_buckets,
                            index_t *OUTi, index_t *OUTj, double *OUTx)
{
  for (int i = 0; i < n_buckets; i++)
  {
    index_t count = buffer[i].row[CACHELINE_SIZE - 1];
    index_t target = count & ~(CACHELINE_SIZE - 1);
    index_t start = buffer[i].col[CACHELINE_SIZE - 1];
    for (index_t j = target + start; j < count; j++)
    {
      OUTi[j] = buffer[i].row[j - target];
      OUTj[j] = buffer[i].col[j - target];
      OUTx[j] = buffer[i].value[j - target];
    }
  }
}

/* Setup the buffer for a new pass */
static inline void wc_half_prime(struct half_cacheline_t *buffer, char *start,
                                 const index_t *COUNT, int n_buckets)
{
  for (int i = 0; i < n_buckets; i++)
  {
    buffer[i].row[CACHELINE_SIZE - 1] = COUNT[i];
    start[i] = COUNT[i] & (CACHELINE_SIZE - 1);
  }
}

/* Transfer data stored in this buffer to the output arrays.
   Assumption: this buckget is filled to the end */
static inline void wc_half_flush(struct half_cacheline_t *self, index_t count,
                                 char start, index_t *OUTi, double *OUTx)
{
  index_t target = count & ~(CACHELINE_SIZE - 1);
  if (start == 0)
  { /* complete cache line flush */
    store_nontemp_int(OUTi + target, self->row);
    store_nontemp_double(OUTx + target, self->value);
  }
  else
  { /* incomplete flush */
    for (int i = start; i < CACHELINE_SIZE; i++)
    {
      OUTi[target + i] = self->row[i];
      OUTx[target + i] = self->value[i];
    }
  }
}

/* push an (i,j) pair into the buffer */
static inline void wc_half_push(index_t i, double x, struct half_cacheline_t *buffer,
                                char *start, index_t bucket_idx, index_t *OUTi, double *OUTx)
{
  struct half_cacheline_t *self = buffer + bucket_idx;
  index_t count = self->row[CACHELINE_SIZE - 1];
  index_t slot = count & (CACHELINE_SIZE - 1);
  self->row[slot] = i;
  self->value[slot] = x;
  if (slot == CACHELINE_SIZE - 1)
  {
    wc_half_flush(self, count, start[bucket_idx], OUTi, OUTx);
    start[bucket_idx] = 0;
  }
  self->row[CACHELINE_SIZE - 1] = count + 1;
}

/* flush all buffer entries to the OUT arrays */
static inline void wc_half_purge(struct half_cacheline_t *buffer, char *start,
                                 int n_buckets, index_t *OUTi, double *OUTx)
{
  for (int i = 0; i < n_buckets; i++)
  {
    index_t count = buffer[i].row[CACHELINE_SIZE - 1];
    index_t target = count & ~(CACHELINE_SIZE - 1);
    index_t start_ptr = start[i];
    for (index_t j = target + start_ptr; j < count; j++)
    {
      OUTi[j] = buffer[i].row[j - target];
      OUTx[j] = buffer[i].value[j - target];
    }
  }
}

/* prepare the ctx object with information needed for all passes */
static void planification(struct ctx_t *ctx, spasm *R, index_t *scratch, double *scratch2)
{
  index_t Rn = R->n;
  index_t nnz = R->nnz_max;
  index_t *Ri = R->j;
  double *Rx = R->x;
  int bits = 0;
  int tmp = Rn;
  while (tmp > 0)
  {
    bits++;
    tmp >>= 1;
  }
  int n = ceil((double)bits / MAX_RADIX_BITS);
  ctx->bits = bits;
  ctx->n_passes = n;

  /* the first pass is done on the most significant bits */
  ctx->radix[0] = ceil((double)bits / n);
  bits -= ctx->radix[0];
  ctx->shift[0] = bits;

  /* other passes are done on least significant bits first */
  int s_shift = 0;
  for (int p = 1; p < n; p++)
  {
    ctx->shift[p] = s_shift;
    int r = ceil((double)bits / (n - p)); /* bits in p-th pass */
    ctx->radix[p] = r;
    bits -= r;
    s_shift += r;
  }
  assert(bits == 0);

  /* common to all passes */
  int s_count = 0;
  for (int p = 0; p < n; p++)
  {
    if (p == 1)
      s_count = 0;
    int size = 1 << ctx->radix[p];
    ctx->n_buckets[p] = size;
    ctx->pCOUNT[p] = s_count;
    s_count += size;
    ctx->mask[p] = size - 1;
    if ((n - p) & 1)
    {
      ctx->OUTi[p] = Ri;
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

/* returns k such that buckets [0:k] are non-empty. */
static int partitioning(const struct ctx_t *ctx, spasm_triplet *A,
                        struct cacheline_t *buffer, index_t *tCOUNT,
                        index_t *gCOUNT)
{
  const index_t *Ai = A->i;
  const index_t *Aj = A->j;
  const double *Ax = A->x;
  index_t nnz = A->nnz;
  index_t mask = ctx->mask[0];
  index_t size = ctx->n_buckets[0];
  index_t shift = ctx->shift[0];
  index_t *OUTi = ctx->OUTi[0];
  index_t *OUTj = ctx->OUTj[0];
  double *OUTx = ctx->OUTx[0];
  assert(size == mask + 1);

  /* parallel partitioning */
  int t = omp_get_thread_num();
  int T = omp_get_num_threads();

  index_t *COUNT = tCOUNT + t * size;
  memset(COUNT, 0, size * sizeof(*COUNT));

#pragma omp for schedule(static)
  for (index_t k = 0; k < nnz; k++)
  {
    index_t j = Aj[k];
    index_t b = (j >> shift) & mask;
    COUNT[b]++;
  }

  index_t last = 0;
// TODO utiliser omp single ?
#pragma omp master
  {
    index_t s = 0;
    for (index_t i = 0; i < size; i++)
    {
      gCOUNT[i] = s;
      for (int t = 0; t < T; t++)
      {
        index_t w = tCOUNT[t * size + i];
        tCOUNT[t * size + i] = s;
        s += w;
      }
      if (s > gCOUNT[i])
        last = i;
    }
    gCOUNT[size] = s;
    assert(s == nnz);
  }

#pragma omp barrier

  wc_prime(buffer, COUNT, size);
#pragma omp for schedule(static) nowait
  for (index_t k = 0; k < nnz; k++)
  {
    index_t i = Ai[k];
    index_t j = Aj[k];
    double x = Ax[k];
    index_t b = (j >> shift) & mask;
    wc_push(i, j, x, buffer, b, OUTi, OUTj, OUTx);
  }
  wc_purge(buffer, size, OUTi, OUTj, OUTx);

  // correctness check
  //   #pragma omp barrier
  //   for (index_t b = 0; b < size; b++)
  //       for (index_t k = gCOUNT[b]; k < gCOUNT[b + 1]; k++) {
  //           index_t i = OUTi[k];
  //           index_t j = OUTj[k];
  //           index_t c = (j >> shift) & mask;
  //           if (c != b) {
  //               printf("discrepancy in bucket %d for index %d: found (%d, %d)
  //   which should be in bucket %d\n", b, k, i, j, c); assert(0);
  //           }
  //       }

  return last + 1;
}

static void histogram(struct ctx_t *ctx, const index_t *Aj, index_t lo,
                      index_t hi, int n, index_t **W)
{
  index_t *shift = ctx->shift;
  index_t *mask = ctx->mask;

  switch (n)
  {
  case 2:
    for (index_t k = lo; k < hi; k++)
    {
      index_t j = Aj[k];
      int q1 = j & mask[1];
      W[1][q1]++;
    }
    break;
  case 3:
    for (index_t k = lo; k < hi; k++)
    {
      index_t j = Aj[k];
      int q1 = j & mask[1];
      int q2 = (j >> shift[2]) & mask[2];
      W[1][q1]++;
      W[2][q2]++;
    }
    break;
  case 4:
    for (index_t k = lo; k < hi; k++)
    {
      index_t j = Aj[k];
      int q1 = j & mask[1];
      int q2 = (j >> shift[2]) & mask[2];
      int q3 = (j >> shift[3]) & mask[3];
      W[1][q1]++;
      W[2][q2]++;
      W[3][q3]++;
    }
    break;
#if __SIZEOF_INDEX__ == 8
  case 5:
#pragma omp for schedule(static)
    for (index_t k = lo; k < hi; k++)
    {
      index_t j = Aj[k];
      int q1 = j & mask[1];
      int q2 = (j >> shift[2]) & mask[2];
      int q3 = (j >> shift[3]) & mask[3];
      int q4 = (j >> shift[4]) & mask[4];
      W[1][q1]++;
      W[2][q2]++;
      W[3][q3]++;
      W[4][q4]++;
    }
    break;
#endif
  default:
    err(1,
        "Ask the programmer to hardcode more passes in (radix) transpose...\n");
  }
}

/* sequentially transpose a single bucket */
void transpose_bucket(struct ctx_t *ctx, struct cacheline_t *buffer, index_t lo,
                      index_t hi)
{
  int n = ctx->n_passes;
  index_t csize = ctx->seq_count_size;
  index_t COUNT[csize];
  index_t *W[n];

  for (int p = 1; p < n; p++)
    W[p] = COUNT + ctx->pCOUNT[p];

  index_t *INi = ctx->OUTi[0];
  index_t *INj = ctx->OUTj[0];
  double *INx = ctx->OUTx[0];

  // printf("histograming [%d:%d] with %d passes\n", lo, hi, n);
  memset(COUNT, 0, csize * sizeof(*COUNT));
  histogram(ctx, INj, lo, hi, n, W);

  for (int p = 1; p < n; p++)
  {
    index_t shift = ctx->shift[p];
    index_t mask = ctx->mask[p];
    index_t *OUTi = ctx->OUTi[p];
    index_t *OUTj = ctx->OUTj[p];
    double *OUTx = ctx->OUTx[p];

    /* prefix-sum */
    index_t s = lo;
    index_t size = ctx->n_buckets[p];
    for (index_t i = 0; i < size; i++)
    {
      index_t w = W[p][i];
      W[p][i] = s;
      s += w;
    }

    if (p < n - 1)
    {
      /* full pass */
      wc_prime(buffer, W[p], size);
      for (index_t k = lo; k < hi; k++)
      {
        index_t i = INi[k];
        index_t j = INj[k];
        double x = INx[k];
        int b = (j >> shift) & mask;
        wc_push(i, j, x, buffer, b, OUTi, OUTj, OUTx);
      }
      wc_purge(buffer, size, OUTi, OUTj, OUTx);
    }
    else
    {
      // TODO pourquoi
      /* last pass: we don't need to write the j values */
      struct half_cacheline_t *half_buffer = (struct half_cacheline_t *)buffer;
      char start[size];
      wc_half_prime(half_buffer, start, W[p], size);
      for (index_t k = lo; k < hi; k++)
      {
        index_t i = INi[k];
        index_t j = INj[k];
        double x = INx[k];
        int b = (j >> shift) & mask;
        wc_half_push(i, x, half_buffer, start, b, OUTi, OUTx);
      }
      wc_half_purge(half_buffer, start, size, OUTi, OUTx);
    }

    /* check
    for (index_t b = 0; b < size - 1; b++)
        for (index_t k = W[p][b]; k < W[p][b + 1]; k++) {
            index_t i = OUTi[k];
            index_t j = OUTj[k];
            index_t c = (j >> shift) & mask;
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
}

void transpose(spasm_triplet *A, spasm *R)
{
  u32 nnz = A->nnz;
  index_t *Ai = A->i;
  index_t *Aj = A->j;
  double *Ax = A->x;
  index_t Rn = R->n;
  index_t Rm = R->m;
  index_t *Rp = R->p;
  index_t *Ri = R->j;
  double *Rx = R->x;
  (void)Rp;

  struct ctx_t ctx;
  index_t *scratch =
      (index_t *)malloc_aligned(3 * ((nnz | 63) + 1) * sizeof(index_t), 64);
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

  int size = ctx.par_count_size;
  int T = omp_get_max_threads();
  index_t tCOUNT[T * size];
  index_t gCOUNT[size + 1];

  int non_empty;

#pragma omp parallel
  {
    struct cacheline_t *buffer = wc_alloc();
#ifdef BIG_BROTHER
    double start = wct_seconds();
#endif
    int tmp = partitioning(&ctx, A, buffer, tCOUNT, gCOUNT);

#pragma omp master
    {
      non_empty = tmp;

#ifdef BIG_BROTHER
      printf("$$$        buckets: %d\n", non_empty);
      printf("$$$        partitioning-wct: %.2f\n", wct_seconds() - start);
      index_t smallest = nnz;
      index_t biggest = 0;
      for (int i = 0; i < non_empty; i++)
      {
        index_t size = gCOUNT[i + 1] - gCOUNT[i];
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

    if (ctx.n_passes > 1)
#pragma omp for schedule(dynamic, 1)
      for (int i = 0; i < non_empty; i++)
      {
        transpose_bucket(&ctx, buffer, gCOUNT[i], gCOUNT[i + 1]);
      }
    free_aligned(buffer);

#pragma omp master
    {
#ifdef BIG_BROTHER
      printf("$$$        buckets-wct: %.2f\n", wct_seconds() - sub_start);
#endif
    }
  }
  free_aligned(scratch);
  free_aligned(scratch2);

  if (ctx.n_passes == 1)
  {
    for (int i = 0; i <= non_empty; i++)
    {
      Rp[i] = gCOUNT[i];
    }
  }
  else
  {
    u32 length = spasm_max(Rn, Rm) + 1;
    index_t *W = (index_t *)malloc_aligned(length * sizeof(*W), 64);
    for (u32 i = 0; i < length; ++i)
      W[i] = 0;
    for (index_t k = 0; k < nnz; k++)
    {
      // W[Ti[k]]++
      index_t i = Aj[k];
      index_t count = W[i];
      count += 1;
      W[i] = count;
    }

    // prefix-sum on W and initializing Ap
    index_t sum = 0;
    for (index_t i = 0; i < Rn; i++)
    {
      // W[i], sum = sum, W[i] + sum
      index_t count = W[i];
      Rp[i] = sum;
      W[i] = sum;
      sum += count;
    }
    Rp[Rn] = sum;
    free_aligned(W);
  }
}
#endif
