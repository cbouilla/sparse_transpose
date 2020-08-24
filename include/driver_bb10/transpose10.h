///
/// \file transpose.h
/// \author Charles Bouillaguet and Jérôme Bonacchi
/// \brief //TODO
/// \date 2020-08-07
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_DRIVER_BB10_TRANSPOSE10_H
#define INCLUDE_DRIVER_BB10_TRANSPOSE10_H

#include <assert.h>
#include <err.h>
#include <stdint.h>
#include <stdlib.h>

#include "sparse.h"
#include "tools.h"

#define ENTRY_SIZE ((u8) (2 * sizeof(u32) + sizeof(double)))
///
/// \brief L1 cache line has size 64 on most CPUs
///
#define CACHELINE_SIZE ((u8)(64 / sizeof(u32))) ///< equal to 16
#define MSD_RADIX ((u8)7) ///< first pass
#define LSD_RADIX ((u8)9) ///< other passes
#define MAX_RADIX ((u8) 9) ///< maximum of the two radices above
#define L3_CACHE_SIZE ((u32) 31719424) // ((u32) (1 << 25)) 
#define MAX_PASSES ((u8) 8)

///
/// \brief A structure to store the information for the radix sort
///
struct ctx_t
{
  u8 bits;
  u8 n_passes; /* last pass is done... first, in parallel. */
  u8 radix[MAX_PASSES];
  u8 shift[MAX_PASSES];
  u32 n_buckets[MAX_PASSES];
  u32 pCOUNT[MAX_PASSES];
  u32 mask[MAX_PASSES];
  u32 *OUTi[MAX_PASSES];
  u32 *OUTj[MAX_PASSES];
  double *OUTx[MAX_PASSES];
  u32 seq_count_size;
  u32 par_count_size;
};

/* cache-resident buffer for (i, j) pairs. One such entry per output bucket.
   Invariants:  row[CACHELINE_SIZE - 1] contains COUNT[...] for this bucket,
        col[CACHELINE_SIZE - 1] contains offset of the first entry in this
   buffer.
     Actually 4 cache lines
   */
typedef struct
{
  u32 row[CACHELINE_SIZE];
  u32 col[CACHELINE_SIZE];
  double value[CACHELINE_SIZE];
} cacheline;

/* copy 64 bytes from src to dst using non-temporal store instructions
   if available (this bypasses the cache). */
static inline void store_nontemp_int(void *dst, const void *src);

/* copy 128 bytes from src to dst using non-temporal store instructions
   if available (this bypasses the cache). */
static inline void store_nontemp_double(void *dst, const void *src);

#if __AVX__
#include <immintrin.h>
static inline void store_nontemp_int(void *dst, const void *src)
{
  register __m256i *d1 = (__m256i *)dst;
  register __m256i s1 = *((__m256i *)src);
  register __m256i *d2 = d1 + 1;
  register __m256i s2 = *(((__m256i *)src) + 1);
  _mm256_stream_si256(d1, s1);
  _mm256_stream_si256(d2, s2);
  /* note : it can also be done using SSE for non-AVX machines */
}

static inline void store_nontemp_double(void *dst, const void *src)
{
  register __m256d *d1 = (__m256d *)dst;
  register __m256d s1 = *((__m256d *)src);
  register __m256d *d2 = d1 + 1;
  register __m256d s2 = *(((__m256d *)src) + 1);
  register __m256d *d3 = d1 + 2;
  register __m256d s3 = *(((__m256d *)src) + 2);
  register __m256d *d4 = d1 + 3;
  register __m256d s4 = *(((__m256d *)src) + 3);
  _mm256_stream_pd((double *)d1, s1);
  _mm256_stream_pd((double *)d2, s2);
  _mm256_stream_pd((double *)d3, s3);
  _mm256_stream_pd((double *)d4, s4);
  /* note : it can also be done using SSE for non-AVX machines */
}
#else
static inline void store_nontemp_int(void *dst, const void *src)
{
  const u32 *in = src;
  u32 *out = dst;
  for (u8 i = 0; i < CACHELINE_SIZE; i++)
    out[i] = in[i];
}

static inline void store_nontemp_double(void *dst, const void *src)
{
  const double *in = src;
  double *out = dst;
  for (u8 i = 0; i < CACHELINE_SIZE; i++)
    out[i] = in[i];
}
#endif

static inline void *malloc_aligned(const size_t size, const size_t alignment)
{
  void *x = aligned_alloc(alignment, size);
  if (x == NULL)
    err(1, "malloc failed");
  return x;
}

/* Allocate the buffer */
static cacheline *wc_alloc()
{
  return (cacheline *)malloc_aligned(
      sizeof(cacheline) * (1 << MAX_RADIX), 64);
}

/* Setup the buffer for a new pass */
static inline void wc_prime(cacheline *buffer, const u32 *COUNT,
                            const u32 n_buckets)
{
  for (u32 i = 0; i < n_buckets; i++)
  {
    buffer[i].row[CACHELINE_SIZE - 1] = COUNT[i];
    buffer[i].col[CACHELINE_SIZE - 1] = COUNT[i] & (CACHELINE_SIZE - 1);
  }
}

/* Transfer data stored in this buffer to the output arrays.
   Assumption: this bucket is filled to the end */
static inline void wc_flush(cacheline *self, const u32 count,
                            const u32 start, u32 *OUTi, u32 *OUTj, double *OUTx)
{
  u32 target = count & ~(CACHELINE_SIZE - 1);
  if (start != 0) // incomplete flush
  {
    for (u8 i = start; i < CACHELINE_SIZE; i++)
    {
      OUTi[target + i] = self->row[i];
      OUTj[target + i] = self->col[i];
      OUTx[target + i] = self->value[i];
    }
  }
  else // complete cache line flush
  {
    store_nontemp_int(OUTi + target, self->row);
    store_nontemp_int(OUTj + target, self->col);
    store_nontemp_double(OUTx + target, self->value);
  }
  self->col[CACHELINE_SIZE - 1] = 0;
}

/* push an (i,j) pair into the buffer */
static inline void wc_push(const u32 i, const u32 j, const double x,
                           cacheline *buffer, const u32 bucket_idx,
                           u32 *OUTi, u32 *OUTj, double *OUTx)
{
  cacheline *self = buffer + bucket_idx;
  u32 count = self->row[CACHELINE_SIZE - 1];
  u32 start = self->col[CACHELINE_SIZE - 1];
  u32 slot = count & (CACHELINE_SIZE - 1);
  self->row[slot] = i;
  self->col[slot] = j;
  self->value[slot] = x;
  if (slot == CACHELINE_SIZE - 1)
    wc_flush(self, count, start, OUTi, OUTj, OUTx);
  self->row[CACHELINE_SIZE - 1] = count + 1;
}

/* flush all buffer entries to the OUT arrays */
static inline void wc_purge(const cacheline *buffer,
                            const u32 n_buckets, u32 *OUTi, u32 *OUTj,
                            double *OUTx)
{
  for (u32 i = 0; i < n_buckets; i++)
  {
    const u32 count = buffer[i].row[CACHELINE_SIZE - 1];
    const u32 target = count & ~(CACHELINE_SIZE - 1);
    const u32 start = buffer[i].col[CACHELINE_SIZE - 1];
    for (u32 j = target + start; j < count; j++)
    {
      OUTi[j] = buffer[i].row[j - target];
      OUTj[j] = buffer[i].col[j - target];
      OUTx[j] = buffer[i].value[j - target];
    }
  }
}

/* prepare the ctx object with information needed for all passes */
void planification(struct ctx_t *ctx, mtx_CSR *R, u32 *scratch,
                   double *scratch2);

/* returns k such that buckets [0:k] are non-empty. */
u32 partitioning(struct ctx_t *ctx, const mtx_COO *A,
                 cacheline *buffer, u32 *tCOUNT, u32 *gCOUNT);

void histogram(const struct ctx_t *ctx, const u32 *Aj, const u32 lo,
               const u32 hi, const u8 n, u32 **W);

/* sequentially transpose a single bucket */
void transpose_bucket(struct ctx_t *ctx, cacheline *buffer,
                      const u32 lo, const u32 hi, mtx_CSR *R, const u32 bucket);

/* converts a sparse matrix in COOrdinate format to the CSR format.
   INPUT:  COO sparse matrix in Ai, Aj (both of size nnz), with n rows

   OUTPUT: CSR sparse matrix in Rp, Rj.

   Rp and Rj MUST be preallocated (of sizes n+1 and nnz, respectively).
   The "row pointers" Rp MUST be already computed.

   Ai, Aj, Rj MUST be aligned on a 64-byte boundary (for good cache behavior).
   The input arrays are expendable (i.e. they might be destroyed).
   The current code only reads them though. */
void transpose(const mtx_COO *A, mtx_CSR *R, const u32 num_threads);

#endif /* INCLUDE_DRIVER_BB10_TRANSPOSE10_H */
