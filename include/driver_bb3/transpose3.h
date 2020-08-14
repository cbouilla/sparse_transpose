///
/// \file transpose.h
/// \author Charles Bouillaguet and Jérôme Bonacchi
/// \brief //TODO
/// \date 2020-08-07
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_DRIVER_BB3_TRANSPOSE3_H
#define INCLUDE_DRIVER_BB3_TRANSPOSE3_H

#include <assert.h>
#include <err.h>
#include <stdint.h>
#include <stdlib.h>

#include "sparse.h"
#include "tools.h"

///
/// \brief L1 cache line has size 64 on most CPUs
///
#define CACHELINE_SIZE ((u8)(64 / sizeof(mtx_entry)))
#define MAX_RADIX_BITS 10 ///< was experimentally found to be OK
#define MAX_PASSES 4

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
  // u32 *OUTi[MAX_PASSES];
  // u32 *OUTj[MAX_PASSES];
  // double *OUTx[MAX_PASSES];
  mtx_entry *OUT[MAX_PASSES];
  u32 seq_count_size;
  u32 par_count_size;
};

/* cache-resident buffer for (i, j) pairs. One such entry per output bucket.
   Invariants:  row[CACHELINE_SIZE - 1] contains COUNT[...] for this bucket,
        col[CACHELINE_SIZE - 1] contains offset of the first entry in this
   buffer.
   */
typedef mtx_entry cacheline[CACHELINE_SIZE];

/* copy 64 bytes from src to dst using non-temporal store instructions
   if available (this bypasses the cache). */
static inline void store_nontemp_cacheline(void *dst, const void *src);

#if __AVX__
#include <immintrin.h>
static inline void store_nontemp_cacheline(void *dst, const void *src)
{
  register __m256i *d1 = (__m256i *)dst;
  register __m256i s1 = *((__m256i *)src);
  register __m256i *d2 = d1 + 1;
  register __m256i s2 = *(((__m256i *)src) + 1);
  _mm256_stream_si256(d1, s1);
  _mm256_stream_si256(d2, s2);
  /* note : it can also be done using SSE for non-AVX machines */
}
#else
static inline void store_nontemp_cacheline(void *dst, const void *src)
{
  const mtx_entry *in = src;
  mtx_entry *out = dst;
  // TODO use memcpy ?
  for (u8 i = 0; i < CACHELINE_SIZE; i++)
  {
    out[i].i = in[i].i;
    out[i].j = in[i].j;
    out[i].x = in[i].x;
  }
}

#endif

static inline void *malloc_aligned(const size_t size, const size_t alignment)
{
  void *x = aligned_alloc(alignment, size);
  if (x == NULL)
    err(1, "malloc failed");
  return x;
}

/* Allocate the buffer */ // TODO allouer moins de mémoire en prenant un
                          // argument
static cacheline *wc_alloc()
{
  return (cacheline *)malloc_aligned(sizeof(cacheline) * (1 << MAX_RADIX_BITS),
                                     64);
}

/* Setup the buffer for a new pass */
static inline void wc_prime(cacheline *buffer, const u32 *COUNT,
                            const u32 n_buckets)
{
  for (u32 i = 0; i < n_buckets; i++)
  {
    buffer[i][CACHELINE_SIZE - 1].i = COUNT[i];
    buffer[i][CACHELINE_SIZE - 1].j = COUNT[i] & (CACHELINE_SIZE - 1);
  }
}

/* Transfer data stored in this buffer to the output arrays.
   Assumption: this bucket is filled to the end */
static inline void wc_flush(cacheline *self, const u32 count, const u32 start,
                            mtx_entry *OUT)
{
  u32 target = count & ~(CACHELINE_SIZE - 1);
  if (start != 0)
  { /* incomplete flush */
    for (u8 i = start; i < CACHELINE_SIZE; i++)
    {
      OUT[target + i].i = (*self)[i].i;
      OUT[target + i].j = (*self)[i].j;
      OUT[target + i].x = (*self)[i].x;
    }
  }
  else
  { /* complete cache line flush */
    store_nontemp_cacheline(OUT + target, self);
  }
  (*self)[CACHELINE_SIZE - 1].j = 0;
}

/* push an (i,j) pair into the buffer */
static inline void wc_push(const mtx_entry *entry, cacheline *buffer,
                           const u32 bucket_idx, mtx_entry *OUT)
{
  cacheline *self = buffer + bucket_idx;
  u32 count = (*self)[CACHELINE_SIZE - 1].i;
  u32 start = (*self)[CACHELINE_SIZE - 1].j;
  u32 slot = count & (CACHELINE_SIZE - 1);
  (*self)[slot] = *entry;
  if (slot == CACHELINE_SIZE - 1)
    wc_flush(self, count, start, OUT);
  (*self)[CACHELINE_SIZE - 1].i = count + 1;
}

/* flush all buffer entries to the OUT arrays */
static inline void wc_purge(const cacheline *buffer, const u32 n_buckets,
                            mtx_entry *OUT)
{
  for (u32 i = 0; i < n_buckets; i++)
  {
    const u32 count = buffer[i][CACHELINE_SIZE - 1].i;
    const u32 target = count & ~(CACHELINE_SIZE - 1);
    const u32 start = buffer[i][CACHELINE_SIZE - 1].j;
    for (u32 j = target + start; j < count; j++)
    {
      OUT[j] = buffer[i][j - target];
    }
  }
}

static inline void wc_flush_last(cacheline *self, const u32 count,
                                 const u32 start, u32 *OUTi, u32 *OUTj,
                                 double *OUTx)
{
  u32 target = count & ~(CACHELINE_SIZE - 1);
  for (u8 i = start; i < CACHELINE_SIZE; i++)
  {
    OUTi[target + i] = (*self)[i].i;
    OUTj[target + i] = (*self)[i].j;
    OUTx[target + i] = (*self)[i].x;
  }
  (*self)[CACHELINE_SIZE - 1].j = 0;
}

/* push an (i,j) pair into the buffer */
static inline void wc_push_last(const mtx_entry *entry, cacheline *buffer,
                                const u32 bucket_idx, u32 *OUTi, u32 *OUTj,
                                double *OUTx)
{
  cacheline *self = buffer + bucket_idx;
  u32 count = (*self)[CACHELINE_SIZE - 1].i;
  u32 start = (*self)[CACHELINE_SIZE - 1].j;
  u32 slot = count & (CACHELINE_SIZE - 1);
  (*self)[slot] = *entry;
  if (slot == CACHELINE_SIZE - 1)
    wc_flush_last(self, count, start, OUTi, OUTj, OUTx);
  (*self)[CACHELINE_SIZE - 1].i = count + 1;
}

/* flush all buffer entries to the OUT arrays */
static inline void wc_purge_last(const cacheline *buffer, const u32 n_buckets,
                                 u32 *OUTi, u32 *OUTj, double *OUTx)
{
  for (u32 i = 0; i < n_buckets; i++)
  {
    const u32 count = buffer[i][CACHELINE_SIZE - 1].i;
    const u32 target = count & ~(CACHELINE_SIZE - 1);
    const u32 start = buffer[i][CACHELINE_SIZE - 1].j;
    for (u32 j = target + start; j < count; j++)
    {
      OUTi[j] = buffer[i][j - target].i;
      OUTj[j] = buffer[i][j - target].j;
      OUTx[j] = buffer[i][j - target].x;
    }
  }
}

static inline void wc_flush_only1(cacheline *self, const u32 count,
                                  const u32 start, u32 *OUTi, double *OUTx)
{
  u32 target = count & ~(CACHELINE_SIZE - 1);
  for (u8 i = start; i < CACHELINE_SIZE; i++)
  {
    OUTi[target + i] = (*self)[i].i;
    OUTx[target + i] = (*self)[i].x;
  }
  (*self)[CACHELINE_SIZE - 1].j = 0;
}

/* push an (i,j) pair into the buffer */
static inline void wc_push_only1(const u32 i, const double x, cacheline *buffer,
                                 const u32 bucket_idx, u32 *OUTi, double *OUTx)
{
  cacheline *self = buffer + bucket_idx;
  u32 count = (*self)[CACHELINE_SIZE - 1].i;
  u32 start = (*self)[CACHELINE_SIZE - 1].j;
  u32 slot = count & (CACHELINE_SIZE - 1);
  (*self)[slot].i = i;
  (*self)[slot].x = x;
  if (slot == CACHELINE_SIZE - 1)
    wc_flush_only1(self, count, start, OUTi, OUTx);
  (*self)[CACHELINE_SIZE - 1].i = count + 1;
}

/* flush all buffer entries to the OUT arrays */
static inline void wc_purge_only1(const cacheline *buffer, const u32 n_buckets,
                                  u32 *OUTi, double *OUTx)
{
  for (u32 i = 0; i < n_buckets; i++)
  {
    const u32 count = buffer[i][CACHELINE_SIZE - 1].i;
    const u32 target = count & ~(CACHELINE_SIZE - 1);
    const u32 start = buffer[i][CACHELINE_SIZE - 1].j;
    for (u32 j = target + start; j < count; j++)
    {
      OUTi[j] = buffer[i][j - target].i;
      OUTx[j] = buffer[i][j - target].x;
    }
  }
}

/* prepare the ctx object with information needed for all passes */
void planification(struct ctx_t *ctx, mtx_CSR *R, mtx_entry *scratch);

/* returns k such that buckets [0:k] are non-empty. */
u32 partitioning(struct ctx_t *ctx, const mtx_COO *A, mtx_CSR *R,
                 cacheline *buffer, u32 *tCOUNT, u32 *gCOUNT);

void histogram(const struct ctx_t *ctx, const mtx_entry *Te, const u32 lo,
               const u32 hi, const u8 n, u32 **W);

/* sequentially transpose a single bucket */
void transpose_bucket(struct ctx_t *ctx, cacheline *buffer, const u32 lo,
                      const u32 hi, mtx_CSR *R, const u32 bucket);

/* converts a sparse matrix in COOrdinate format to the CSR format.
   INPUT:  COO sparse matrix in Ai, Aj (both of size nnz), with n rows

   OUTPUT: CSR sparse matrix in Rp, Rj.

   Rp and Rj MUST be preallocated (of sizes n+1 and nnz, respectively).
   The "row pointers" Rp MUST be already computed.

   Ai, Aj, Rj MUST be aligned on a 64-byte boundary (for good cache behavior).
   The input arrays are expendable (i.e. they might be destroyed).
   The current code only reads them though. */
void transpose(const mtx_COO *A, mtx_CSR *R, const u32 num_threads);

#endif /* INCLUDE_DRIVER_BB3_TRANSPOSE3_H */
