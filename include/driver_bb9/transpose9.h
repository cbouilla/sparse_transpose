///
/// \file transpose.h
/// \author Charles Bouillaguet (Github: cbouilla) and Jérôme Bonacchi (Github:
/// MarsParallax)
/// \brief This file contains functions to transpose matrices in CSR format
/// using radix sort.
/// \date 2020
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_DRIVER_BB8_TRANSPOSE8_H
#define INCLUDE_DRIVER_BB8_TRANSPOSE8_H

#include <assert.h>
#include <err.h>
#include <stdint.h>
#include <stdlib.h>

#include "sparse.h"
#include "tools.h"

///
/// \brief Number of variables in L1 cache line. L1 cache line has size 64 on
/// most CPUs.
///
#define CACHELINE_SIZE ((u8)(64 / sizeof(mtx_entry)))
#define L1_CACHE_SIZE ((u8)9)  // ((u32)((1 << 15) / (8 * sizeof(cacheline))))
#define L2_CACHE_SIZE ((u8)14) // ((u32)((1 << 20) / (8 * sizeof(cacheline))))
#define L3_CACHE_SIZE ((u8)19) // ((u32)((1 << 24) / (8 * sizeof(cacheline))))
#define MAX_RADIX_BITS 8
#define MAX_PASSES ((u8)(8 * sizeof(u32) / sizeof(MAX_RADIX_BITS)) + 1) // 4

///
/// \brief A structure to store the information for the radix sort
///
struct ctx_t
{
  u8 bits;                   ///< length in bits of the numbers
  u8 n_passes;               ///< number of passes (first in parallel)
  u8 radix[MAX_PASSES];      ///< radices used in each pass
  u8 shift[MAX_PASSES];      ///< shifts used in each pass
  u32 n_buckets[MAX_PASSES]; ///< number of buckets used in each pass
  u32 pCOUNT[MAX_PASSES];
  u32 mask[MAX_PASSES];       ///< masks used in each pass
  mtx_entry *OUT[MAX_PASSES]; ///< triplets
  u32 seq_count_size;
  u32 par_count_size;
};

///
/// \brief cache-resident buffer for triplets. One such entry per output bucket.
///
typedef mtx_entry cacheline[CACHELINE_SIZE];

///
/// \brief Copies 64 bytes from `src` to `dst` using non-temporal store
/// instructions if available (this bypasses the cache).
///
/// \param[out] dst the destination array
/// \param[in] src the source array
///
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
  // note : it can also be done using SSE for non-AVX machines
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

#endif // __AVX__

///
/// \brief Allocates aligned memory. Wrapper of `aligned_alloc`.
///
/// \param[in] size
/// \param[in] alignment
/// \return void* the pointer to allocated memory
///
static inline void *malloc_aligned(const size_t size, const size_t alignment)
{
  void *x = aligned_alloc(alignment, size);
  if (x == NULL)
    err(1, "malloc failed");
  return x;
}

///
/// \brief Allocates the software write-combining buffer.
///
/// \return struct cacheline_t* a list of `cacheline_t`, i.e. buckets
///
static cacheline *wc_alloc()
{
  return (cacheline *)malloc_aligned(sizeof(cacheline) * (1 << MAX_RADIX_BITS),
                                     64);
}

///
/// \brief Sets up the buffer for a new pass.
///
/// \param[out] buffer the buffer, i.e. a list of buckets
/// \param[in] COUNT
/// \param[in] n_buckets the number of buckets
///
static inline void wc_prime(cacheline *buffer, const u32 *COUNT,
                            const u32 n_buckets)
{
  for (u32 i = 0; i < n_buckets; i++)
  {
    buffer[i][CACHELINE_SIZE - 1].i = COUNT[i];
    buffer[i][CACHELINE_SIZE - 1].j = COUNT[i] & (CACHELINE_SIZE - 1);
  }
}

///
/// \brief Transfers data stored in the bucket to the output arrays.
/// Assumption: this bucket is filled to the end
///
/// \param[in] self the bucket
/// \param[in] count
/// \param[in] start the index where to start in the buffer
/// \param[out] OUT the triplets
///
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

///
/// \brief Pushes an (i, j, x) triplet into the buffer.
///
/// \param[in] entry the triplet
/// \param[inout] buffer the buffer
/// \param[in] bucket_idx the index of the bucket in the buffer
/// \param[out] OUT the triplets
///
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

///
/// \brief Flushes all buffer entries to the OUT arrays.
///
/// \param[in] buffer the buffer
/// \param[in] n_buckets the number of buckets
/// \param[out] OUT the triplets
///
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

///
/// \brief Prepares the ctx object with information needed for all passes.
///
/// \param[out] ctx
/// \param[in] R the output matrix in CSR format
/// \param[in] scratch a scratch space for triplets
///
void planification(struct ctx_t *ctx, mtx_CSR *R, mtx_entry *scratch);

///
/// \brief Does the first pass of the radix sort. Returns k such that buckets
/// [0:k] are non-empty.
///
/// \param[inout] ctx
/// \param[in] A the input matrix
/// \param[inout] buffer the buffer
/// \param[out] tCOUNT
/// \param[out] gCOUNT
/// \return u32 k such that buckets [0:k] are non-empty
///
u32 partitioning(struct ctx_t *ctx, const mtx_COO *A, cacheline *buffer,
                 u32 *tCOUNT, u32 *gCOUNT);

///
/// \brief Computes the histogram for all the passes of the radix sort for the
/// indices between `lo` and `hi`.
///
/// \param[in] ctx
/// \param[in] Aj the input triplets
/// \param[in] lo the lower bound of column indices
/// \param[in] hi the upper bound of column indices
/// \param[in] n the number of passes
/// \param[out] W the histogram
///
void histogram(const struct ctx_t *ctx, const mtx_entry *Te, const u32 lo,
               const u32 hi, const u8 n, u32 **W);

///
/// \brief Sequentially transposes a single bucket.
///
/// \param[inout] ctx
/// \param[in] buffer the bucket
/// \param[in] lo the lower bound of column indices
/// \param[in] hi the upper bound of column indices
/// \param[out] R the output matrix in CSR format
/// \param[in] the bucket index in the buffer
///
void transpose_bucket(struct ctx_t *ctx, cacheline *buffer, const u32 lo,
                      const u32 hi, mtx_CSR *R, const u32 bucket);

///
/// \brief Converts a sparse matrix in COOrdinate format to the CSR format.
/// Rp and Rj MUST be preallocated (of sizes n+1 and nnz, respectively).
/// The "row pointers" Rp MUST be already computed.
/// Ai, Aj, Rj MUST be aligned on a 64-byte boundary (for good cache behavior).
/// The input arrays are expendable (i.e. they might be destroyed).
/// The current code only reads them though.
///
/// \param[in] A the input matrix
/// \param[out] R the output matrix
/// \param[in] num_threads the number of threads to use
///
void transpose(const mtx_COO *A, mtx_CSR *R, const u32 num_threads);

#endif /* INCLUDE_DRIVER_BB8_TRANSPOSE8_H */
