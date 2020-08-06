#ifndef INCLUDE_DRIVER_BB_TRANSPOSE_H
#define INCLUDE_DRIVER_BB_TRANSPOSE_H

#include <assert.h>
#include <err.h>
#include <stdint.h>
#include <stdlib.h>


// #include "typedefs.h"

// #define TRANSPOSE_EASY_WAY

/* L1 cache line has size 64 on most CPUs */
#define CACHELINE_SIZE ((int)(64 / sizeof(index_t)))
#define MAX_RADIX_BITS 10 /* was experimentally found to be OK */
#define MAX_PASSES 4

typedef uint32_t index_t;

struct ctx_t
{
  int bits;
  int n_passes; /* last pass is done... first, in parallel. */
  int radix[MAX_PASSES];
  index_t shift[MAX_PASSES];
  int n_buckets[MAX_PASSES];
  int pCOUNT[MAX_PASSES];
  index_t mask[MAX_PASSES];
  index_t *OUTi[MAX_PASSES];
  index_t *OUTj[MAX_PASSES];
  int seq_count_size;
  int par_count_size;
};

/* cache-resident buffer for (i, j) pairs. One such entry per output bucket.
   Invariants:  row[CACHELINE_SIZE - 1] contains COUNT[...] for this bucket,
        col[CACHELINE_SIZE - 1] contains offset of the first entry in this
   buffer.
   */
struct cacheline_t
{
  index_t row[CACHELINE_SIZE];
  index_t col[CACHELINE_SIZE];
};

struct half_cacheline_t
{
  index_t row[CACHELINE_SIZE];
};

/* copy 64 bytes from src to dst using non-temporal store instructions
   if available (this bypasses the cache). */
static inline void store_nontemp_64B(void *dst, void *src);

/* converts a sparse matrix in COOrdinate format to the CSR format.
   INPUT:  COO sparse matrix in Ai, Aj (both of size nnz), with n rows

   OUTPUT: CSR sparse matrix in Rp, Ri.

   Rp and Ri MUST be preallocated (of sizes n+1 and nnz, respectively).
   The "row pointers" Rp MUST be already computed.

   Ai, Aj, Ri MUST be aligned on a 64-byte boundary (for good cache behavior).
   The input arrays are expendable (i.e. they might be destroyed).
   The current code only reads them though. */
void transpose(uint64_t nnz, index_t *Ai, index_t *Aj, double *Ax, index_t n, index_t *Rp,
               index_t *Ri, double *Rx);

#if __AVX__
#include <immintrin.h>
static inline void store_nontemp_64B(void *dst, void *src)
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
static inline void store_nontemp_64B(void *dst, void *src)
{
  index_t *in = src;
  index_t *out = dst;
  for (int i = 0; i < CACHELINE_SIZE; i++)
    out[i] = in[i];
}
#endif

static inline void* malloc_aligned(int size, int alignment) {
  void *x = aligned_alloc(alignment, size);
  if (x == NULL)
    err(1, "malloc failed");
  return x;
}

static inline void free_aligned(void* ptr)
{
  free(ptr);
}

#endif /* INCLUDE_DRIVER_BB_TRANSPOSE_H */
