#ifndef INCLUDE_DRIVER_BB_TRANSPOSE_H
#define INCLUDE_DRIVER_BB_TRANSPOSE_H

#include <assert.h>
#include <err.h>
#include <stdint.h>
#include <stdlib.h>

#include "mini_spasm.h"

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
  double *OUTx[MAX_PASSES];
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
  double value[CACHELINE_SIZE];
};

struct half_cacheline_t
{
  index_t row[CACHELINE_SIZE];
  double value[CACHELINE_SIZE];
};

/* copy 64 bytes from src to dst using non-temporal store instructions
   if available (this bypasses the cache). */
static inline void store_nontemp_int(void *dst, void *src);

/* copy 128 bytes from src to dst using non-temporal store instructions
   if available (this bypasses the cache). */
static inline void store_nontemp_double(void *dst, void *src);

/* converts a sparse matrix in COOrdinate format to the CSR format.
   INPUT:  COO sparse matrix in Ai, Aj (both of size nnz), with n rows

   OUTPUT: CSR sparse matrix in Rp, Ri.

   Rp and Ri MUST be preallocated (of sizes n+1 and nnz, respectively).
   The "row pointers" Rp MUST be already computed.

   Ai, Aj, Ri MUST be aligned on a 64-byte boundary (for good cache behavior).
   The input arrays are expendable (i.e. they might be destroyed).
   The current code only reads them though. */
void transpose(spasm_triplet *A, spasm *R);

#if __AVX__
#include <immintrin.h>
static inline void store_nontemp_int(void *dst, void *src)
{
  register __m256i *d1 = (__m256i *)dst;
  register __m256i s1 = *((__m256i *)src);
  register __m256i *d2 = d1 + 1;
  register __m256i s2 = *(((__m256i *)src) + 1);
  _mm256_stream_si256(d1, s1);
  _mm256_stream_si256(d2, s2);
  /* note : it can also be done using SSE for non-AVX machines */
}

static inline void store_nontemp_double(void *dst, void *src)
{
  register __m256d *d1 = (__m256d *)dst;
  register __m256d s1 = *((__m256d *)src);
  register __m256d *d2 = d1 + 1;
  register __m256d s2 = *(((__m256d *)src) + 1);
  register __m256d *d3 = d1 + 2;
  register __m256d s3 = *(((__m256d *)src) + 2);
  register __m256d *d4 = d1 + 3;
  register __m256d s4 = *(((__m256d *)src) + 3);
  _mm256_stream_sd(d1, s1);
  _mm256_stream_sd(d2, s2);
  _mm256_stream_sd(d3, s3);
  _mm256_stream_sd(d4, s4);
  /* note : it can also be done using SSE for non-AVX machines */
}
#else
static inline void store_nontemp_int(void *dst, void *src)
{
  index_t *in = src;
  index_t *out = dst;
  for (int i = 0; i < CACHELINE_SIZE; i++)
    out[i] = in[i];
}

static inline void store_nontemp_double(void *dst, void *src)
{
  double *in = src;
  double *out = dst;
  for (int i = 0; i < CACHELINE_SIZE; i++)
    out[i] = in[i];
}
#endif

static inline void *malloc_aligned(int size, int alignment)
{
  void *x = aligned_alloc(alignment, size);
  if (x == NULL)
    err(1, "malloc failed");
  return x;
}

static inline void free_aligned(void *ptr) { free(ptr); }

#endif /* INCLUDE_DRIVER_BB_TRANSPOSE_H */
