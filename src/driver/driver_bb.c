
///
/// \file driver_bb.c
/// \author Charles Bouillaguet and Jérôme Bonacchi
/// \brief This files runs benchmarks.
/// \date 2020-08-06
///
/// @copyright Copyright (c) 2020
///
///

#include <err.h>
#include <stdio.h>

#include "classical_sort.h"
#include "mini_spasm.h"
#include "mmio.h"
#include "transpose.h"

void run_test(const char *matrix_filename, const char *output_filename)
{
  FILE *f = fopen(matrix_filename, "r");
  if (f == NULL)
    err(1, "impossible to open %s", matrix_filename);

  // Loading matrix
  spasm_triplet *T = spasm_load_mm(f);
  fclose(f);

  spasm *A = spasm_csr_alloc(T->n, T->m, T->nnz);
  spasm *B = spasm_csr_alloc(T->m, T->n, T->nnz);
  spasm *C = spasm_csr_alloc(T->m, T->n, T->nnz);
  C->x = T->x; // TODO
  u32 *W = (u32 *)spasm_malloc((spasm_max(T->m, T->n) + 1) * sizeof(*W));
  classical_compress(T, A, W);
  classical_transpose(A, B, W);

  transpose(T->nnz, T->i, T->j, T->x, C->n, C->p, C->j, C->x);

  for (u32 i = 0; i < C->n; i++)
  {
    for (u32 j = C->p[i]; j < C->p[i + 1]; j++)
    {
      assert(B->j[j] == C->j[j]);
      // printf("(%d, %d) %f\t| (%d, %d) %f\n", i, B->j[j], B->x[j], i, C->j[j], C->x[j]);
    }
  }

  spasm_triplet_free(T);
  spasm_csr_free(A);
  spasm_csr_free(B);
  spasm_csr_free(C);
  free(W);
  // Transposing
  // spasm_triplet *R = malloc(sizeof(spasm_triplet));
  // spasm_triplet_transpose(T, R);

  // spasm_triplet_free(T);
  // free(R);
}

int main(void)
{
  run_test("../matrices/pre_transpose12.mtx", NULL);
  // spasm_triplet *T = spasm_triplet_alloc(8);
  // T->nnz = 8;
  // T->n = 3;
  // T->m = 5;

  // index_t tmp[] = {0, 0, 0, 1, 1, 2, 2, 2};
  // T->i = tmp;
  // index_t tmp2[] = {0, 1, 2, 2, 3, 0, 1, 4};
  // T->j = tmp2;
  // double tmp3[] = {0, 1, 2, 3, 4, 5, 6, 7};
  // T->x = tmp3;
  // spasm_triplet *R = malloc(sizeof(spasm_triplet));
  // spasm_triplet_transpose(T, R);

  return 0;
}
