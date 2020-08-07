
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

#include "mini_spasm.h"
#include "mmio.h"
#include "transpose.h"
#include "tools.h"

void run_test(const char *matrix_filename, const char *output_filename)
{
  FILE *f = fopen(matrix_filename, "r");
  if (f == NULL)
    err(1, "impossible to open %s", matrix_filename);

  // Loading matrix
  spasm_triplet *T = spasm_load_mm(f);
  fclose(f);

  // Transposing
  spasm_triplet *R = malloc(sizeof(spasm_triplet));
  spasm_triplet_transpose(T, R);

  spasm *A = spasm_csr_alloc(T->m, T->n, T->nnz); // T transposed, in CSR

  transpose(T, A);
  check(R, A);

  spasm_triplet_free(T);
  free(R);
  spasm_csr_free(A);
}

int main(void)
{
// 	#ifdef BENCHMARK_LARGE_MATRICES
//   for (u32 i = 1; i <= N_LARGE_MATRICES; i++)
//   { // just this one is enough to exhibit the crash
//     char matrix_filename[FILENAME_MAX];
//     sprintf(matrix_filename, "%s/RSA.ok/pre_transpose%d.mtx", MATRIX_PATH, i);

//     printf("#---------------------------------------- %s\n", matrix_filename);
//     run_test(matrix_filename, NULL);
//     fprintf(stderr, "\n");
//   }
// #endif // BENCHMARK_LARGE_MATRICES

  run_test("../matrices/language.mtx", NULL);

  return 0;
}
