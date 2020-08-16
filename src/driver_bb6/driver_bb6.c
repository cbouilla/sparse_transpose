
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
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "sparse.h"
#include "mmio.h"
#include "tools.h"
#include "transpose5.h"

///
/// \brief Benchmarks the radix sort algorithm.
///
/// \param[in] T the matrix in COO format
/// \param[in] R the transposed matrix in COO format
/// \param[in, out] duration the duration of the algorithms
///
void run_test_radixsort1(const mtx_COO *T, const mtx_COO *R,
                         algorithm_times *duration, const u32 num_threads)
{
  double start, stop;
  const u32 n = T->n;
  const u32 m = T->m;
  const u32 nnz = T->nnz;

  mtx_CSR *A = mtx_CSR_alloc(m, n, nnz);
  mtx_CSR *B = mtx_CSR_alloc(n, m, nnz);

  start = spasm_wtime();
  transpose(T, A, num_threads);
  stop = spasm_wtime();
  check(R, A);
  duration->transpose = stop - start;
  total[RADIXSORT].transpose += duration->transpose;
  fprintf(stderr, "-- radix sort 6 (%d threads) transpose [COO->CSR']: %.3fs\n",
          num_threads, duration->transpose);

  start = spasm_wtime();
  transpose(R, B, num_threads);
  stop = spasm_wtime();
  check(T, B);
  duration->transpose_tr = stop - start;
  total[RADIXSORT].transpose_tr += duration->transpose_tr;
  fprintf(stderr, "-- radix sort 6 (%d threads) transpose [COO'->CSR]: %.3fs\n",
          num_threads, duration->transpose_tr);

  mtx_CSR_free(A);
  mtx_CSR_free(B);
}

///
/// \brief Writes the execution durations of the radix sort algorithms.
///
/// \param[in] output_filename the filename where to write
/// \param[in] matrix_filename the name of the matrix studied
/// \param[in] duration the durations to write
/// \param[in] num_threads the number of threads used
///
void write_test_radixsort1(const char *output_filename,
                           const char *matrix_filename,
                           algorithm_times *duration, const u32 num_threads)
{
  FILE *file = fopen(output_filename, "a");
  if (file == NULL)
    err(1, "impossible to open %s", output_filename);
  const char *name = "radix::sort::6";
  for (unsigned short i = 0; i < N_REPEAT; i++)
  {
    fprintf(file, OUTPUT_FORMAT, name, CFLAGS, CXXFLAGS, num_threads,
            matrix_filename, 0.0, duration[i].transpose, 0.0,
            duration[i].transpose_tr, i, N_REPEAT);
  }
  fclose(file);
}

///
/// \brief Runs benchmarks on several algorithms and writes the durations in a
/// file.
///
/// \param[in] matrix_filename the name of the matrix studied
/// \param[in] output_filename the filename where to write
///
void run_test(const char *matrix_filename, const char *output_filename)
{
#ifdef _OPENMP
  FILE *f = fopen(matrix_filename, "r");
  if (f == NULL)
    err(1, "impossible to open %s", matrix_filename);

  u32 max_num_threads = 1;
#pragma omp parallel
  max_num_threads = omp_get_num_threads();

  algorithm_times duration[N_REPEAT];
  for (u32 i = 0; i < N_REPEAT; i++)
  {
    clear_times(&duration[i]);
  }

  // Loading matrix
  mtx_COO *T = mtx_load_mm(f);
  fclose(f);

  // Transposing
  mtx_COO *R = malloc(sizeof(mtx_COO));
  mtx_COO_transpose(T, R);

  // Running radix sort 1
  for (u32 i_thread = 1; i_thread <= max_num_threads; i_thread++)
  {
    for (u32 i = 0; i < N_REPEAT; ++i)
    {
      fprintf(stderr, "-- Step %d/%d:\n", i + 1, N_REPEAT);
      run_test_radixsort1(T, R, &duration[i], i_thread);
    }
    write_test_radixsort1(output_filename, matrix_filename, duration, i_thread);
  }
  for (u32 i = 0; i < N_REPEAT; i++)
  {
    clear_times(&duration[i]);
  }

  fflush(stdout);
  mtx_COO_free(T);
  free(R);
#endif // _OPENMP
}

void show_grand_totals(void)
{
#ifdef _OPENMP
  fprintf(stderr, "\nGRAND TOTALS:\n");
  fprintf(stderr, "  Radix sort 6:\n");
  fprintf(stderr, "    transpsose:    %.3fs\n", total[RADIXSORT].transpose);
  fprintf(stderr, "    transpsose_tr: %.3fs\n", total[RADIXSORT].transpose_tr);
#endif // _OPENMP
}

int main(int argc, char **argv)
{
  char output_filename[FILENAME_MAX];
  if (argc == 1)
  {
    strcpy(output_filename, OUTPUT_FILENAME);
  }
  else if (argc == 2)
  {
    strcpy(output_filename, argv[1]);
  }
  else
  {
    err(1,
        "Usage: ./driver_bb [output_filename]\nThe default filename is '%s'\n",
        OUTPUT_FILENAME);
  }

  printf("Output file: %s\n", output_filename);

  for (int i = 0; i < N_METHOD; i++)
  {
    clear_times(&total[i]);
  }

#ifdef BENCHMARK_SMALL_MATRICES
  for (u32 i = 0; i < N_SMALL_MATRICES; i++)
  {
    char matrix_filename[FILENAME_MAX];
    sprintf(matrix_filename, "%s/%s.mtx", MATRIX_PATH, matrices[i]);

    printf("#---------------------------------------- %s\n", matrix_filename);
    run_test(matrix_filename, output_filename);
    fprintf(stderr, "\n");
  }
  show_grand_totals();
#endif // BENCHMARK_SMALL_MATRICES

#ifdef BENCHMARK_LARGE_MATRICES
  for (u32 i = 0; i < N_LARGE_MATRICES; i++)
  {
    char matrix_filename[FILENAME_MAX];
    sprintf(matrix_filename, "%s/RSA.ok/pre_transpose%d.mtx", MATRIX_PATH, pre_transpose[i]);

    printf("#---------------------------------------- %s\n", matrix_filename);
    run_test(matrix_filename, output_filename);
    fprintf(stderr, "\n");
  }
  show_grand_totals();
#endif // BENCHMARK_LARGE_MATRICES

  // run_test("../matrices/pre_transpose12.mtx", "tmp.csv");
  // run_test("../matrices/language.mtx", "tmp.csv");

  return 0;
}
