///
/// \file driver_wang.cpp
/// \author Charles Bouillaguet and Jérôme Bonacchi
/// \brief This files implements the ScanTrans and MergeTrans algorithms from
/// Wang et al. (2016).
/// The `run_test_scanTrans` and `run_test_mergeTrans` functions are heavily
/// based on their code (`main.cpp` ). \date 2020-07-23
///
/// @copyright Copyright (c) 2020
///

/*
 * (c) 2017 Virginia Polytechnic Institute & State University (Virginia Tech)
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License Version 2.1.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   LICENSE in the root of the repository for details.
 *
 */

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <err.h>
#include <stdlib.h>
#include <sys/time.h>

#include "matio.h"
#include "sptrans.h"
#include "tools.h"

double spasm_wtime()
{
  struct timeval ts;
  gettimeofday(&ts, NULL);
  return (double)ts.tv_sec + ts.tv_usec / 1E6;
}

#ifdef _OPENMP

void run_test_scanTrans(const char *filename, algorithm_times *duration,
                        int num_threads)
{
  omp_set_num_threads(num_threads);
  // input
  int m;
  int n;
  int nnzA;
  int *csrRowPtrA;
  int *csrColIdxA;
  double *csrValA;
  int retCode =
      read_mtx_mat(m, n, nnzA, csrRowPtrA, csrColIdxA, csrValA, filename);
  if (retCode != 0)
  {
    err(1, "Failed to read the matrix from %s.", filename);
  }

  double start, stop;

  int *cscRowIdxA = (int *)malloc(nnzA * sizeof(int));
  int *cscColPtrA = (int *)malloc((n + 1) * sizeof(int));
  double *cscValA = (double *)malloc(nnzA * sizeof(double));
  // clear the buffers
  std::fill_n(cscRowIdxA, nnzA, 0);
  std::fill_n(cscValA, nnzA, 0);
  std::fill_n(cscColPtrA, n + 1, 0);

  start = spasm_wtime();
  sptrans_scanTrans<int, double>(m, n, nnzA, csrRowPtrA, csrColIdxA, csrValA,
                                 cscRowIdxA, cscColPtrA, cscValA);
  stop = spasm_wtime();
  duration->transpose = stop - start;
  total[SCANTRANS].transpose += duration->transpose;
  fprintf(stderr, "-- ScanTrans transpose (%d threads) [CSR->CSR']: %.3fs\n",
          num_threads, duration->transpose);

  // clear the buffers
  std::fill_n(csrColIdxA, nnzA, 0);
  std::fill_n(csrValA, nnzA, 0);
  std::fill_n(csrRowPtrA, n + 1, 0);

  start = spasm_wtime();
  sptrans_scanTrans<int, double>(n, m, nnzA, cscColPtrA, cscRowIdxA, cscValA,
                                 csrColIdxA, csrRowPtrA, csrValA);
  stop = spasm_wtime();
  duration->transpose_tr = stop - start;
  total[SCANTRANS].transpose_tr += duration->transpose_tr;
  fprintf(stderr, "-- ScanTrans transpose (%d threads) [CSR'->CSR]: %.3fs\n",
          num_threads, duration->transpose_tr);

  free(csrRowPtrA);
  free(csrColIdxA);
  free(csrValA);
  free(cscRowIdxA);
  free(cscColPtrA);
  free(cscValA);
}

void run_test_mergeTrans(const char *filename, algorithm_times *duration,
                         int num_threads)
{
  omp_set_num_threads(num_threads);

  // input
  int m, n, nnzA;
  int *csrRowPtrA;
  int *csrColIdxA;
  double *csrValA;
  int retCode =
      read_mtx_mat(m, n, nnzA, csrRowPtrA, csrColIdxA, csrValA, filename);
  if (retCode != 0)
  {
    err(1, "Failed to read the matrix from %s.", filename);
  }

  double start, stop;

  int *cscRowIdxA = (int *)malloc(nnzA * sizeof(int));
  int *cscColPtrA = (int *)malloc((n + 1) * sizeof(int));
  double *cscValA = (double *)malloc(nnzA * sizeof(double));
  // clear the buffers
  std::fill_n(cscRowIdxA, nnzA, 0);
  std::fill_n(cscValA, nnzA, 0);
  std::fill_n(cscColPtrA, n + 1, 0);

  start = spasm_wtime();
  sptrans_mergeTrans<int, double>(m, n, nnzA, csrRowPtrA, csrColIdxA, csrValA,
                                  cscRowIdxA, cscColPtrA, cscValA);
  stop = spasm_wtime();
  duration->transpose = stop - start;
  total[MERGETRANS].transpose += duration->transpose;
  fprintf(stderr, "-- MergeTrans transpose (%d threads) [CSR->CSR']: %.3fs\n",
          num_threads, duration->transpose);

  // clear the buffers
  std::fill_n(csrColIdxA, nnzA, 0);
  std::fill_n(csrValA, nnzA, 0);
  std::fill_n(csrRowPtrA, n + 1, 0);

  start = spasm_wtime();
  sptrans_mergeTrans<int, double>(n, m, nnzA, cscColPtrA, cscRowIdxA, cscValA,
                                  csrColIdxA, csrRowPtrA, csrValA);
  stop = spasm_wtime();
  duration->transpose_tr = stop - start;
  total[MERGETRANS].transpose_tr += duration->transpose_tr;
  fprintf(stderr, "-- MergeTrans transpose (%d threads) [CSR'->CSR]: %.3fs\n",
          num_threads, duration->transpose_tr);

  free(csrRowPtrA);
  free(csrColIdxA);
  free(csrValA);
  free(cscRowIdxA);
  free(cscColPtrA);
  free(cscValA);
}

#endif // _OPENMP

void write_test_scanTrans(const char *output_filename,
                          const char *matrix_filename,
                          algorithm_times *duration, int num_threads)
{
  std::FILE *file = fopen(output_filename, "a");
  if (file == NULL)
    err(1, "impossible to open %s", output_filename);
  const char *name = "ScanTrans";
  for (unsigned short i = 0; i < N_REPEAT; i++)
  {
    fprintf(file, OUTPUT_FORMAT, name, CFLAGS, CXXFLAGS, num_threads,
            matrix_filename, 0.0, duration[i].transpose, 0.0,
            duration[i].transpose_tr, i, N_REPEAT);
  }
  fclose(file);
}

void write_test_mergeTrans(const char *output_filename,
                           const char *matrix_filename,
                           algorithm_times *duration, int num_threads)
{
  std::FILE *file = fopen(output_filename, "a");
  if (file == NULL)
    err(1, "impossible to open %s", output_filename);
  const char *name = "MergeTrans";
  for (unsigned short i = 0; i < N_REPEAT; i++)
  {
    fprintf(file, OUTPUT_FORMAT, name, CFLAGS, CXXFLAGS, num_threads,
            matrix_filename, 0.0, duration[i].transpose, 0.0,
            duration[i].transpose_tr, i, N_REPEAT);
  }
  fclose(file);
}

void run_test(const char *matrix_filename, const char *output_filename)
{
  algorithm_times duration[N_REPEAT];
  for (int i = 0; i < N_REPEAT; i++)
  {
    clear_times(&duration[i]);
  }

#ifdef _OPENMP
  int max_num_threads = 1;
#pragma omp parallel
  max_num_threads = omp_get_num_threads();

  for (int i_thread = 1; i_thread <= max_num_threads; i_thread++)
  {
    for (int i = 0; i < N_REPEAT; ++i)
    {
      fprintf(stderr, "-- Step %d/%d:\n", i + 1, N_REPEAT);
      run_test_scanTrans(matrix_filename, &duration[i], i_thread);
    }
    write_test_scanTrans(output_filename, matrix_filename, duration, i_thread);
  }
  for (int i = 0; i < N_REPEAT; i++)
  {
    clear_times(&duration[i]);
  }

  for (int i_thread = 1; i_thread <= max_num_threads; i_thread++)
  {
    for (int i = 0; i < N_REPEAT; ++i)
    {
      fprintf(stderr, "-- Step %d/%d:\n", i + 1, N_REPEAT);
      run_test_mergeTrans(matrix_filename, &duration[i], i_thread);
    }
    write_test_mergeTrans(output_filename, matrix_filename, duration, i_thread);
  }
  for (int i = 0; i < N_REPEAT; i++)
  {
    clear_times(&duration[i]);
  }

#endif // _OPENMP
}

void show_grand_totals()
{
#ifdef _OPENMP
  fprintf(stderr, "\nGRAND TOTALS:\n");
  fprintf(stderr, "  ScanTrans:\n");
  fprintf(stderr, "    transpsose:    %.3fs\n", total[SCANTRANS].transpose);
  fprintf(stderr, "    transpsose_tr: %.3fs\n", total[SCANTRANS].transpose_tr);

  fprintf(stderr, "  MergeTrans:\n");
  fprintf(stderr, "    transpsose:    %.3fs\n", total[MERGETRANS].transpose);
  fprintf(stderr, "    transpsose_tr: %.3fs\n", total[MERGETRANS].transpose_tr);
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
    fprintf(stderr,
            "Usage: ./driver [output_filename]\nThe default filename is '%s'",
            OUTPUT_FILENAME);
    return EXIT_FAILURE;
  }

  printf("Output file: %s\n", output_filename);

  for (int i = 0; i < N_METHOD; i++)
  {
    clear_times(&total[i]);
  }

#ifdef BENCHMARK_SMALL_MATRICES
  for (int i = 0; i < N_SMALL_MATRICES; i++)
  {
    char matrix_filename[FILENAME_MAX];
    sprintf(matrix_filename, "%s/%s.mtx", MATRIX_PATH, matrices[i]);

    printf("#---------------------------------------- %s\n", matrix_filename);
    run_test(matrix_filename, output_filename);
    fprintf(stderr, "\n");
  }

  show_grand_totals();
#endif // BENCHMARK_SMAL_MATRICES

#ifdef BENCHMARK_LARGE_MATRICES
  for (int i = 1; i <= N_LARGE_MATRICES; i++)
  {
    char matrix_filename[FILENAME_MAX];
    sprintf(matrix_filename, "%s/RSA.ok/pre_transpose%d.mtx", MATRIX_PATH, i);

    printf("#---------------------------------------- %s\n", matrix_filename);
    run_test(matrix_filename, output_filename);
    fprintf(stderr, "\n");
  }
  show_grand_totals();
#endif // BENCHMARK_LARGE_MATRICES

  return EXIT_SUCCESS;
}
