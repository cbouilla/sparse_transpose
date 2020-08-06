///
/// \file driver.c
/// \author Charles Bouillaguet and Jérôme Bonacchi
/// \brief This files runs benchmarks.
/// \date 2020-07-09
///
/// @copyright Copyright (c) 2020
///

#define _XOPEN_SOURCE

#include <assert.h>
#include <err.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_MKL
#include <mkl.h>
#include <mkl_spblas.h>
#endif // HAVE_MKL

#include "classical_sort.h"
#include "mini_spasm.h"
#include "simple_sort.h"
#include "tools.h"

///
/// \brief Benchmarks the "classical" algorithm.
///
/// \param[in] T the matrix in triplet format
/// \param[in] R the transposed matrix in triplet format
/// \param[in, out] duration the duration of the algorithms
///
void run_test_classical(const spasm_triplet *T, const spasm_triplet *R,
                        algorithm_times *duration)
{
  double start, stop;
  const u32 n = T->n;
  const u32 m = T->m;
  const u32 nnz = T->nnz;

  spasm *A = spasm_csr_alloc(n, m, nnz);
  spasm *B = spasm_csr_alloc(m, n, nnz);
  spasm *C = spasm_csr_alloc(m, n, nnz);
  spasm *D = spasm_csr_alloc(n, m, nnz);
  u32 *W = (u32 *)spasm_malloc((spasm_max(n, m) + 1) * sizeof(*W));

  start = spasm_wtime();
  classical_compress(T, A, W);
  stop = spasm_wtime();
  check(T, A);
  duration->compress = stop - start;
  total[GUSTAVSON].compress += duration->compress;
  fprintf(stderr, "-- Gustavson compress [COO->CSR]: %.3fs\n",
          duration->compress);

  start = spasm_wtime();
  classical_transpose(A, B, W);
  stop = spasm_wtime();
  check(R, B);
  duration->transpose = stop - start;
  total[GUSTAVSON].transpose += duration->transpose;
  fprintf(stderr, "-- Gustavson transpose [CSR->CSR']: %.3fs\n",
          duration->transpose);

  start = spasm_wtime();
  classical_compress(R, C, W);
  stop = spasm_wtime();
  check(R, C);
  duration->compress_tr = stop - start;
  total[GUSTAVSON].compress_tr += duration->compress_tr;
  fprintf(stderr, "-- Gustavson compress [COO'->CSR']: %.3fs\n",
          duration->compress_tr);

  start = spasm_wtime();
  classical_transpose(C, D, W);
  stop = spasm_wtime();
  check(T, D);
  duration->transpose_tr = stop - start;
  total[GUSTAVSON].transpose_tr += duration->transpose_tr;
  fprintf(stderr, "-- Gustavson transpose [CSR'->CSR]: %.3fs\n",
          duration->transpose_tr);

  spasm_csr_free(A);
  spasm_csr_free(B);
  spasm_csr_free(C);
  spasm_csr_free(D);
  free(W);
}

///
/// \brief Benchmarks the std::sort algorithm.
///
/// \param[in] T the matrix in triplet format
/// \param[in] R the transposed matrix in triplet format
/// \param[in, out] duration the duration of the algorithms
///
void run_test_stdsort(const spasm_triplet *T, const spasm_triplet *R,
                      algorithm_times *duration)
{
  double start, stop;
  const u32 n = T->n;
  const u32 m = T->m;
  const u32 nnz = T->nnz;

  spasm *A = spasm_csr_alloc(n, m, nnz);
  spasm *B = spasm_csr_alloc(m, n, nnz);
  spasm *C = spasm_csr_alloc(m, n, nnz);
  spasm *D = spasm_csr_alloc(n, m, nnz);
  struct matrix_entry_t *Te =
      (struct matrix_entry_t *)spasm_malloc(nnz * sizeof(*Te));

  start = spasm_wtime();
  stdsort_compress(T, A, Te);
  stop = spasm_wtime();
  check(T, A);
  duration->compress = stop - start;
  total[STDSORT].compress += duration->compress;
  fprintf(stderr, "-- std::sort compress [COO->CSR]: %.3fs\n",
          duration->compress);

  start = spasm_wtime();
  stdsort_transpose(A, B, Te);
  stop = spasm_wtime();
  check(R, B);
  duration->transpose = stop - start;
  total[STDSORT].transpose += duration->transpose;
  fprintf(stderr, "-- std::sort transpose [CSR->CSR']: %.3fs\n",
          duration->transpose);

  start = spasm_wtime();
  stdsort_compress(R, C, Te);
  stop = spasm_wtime();
  check(R, C);
  duration->compress_tr = stop - start;
  total[STDSORT].compress_tr += duration->compress_tr;
  fprintf(stderr, "-- std::sort compress [COO'->CSR']: %.3fs\n",
          duration->compress_tr);

  start = spasm_wtime();
  stdsort_transpose(C, D, Te);
  stop = spasm_wtime();
  check(T, D);
  duration->transpose_tr = stop - start;
  total[STDSORT].transpose_tr += duration->transpose_tr;
  fprintf(stderr, "-- std::sort transpose [CSR'->CSR]: %.3fs\n",
          duration->transpose_tr);

  spasm_csr_free(A);
  spasm_csr_free(B);
  spasm_csr_free(C);
  spasm_csr_free(D);
  free(Te);
}

#ifdef HAVE_TBB

///
/// \brief Benchmarks the tbb::parallel_sort algorithm.
///
/// \param[in] T the matrix in triplet format
/// \param[in] R the transposed matrix in triplet format
/// \param[in, out] duration the duration of the algorithms
///
void run_test_tbbsort(const spasm_triplet *T, const spasm_triplet *R,
                      algorithm_times *duration, const u32 num_threads)
{
  double start, stop;
  const u32 n = T->n;
  const u32 m = T->m;
  const u32 nnz = T->nnz;

  spasm *A = spasm_csr_alloc(n, m, nnz);
  spasm *B = spasm_csr_alloc(m, n, nnz);
  spasm *C = spasm_csr_alloc(m, n, nnz);
  spasm *D = spasm_csr_alloc(n, m, nnz);
  struct matrix_entry_t *Te =
      (struct matrix_entry_t *)spasm_malloc(nnz * sizeof(*Te));

  start = spasm_wtime();
  tbbsort_compress(T, A, Te, num_threads);
  stop = spasm_wtime();
  check(T, A);
  duration->compress = stop - start;
  total[TBBSORT].compress += duration->compress;
  fprintf(stderr,
          "-- tbb::parallel_sort (%d threads) compress [COO->CSR]: %.3fs\n",
          num_threads, duration->compress);

  start = spasm_wtime();
  tbbsort_transpose(A, B, Te, num_threads);
  stop = spasm_wtime();
  check(R, B);
  duration->transpose = stop - start;
  total[TBBSORT].transpose += duration->transpose;
  fprintf(stderr,
          "-- tbb::parallel_sort (%d threads) transpose [CSR->CSR']: %.3fs\n",
          num_threads, duration->transpose);

  start = spasm_wtime();
  tbbsort_compress(R, C, Te, num_threads);
  stop = spasm_wtime();
  check(R, C);
  duration->compress_tr = stop - start;
  total[TBBSORT].compress_tr += duration->compress_tr;
  fprintf(stderr,
          "-- tbb::parallel_sort (%d threads) compress [COO'->CSR']: %.3fs\n",
          num_threads, duration->compress_tr);

  start = spasm_wtime();
  tbbsort_transpose(C, D, Te, num_threads);
  stop = spasm_wtime();
  check(T, D);
  duration->transpose_tr = stop - start;
  total[TBBSORT].transpose_tr += duration->transpose_tr;
  fprintf(stderr,
          "-- tbb::parallel_sort (%d threads) transpose [CSR'->CSR]: %.3fs\n",
          num_threads, duration->transpose_tr);

  spasm_csr_free(A);
  spasm_csr_free(B);
  spasm_csr_free(C);
  spasm_csr_free(D);
  free(Te);
}
#endif // HAVE_TBB

#ifdef HAVE_MKL

///
/// \brief Benchmarks the MKL sort algorithm.
///
/// \param[in] T the matrix in triplet format
/// \param[in] R the transposed matrix in triplet format
/// \param[in, out] duration the duration of the algorithms
///
void run_test_MKL(const spasm_triplet *T, algorithm_times *duration,
                  const u32 num_threads)
{
  double start, stop;
  sparse_matrix_t mkl_T, mkl_A, mkl_B, mkl_C, mkl_D;
  sparse_status_t status;

  assert(sizeof(MKL_INT) == sizeof(int));

#ifdef _OPENMP
  mkl_set_num_threads(num_threads);
#endif // _OPENMP

  status = mkl_sparse_d_create_coo(&mkl_T, SPARSE_INDEX_BASE_ZERO, T->n, T->m,
                                   T->nnz, (int *)T->i, (int *)T->j, T->x);
  if (status != SPARSE_STATUS_SUCCESS)
    err(1, "MKL create_coo (T) failed");

  // COO --> CSR
  start = spasm_wtime();
  status =
      mkl_sparse_convert_csr(mkl_T, SPARSE_OPERATION_NON_TRANSPOSE, &mkl_A);
  stop = spasm_wtime();
  if (status != SPARSE_STATUS_SUCCESS)
    err(1, "MKL convert_csr (coo->csr) failed");
  duration->compress = stop - start;
  total[MKL].compress += duration->compress;
  fprintf(stderr, "-- MKL compress (%s, %d threads) [COO->CSR]: %.3fs\n",
          HAVE_MKL, num_threads, duration->compress);

  start = spasm_wtime();
  status = mkl_sparse_convert_csr(mkl_A, SPARSE_OPERATION_TRANSPOSE, &mkl_B);
  stop = spasm_wtime();
  if (status != SPARSE_STATUS_SUCCESS)
    err(1, "MKL convert_csr (csr->csr) failed");
  duration->transpose = stop - start;
  total[MKL].transpose += duration->transpose;
  fprintf(stderr, "-- MKL transpose (%s, %d threads) [CSR->CSR']: %.3fs\n",
          HAVE_MKL, num_threads, duration->transpose);

  start = spasm_wtime();
  status = mkl_sparse_convert_csr(mkl_T, SPARSE_OPERATION_TRANSPOSE, &mkl_C);
  stop = spasm_wtime();
  if (status != SPARSE_STATUS_SUCCESS)
    err(1, "MKL convert_csr (coo->csr) failed");
  duration->compress_tr = stop - start;
  total[MKL].compress_tr += duration->compress_tr;
  fprintf(stderr, "-- MKL compress (%s, %d threads) [COO'->CSR']: %.3fs\n",
          HAVE_MKL, num_threads, duration->compress_tr);

  start = spasm_wtime();
  status = mkl_sparse_convert_csr(mkl_C, SPARSE_OPERATION_TRANSPOSE, &mkl_D);
  stop = spasm_wtime();
  if (status != SPARSE_STATUS_SUCCESS)
    err(1, "MKL convert_csr (csr->csr) failed");
  duration->transpose_tr = stop - start;
  total[MKL].transpose_tr += duration->transpose_tr;
  fprintf(stderr, "-- MKL transpose (%s, %d threads) [CSR'->CSR]: %.3fs\n",
          HAVE_MKL, num_threads, duration->transpose_tr);

  status = mkl_sparse_destroy(mkl_T);
  if (status != SPARSE_STATUS_SUCCESS)
    err(1, "MKL destroy (T) failed");
  status = mkl_sparse_destroy(mkl_A);
  if (status != SPARSE_STATUS_SUCCESS)
    err(1, "MKL destroy (A) failed");
  status = mkl_sparse_destroy(mkl_B);
  if (status != SPARSE_STATUS_SUCCESS)
    err(1, "MKL destroy (B) failed");
  status = mkl_sparse_destroy(mkl_C);
  if (status != SPARSE_STATUS_SUCCESS)
    err(1, "MKL destroy (C) failed");
  status = mkl_sparse_destroy(mkl_D);
  if (status != SPARSE_STATUS_SUCCESS)
    err(1, "MKL destroy (D) failed");
}

///
/// \brief Prints the version of the Intel's MKL, the MKL option used a compile
/// time and the number of threads available.
///
void mkl_version()
{
  const int len = 198;
  char buf[len];
  mkl_get_version_string(buf, len);
  printf("MKL version: %s\n", buf);
  const int max_num_threads = mkl_get_max_threads();
  printf("MKL: %s, maximum %d threads\n", HAVE_MKL, max_num_threads);
}
#endif // HAVE_MKL

///
/// \brief Writes the execution durations of the classical algorithms.
///
/// \param[in] output_filename the filename where to write
/// \param[in] matrix_filename the name of the matrix studied
/// \param[in] duration the durations to write
///
void write_test_classical(const char *output_filename,
                          const char *matrix_filename,
                          algorithm_times *duration)
{
  FILE *file = fopen(output_filename, "a");
  if (file == NULL)
    err(1, "impossible to open %s", output_filename);
  const u32 num_threads = 1;
  const char *name = "Gustavson";
  for (unsigned short i = 0; i < N_REPEAT; i++)
  {
    fprintf(file, OUTPUT_FORMAT, name, CFLAGS, CXXFLAGS, num_threads,
            matrix_filename, duration[i].compress, duration[i].transpose,
            duration[i].compress_tr, duration[i].transpose_tr, i, N_REPEAT);
  }
  fclose(file);
}

///
/// \brief Writes the execution durations of the std::sort algorithms.
///
/// \param[in] output_filename the filename where to write
/// \param[in] matrix_filename the name of the matrix studied
/// \param[in] duration the durations to write
///
void write_test_stdsort(const char *output_filename,
                        const char *matrix_filename, algorithm_times *duration)
{
  FILE *file = fopen(output_filename, "a");
  if (file == NULL)
    err(1, "impossible to open %s", output_filename);
  const u32 num_threads = 1;
  const char *name = "STL::sort";
  for (unsigned short i = 0; i < N_REPEAT; i++)
  {
    fprintf(file, OUTPUT_FORMAT, name, CFLAGS, CXXFLAGS, num_threads,
            matrix_filename, duration[i].compress, duration[i].transpose,
            duration[i].compress_tr, duration[i].transpose_tr, i, N_REPEAT);
  }
  fclose(file);
}

///
/// \brief Writes the execution durations of the tbb::parallel_sort algorithms.
///
/// \param[in] output_filename the filename where to write
/// \param[in] matrix_filename the name of the matrix studied
/// \param[in] duration the durations to write
/// \param[in] num_threads the number of threads used
///
void write_test_tbbsort(const char *output_filename,
                        const char *matrix_filename, algorithm_times *duration,
                        const u32 num_threads)
{
  FILE *file = fopen(output_filename, "a");
  if (file == NULL)
    err(1, "impossible to open %s", output_filename);
  const char *name = "TBB::sort";
  for (unsigned short i = 0; i < N_REPEAT; i++)
  {
    fprintf(file, OUTPUT_FORMAT, name, CFLAGS, CXXFLAGS, num_threads,
            matrix_filename, duration[i].compress, duration[i].transpose,
            duration[i].compress_tr, duration[i].transpose_tr, i, N_REPEAT);
  }
  fclose(file);
}

///
/// \brief Writes the execution durations of the Intel's MKL algorithms.
///
/// \param[in] output_filename the filename where to write
/// \param[in] matrix_filename the name of the matrix studied
/// \param[in] duration the durations to write
/// \param[in] num_threads the number of threads used
///
void write_test_MKL(const char *output_filename, const char *matrix_filename,
                    algorithm_times *duration, const u32 num_threads)
{
  FILE *file = fopen(output_filename, "a");
  if (file == NULL)
    err(1, "impossible to open %s", output_filename);
  char name[15] = "MKL";
#ifdef HAVE_MKL
  sprintf(name, "%s %s", name, HAVE_MKL);
#endif // HAVE_MKL
  for (unsigned short i = 0; i < N_REPEAT; i++)
  {
    fprintf(file, OUTPUT_FORMAT, name, CFLAGS, CXXFLAGS, num_threads,
            matrix_filename, duration[i].compress, duration[i].transpose,
            duration[i].compress_tr, duration[i].transpose_tr, i, N_REPEAT);
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
  FILE *f = fopen(matrix_filename, "r");
  if (f == NULL)
    err(1, "impossible to open %s", matrix_filename);

    // Even with TBB, OpenMP is needed to manage the number of threads used
#if defined _OPENMP && (defined HAVE_TBB || defined HAVE_MKL)
  u32 max_num_threads = 1;
#pragma omp parallel
  max_num_threads = omp_get_num_threads();
#endif // _OPENMP

  algorithm_times duration[N_REPEAT];
  for (u32 i = 0; i < N_REPEAT; i++)
  {
    clear_times(&duration[i]);
  }

  // Loading matrix
  spasm_triplet *T = spasm_load_mm(f);
  fclose(f);

  // Transposing
  spasm_triplet *R = malloc(sizeof(spasm_triplet));
  spasm_triplet_transpose(T, R);

  // Running classical
  for (u32 i = 0; i < N_REPEAT; ++i)
  {
    fprintf(stderr, "-- Step %d/%d:\n", i + 1, N_REPEAT);
    run_test_classical(T, R, &duration[i]);
  }
  write_test_classical(output_filename, matrix_filename, duration);
  for (u32 i = 0; i < N_REPEAT; i++)
  {
    clear_times(&duration[i]);
  }

  // Running std::sort
  for (u32 i = 0; i < N_REPEAT; ++i)
  {
    fprintf(stderr, "-- Step %d/%d:\n", i + 1, N_REPEAT);
    run_test_stdsort(T, R, &duration[i]);
  }
  write_test_stdsort(output_filename, matrix_filename, duration);
  for (u32 i = 0; i < N_REPEAT; i++)
  {
    clear_times(&duration[i]);
  }

  // Running TBB
#ifdef HAVE_TBB

  for (u32 i_thread = 1; i_thread <= max_num_threads; i_thread++)
  {
    for (u32 i = 0; i < N_REPEAT; ++i)
    {
      fprintf(stderr, "-- Step %d/%d:\n", i + 1, N_REPEAT);
      run_test_tbbsort(T, R, &duration[i], i_thread);
    }
    write_test_tbbsort(output_filename, matrix_filename, duration, i_thread);
  }
  for (u32 i = 0; i < N_REPEAT; i++)
  {
    clear_times(&duration[i]);
  }

#endif // HAVE_TBB

  // Running MKL
#ifdef HAVE_MKL

#ifdef _OPENMP
  for (u32 i_thread = 1; i_thread <= max_num_threads; i_thread++)
  {
    for (u32 i = 0; i < N_REPEAT; ++i)
    {
      fprintf(stderr, "-- Step %d/%d:\n", i + 1, N_REPEAT);
      run_test_MKL(T, &duration[i], i_thread);
    }
    write_test_MKL(output_filename, matrix_filename, duration, i_thread);
  }
#else
  for (u32 i = 0; i < N_REPEAT; ++i)
  {
    fprintf(stderr, "-- Step %d/%d:\n", i + 1, N_REPEAT);
    run_test_MKL(T, &duration[i], 1);
  }
  write_test_MKL(output_filename, matrix_filename, duration, 1);
#endif // _OPENMP

  for (u32 i = 0; i < N_REPEAT; i++)
  {
    clear_times(&duration[i]);
  }

#endif // HAVE_MKL

  fflush(stdout);
  spasm_triplet_free(T);
  free(R);
}

///
/// \brief Prints the values of the global array `total`.
///
void show_grand_totals(void)
{
  fprintf(stderr, "\nGRAND TOTALS:\n");
  fprintf(stderr, "  Gustavson:\n");
  fprintf(stderr, "    compress:      %.3fs\n", total[GUSTAVSON].compress);
  fprintf(stderr, "    compress_tr:   %.3fs\n", total[GUSTAVSON].compress_tr);
  fprintf(stderr, "    transpsose:    %.3fs\n", total[GUSTAVSON].transpose);
  fprintf(stderr, "    transpsose_tr: %.3fs\n", total[GUSTAVSON].transpose_tr);
  fprintf(stderr, "  std::sort:\n");
  fprintf(stderr, "    compress:      %.3fs\n", total[STDSORT].compress);
  fprintf(stderr, "    compress_tr:   %.3fs\n", total[STDSORT].compress_tr);
  fprintf(stderr, "    transpsose:    %.3fs\n", total[STDSORT].transpose);
  fprintf(stderr, "    transpsose_tr: %.3fs\n", total[STDSORT].transpose_tr);
#ifdef HAVE_TBB
  fprintf(stderr, "  tbb::parallel_sort:\n");
  fprintf(stderr, "    compress:      %.3fs\n", total[TBBSORT].compress);
  fprintf(stderr, "    compress_tr:   %.3fs\n", total[TBBSORT].compress_tr);
  fprintf(stderr, "    transpsose:    %.3fs\n", total[TBBSORT].transpose);
  fprintf(stderr, "    transpsose_tr: %.3fs\n", total[TBBSORT].transpose_tr);
#endif // HAVE_TBB

#ifdef HAVE_MKL
  fprintf(stderr, "  MKL:\n");
  fprintf(stderr, "    compress:      %.3fs\n", total[MKL].compress);
  fprintf(stderr, "    compress_tr:   %.3fs\n", total[MKL].compress_tr);
  fprintf(stderr, "    transpsose:    %.3fs\n", total[MKL].transpose);
  fprintf(stderr, "    transpsose_tr: %.3fs\n", total[MKL].transpose_tr);
#endif // HAVE_MKL
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
    err(1, "Usage: ./driver [output_filename]\nThe default filename is '%s'\n",
        OUTPUT_FILENAME);
  }

  printf("Output file: %s\n", output_filename);

  for (int i = 0; i < N_METHOD; i++)
  {
    clear_times(&total[i]);
  }

#ifdef HAVE_TBB
  tbb_version();
#endif // HAVE_TBB
#ifdef HAVE_MKL
  mkl_version();
#endif // HAVE_MKL

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
#endif // BENCHMARK_SMAL_MATRICES

#ifdef BENCHMARK_LARGE_MATRICES
  for (u32 i = 1; i <= N_LARGE_MATRICES; i++)
  { // just this one is enough to exhibit the crash
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
