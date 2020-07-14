#include <cstdio>
#include <err.h>
#include <algorithm>
#include <stdlib.h>

#include <sys/time.h> // timing
#include "matio.h"
#include "sptrans.h"
#include "tools.h"

double spasm_wtime() {
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1E6;
}

#ifdef _OPENMP

void run_test_scanTrans(const char *filename, struct bench_time *duration, int num_threads)
{
	omp_set_num_threads(num_threads);
// input
    int m, n, nnzA;
    int *csrRowPtrA;
    int *csrColIdxA;
    double *csrValA;
    int retCode = read_mtx_mat(m, n, nnzA, csrRowPtrA, csrColIdxA, csrValA, filename);
    if(retCode != 0)
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
    std::fill_n(cscColPtrA, n+1, 0);

	start = spasm_wtime();
	sptrans_scanTrans<int, double>(m, n, nnzA, csrRowPtrA, csrColIdxA, csrValA, cscRowIdxA, cscColPtrA, cscValA);
	stop = spasm_wtime();
	duration->transpose = stop - start;
	total[SCANTRANS].transpose += duration->transpose;
	fprintf(stderr, "-- ScanTrans transpose (%d threads) [CSR->CSR']: %.3fs\n", num_threads, duration->transpose);

    // clear the buffers
    std::fill_n(csrColIdxA, nnzA, 0);
    std::fill_n(csrValA, nnzA, 0);
    std::fill_n(csrRowPtrA, n+1, 0);

	start = spasm_wtime();
	sptrans_scanTrans<int, double>(n, m, nnzA, cscColPtrA, cscRowIdxA, cscValA, csrColIdxA, csrRowPtrA, csrValA);
	stop = spasm_wtime();		
	duration->transpose_tr = stop - start;
	total[SCANTRANS].transpose_tr += duration->transpose_tr;
	fprintf(stderr, "-- ScanTrans transpose (%d threads) [CSR'->CSR]: %.3fs\n", num_threads, duration->transpose_tr);

    free(csrRowPtrA); 
    free(csrColIdxA); 
    free(csrValA);
    free(cscRowIdxA);
    free(cscColPtrA);
    free(cscValA);
}

void run_test_mergeTrans(const char* filename, struct bench_time *duration, int num_threads)
{
		omp_set_num_threads(num_threads);

// input
    int m, n, nnzA;
    int *csrRowPtrA;
    int *csrColIdxA;
    double *csrValA;
    int retCode = read_mtx_mat(m, n, nnzA, csrRowPtrA, csrColIdxA, csrValA, filename);
    if(retCode != 0)
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
    std::fill_n(cscColPtrA, n+1, 0);

	start = spasm_wtime();
    sptrans_mergeTrans<int, double>(m, n, nnzA, csrRowPtrA, csrColIdxA, csrValA, cscRowIdxA, cscColPtrA, cscValA);
	stop = spasm_wtime();
	duration->transpose = stop - start;
	total[MERGETRANS].transpose += duration->transpose;
	fprintf(stderr, "-- MergeTrans transpose (%d threads) [CSR->CSR']: %.3fs\n", num_threads, duration->transpose);

    // clear the buffers
    std::fill_n(csrColIdxA, nnzA, 0);
    std::fill_n(csrValA, nnzA, 0);
    std::fill_n(csrRowPtrA, n+1, 0);

	start = spasm_wtime();
	sptrans_mergeTrans<int, double>(n, m, nnzA, cscColPtrA, cscRowIdxA, cscValA, csrColIdxA, csrRowPtrA, csrValA);
	stop = spasm_wtime();		
	duration->transpose_tr = stop - start;
	total[MERGETRANS].transpose_tr += duration->transpose_tr;
	fprintf(stderr, "-- MergeTrans transpose (%d threads) [CSR'->CSR]: %.3fs\n", num_threads, duration->transpose_tr);

    free(csrRowPtrA); 
    free(csrColIdxA); 
    free(csrValA);
    free(cscRowIdxA);
    free(cscColPtrA);
    free(cscValA);
}

#endif // _OPENMP

void write_test_scanTrans(const char *output_filename, const char *matrix_filename, struct bench_time *duration, int num_thread)
{
	std::FILE* file = fopen(output_filename, "a");
	if (file == NULL)
		err(1, "impossible to open %s", output_filename);
	const char * name = "ScanTrans";
	for (size_t i = 0; i < N_REPEAT; i++)
	{
		fprintf(file, output_format, name, matrix_filename, num_thread, 0, duration[i].transpose, 0, duration[i].transpose_tr, i);
	}
	fclose(file);
}

void write_test_mergeTrans(const char *output_filename, const char *matrix_filename, struct bench_time *duration, int num_thread)
{
	std::FILE* file = fopen(output_filename, "a");
	if (file == NULL)
		err(1, "impossible to open %s", output_filename);
	const char * name = "MergeTrans";
	for (size_t i = 0; i < N_REPEAT; i++)
	{
		fprintf(file, output_format, name, matrix_filename, num_thread, 0, duration[i].transpose, 0, duration[i].transpose_tr, i);
	}
	fclose(file);
}

void run_test(const char *matrix_filename, const char * output_filename)
{
	struct bench_time duration[N_REPEAT];
	for (int i = 0; i < N_REPEAT; i++)
	{
		clear_bench_time(&duration[i]);
	}

#ifdef _OPENMP
	int max_num_threads;
	#pragma omp parallel
	max_num_threads = omp_get_num_threads();
	max_num_threads = 1;

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
		clear_bench_time(&duration[i]);
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
		clear_bench_time(&duration[i]);
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

int main()
{
	const char *output_filename = OUTPUT_FILENAME;
	//printf("Output file: %s\n", output_filename);

	for (int i = 0; i < N_METHOD; i++)
	{
	  clear_bench_time(&total[i]);
	}

#ifdef BENCHMARK_SMALL_MATRICES
	for (int i = 0; i < N; i++) {
	 	char matrix_filename[128];
	 	sprintf(matrix_filename, "%s/%s.mtx", MATRIX_PATH, matrices[i]);
		
		printf("#---------------------------------------- %s\n", matrix_filename);
		run_test(matrix_filename, output_filename);
		fprintf(stderr, "\n");
	}
	
	show_grand_totals();
#endif // BENCHMARK_SMAL_MATRICES

	return EXIT_SUCCESS;
}
