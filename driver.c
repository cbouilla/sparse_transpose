#define _XOPEN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#include <err.h>

#include "mini_spasm.h"
#include <omp.h>

void classical_compress(const spasm_triplet * T, spasm * A, int *W);
void classical_transpose(const spasm * A, spasm * R, int *W);

void stdsort_compress(const spasm_triplet * T, spasm * A, struct matrix_entry_t * Te);
void stdsort_transpose(const spasm * A, spasm * R, struct matrix_entry_t * Te);

void tbbsort_compress(const spasm_triplet * T, spasm * A, struct matrix_entry_t * Te, int num_threads);
void tbbsort_transpose(const spasm * A, spasm * R, struct matrix_entry_t * Te, int num_threads);
void tbb_version();

#ifdef HAVE_MKL
#include <mkl.h>
#include <mkl_spblas.h>
void mkl_version()
{
	const int len = 198;
	char buf[len];
	mkl_get_version_string(buf, len);
	printf("MKL version: %s\n", buf);
	int num_threads = mkl_get_max_threads();
	if (num_threads == 1)
	{
		printf("MKL: sequential\n");
	} else
	{
		printf("MKL: parallel, %d threads\n", num_threads);
	}
}
#endif

// #define BENCHMARK_SMALL_MATRICES
#define BENCHMARK_LARGE_MATRICES

/* offsets in the "total" array below */
#define GUSTAVSON 0      
#define MKL 1
#define STDSORT 2
#define TBBSORT 3
struct {
	double compress;           // time to "compress" the matrix    (convert   COO  -> CSR) 
	double compress_tr;        // time to "compress" the transpose (convert   COO' -> CSR') 
	double transpose;          // time to transpose the matrix     (transpose CSR  -> CSR')
	double transpose_tr;       // time to transpose the transpose  (transpose CSR' -> CSR)

} total[4];

const char * MATRIX_PATH="/Infos/lmd/2019/master/ue/MU4IN903-2020fev";

/* the matrices in "Parallel Transposition of Sparse Data Structures" by Wang, Liu, Hou and Feng */
#define N 21
const char *matrices[N] = {
	"language",
	"ASIC_680k",
	"circuit5M",
	"flickr",
	"memchip",
	"rajat21",
	"sme3Dc",
	"stomach",
	"transient",
	"webbase-1M",
	"wiki-Talk",
	"cage14",
	"eu-2005",
	"FullChip",
	"mac_econ_fwd500",
	"para-4",
	"rajat29",
	"Stanford_Berkeley",
	"torso1",
	"venkat01",
	"web-Google"
};

/* asserts that T==A by doing a matrix-vector product */
void check(const spasm_triplet *T, const spasm *A)
{
	int n = T->n;
	int m = T->m; 
	int nnz = T->nz;
	
	// safety checks
	assert(A->n == n);
	assert(A->m == m);
	const int *Ap = A->p;
	const int *Aj = A->j;
	assert(Ap[0] == 0);
	for (int i = 0; i < n; i++)
		assert(Ap[i] <= Ap[i + 1]);
	assert(Ap[n] == nnz);
	for (int k = 0; k < nnz; k++) {
		assert(0 <= Aj[k]);
		assert(Aj[k] < m);
	}

	// matrix-vector product
	double * X = malloc(n * sizeof(double));
	double * Ya = malloc(m * sizeof(double));
	double * Yb = malloc(m * sizeof(double));
	for (int i = 0; i < n; i++)
		X[i] = drand48();
	spasm_triplet_gemv(T, X, Ya);
	spasm_csr_gemv(A, X, Yb);

	double error = 0;
	for (int j = 0; j < m; j++) {
		double x = Ya[j] - Yb[j];
		error += x * x;
	}
	if (error > 1e-3)
		fprintf(stderr, "error = %f\n", error);
	assert(error < 1e-3);
	free(X);
	free(Ya);
	free(Yb);
}

/* benchmark the "classical" algorithm */
void run_test_classical(const spasm_triplet *T, const spasm_triplet *R)
{
	double start, stop;
	int n = T->n;
	int m = T->m;
	int nnz = T->nz;

	spasm *A = spasm_csr_alloc(n, m, nnz);
	spasm *B = spasm_csr_alloc(m, n, nnz);
	spasm *C = spasm_csr_alloc(m, n, nnz);
	spasm *D = spasm_csr_alloc(n, m, nnz);
	int *W = spasm_malloc((spasm_max(n, m) + 1) * sizeof(*W));

	start = spasm_wtime();
	classical_compress(T, A, W);
	stop = spasm_wtime();
	check(T, A);
	total[GUSTAVSON].compress += stop - start;
	fprintf(stderr, "-- Gustavson compress [COO->CSR]: %.3fs\n", stop - start);

	start = spasm_wtime();
	classical_transpose(A, B, W);
	stop = spasm_wtime();
	check(R, B);
	total[GUSTAVSON].transpose += stop - start;
	fprintf(stderr, "-- Gustavson transpose [CSR->CSR']: %.3fs\n", stop - start);

	start = spasm_wtime();
	classical_compress(R, C, W);
	stop = spasm_wtime();
		
	check(R, C);
	total[GUSTAVSON].compress_tr += stop - start;
	fprintf(stderr, "-- Gustavson compress [COO'->CSR']: %.3fs\n", stop - start);

	start = spasm_wtime();
	classical_transpose(C, D, W);
	stop = spasm_wtime();
		
	check(T, D);
	total[GUSTAVSON].transpose_tr += stop - start;
	fprintf(stderr, "-- Gustavson transpose [CSR'->CSR]: %.3fs\n", stop - start);

	spasm_csr_free(A);
	spasm_csr_free(B);
	spasm_csr_free(C);
	spasm_csr_free(D);
	free(W);
}

void run_test_stdsort(const spasm_triplet *T, const spasm_triplet *R)
{
	double start, stop;
	int n = T->n;
	int m = T->m;
	int nnz = T->nz;

	spasm *A = spasm_csr_alloc(n, m, nnz);
	spasm *B = spasm_csr_alloc(m, n, nnz);
	spasm *C = spasm_csr_alloc(m, n, nnz);
	spasm *D = spasm_csr_alloc(n, m, nnz);
	struct matrix_entry_t * Te = (struct matrix_entry_t *) spasm_malloc(nnz * sizeof(*Te));

	start = spasm_wtime();
	stdsort_compress(T, A, Te);
	stop = spasm_wtime();
	check(T, A);
	total[STDSORT].compress += stop - start;
	fprintf(stderr, "-- std::sort compress [COO->CSR]: %.3fs\n", stop - start);

	start = spasm_wtime();
	stdsort_transpose(A, B, Te);
	stop = spasm_wtime();
	check(R, B);
	total[STDSORT].transpose += stop - start;
	fprintf(stderr, "-- std::sort transpose [CSR->CSR']: %.3fs\n", stop - start);

	start = spasm_wtime();
	stdsort_compress(R, C, Te);
	stop = spasm_wtime();
		
	check(R, C);
	total[STDSORT].compress_tr += stop - start;
	fprintf(stderr, "-- std::sort compress [COO'->CSR']: %.3fs\n", stop - start);

	start = spasm_wtime();
	stdsort_transpose(C, D, Te);
	stop = spasm_wtime();	
	check(T, D);
	total[STDSORT].transpose_tr += stop - start;
	fprintf(stderr, "-- std::sort transpose [CSR'->CSR]: %.3fs\n", stop - start);

	spasm_csr_free(A);
	spasm_csr_free(B);
	spasm_csr_free(C);
	spasm_csr_free(D);
}

#ifdef HAVE_TBB
void run_test_tbbsort(const spasm_triplet *T, const spasm_triplet *R)
{
	double start, stop;
	int n = T->n;
	int m = T->m;
	int nnz = T->nz;

	spasm *A = spasm_csr_alloc(n, m, nnz);
	spasm *B = spasm_csr_alloc(m, n, nnz);
	spasm *C = spasm_csr_alloc(m, n, nnz);
	spasm *D = spasm_csr_alloc(n, m, nnz);
	struct matrix_entry_t * Te = (struct matrix_entry_t *) spasm_malloc(nnz * sizeof(*Te));

	start = spasm_wtime();
	tbbsort_compress(T, A, Te, -1);
	stop = spasm_wtime();
	check(T, A);
	total[TBBSORT].compress += stop - start;
	fprintf(stderr, "-- tbb::parallel_sort compress [COO->CSR]: %.3fs\n", stop - start);

	start = spasm_wtime();
	tbbsort_transpose(A, B, Te, -1);
	stop = spasm_wtime();
	check(R, B);
	total[TBBSORT].transpose += stop - start;
	fprintf(stderr, "-- tbb::parallel_sort transpose [CSR->CSR']: %.3fs\n", stop - start);

	start = spasm_wtime();
	tbbsort_compress(R, C, Te, -1);
	stop = spasm_wtime();
		
	check(R, C);
	total[TBBSORT].compress_tr += stop - start;
	fprintf(stderr, "-- tbb::parallel_sort compress [COO'->CSR']: %.3fs\n", stop - start);

	start = spasm_wtime();
	tbbsort_transpose(C, D, Te, -1);
	stop = spasm_wtime();	
	check(T, D);
	total[TBBSORT].transpose_tr += stop - start;
	fprintf(stderr, "-- tbb::parallel_sort transpose [CSR'->CSR]: %.3fs\n", stop - start);

	spasm_csr_free(A);
	spasm_csr_free(B);
	spasm_csr_free(C);
	spasm_csr_free(D);
	free(Te);
}
#endif

#ifdef HAVE_MKL
void run_test_MKL(const spasm_triplet *T)
{
	double start, stop;
	sparse_matrix_t mkl_T, mkl_A, mkl_B, mkl_C, mkl_D;
	sparse_status_t status;

	assert(sizeof(MKL_INT) == sizeof(int));

	status = mkl_sparse_d_create_coo(&mkl_T, SPARSE_INDEX_BASE_ZERO, T->n, T->m, T->nz, T->i, T->j, T->x);
	if (status != SPARSE_STATUS_SUCCESS)
		errx(1, "MKL create_coo (T) failed");

	// COO --> CSR
	start = spasm_wtime();
	status = mkl_sparse_convert_csr(mkl_T, SPARSE_OPERATION_NON_TRANSPOSE, &mkl_A);
	stop = spasm_wtime();
	if (status != SPARSE_STATUS_SUCCESS)
		errx(1, "MKL convert_csr (coo->csr) failed");
	total[MKL].compress += stop - start;
	fprintf(stderr, "-- MKL compress [COO->CSR]: %.3fs\n", stop - start);

	start = spasm_wtime();
	status = mkl_sparse_convert_csr(mkl_A, SPARSE_OPERATION_TRANSPOSE, &mkl_B);
	stop = spasm_wtime();
	if (status != SPARSE_STATUS_SUCCESS)
		errx(1, "MKL convert_csr (csr->csr) failed");
	total[MKL].transpose += stop - start;
	fprintf(stderr, "-- MKL transpose [CSR->CSR']: %.3fs\n", stop - start);

	start = spasm_wtime();
	status = mkl_sparse_convert_csr(mkl_T, SPARSE_OPERATION_TRANSPOSE, &mkl_C);
	stop = spasm_wtime();
	if (status != SPARSE_STATUS_SUCCESS)
		errx(1, "MKL convert_csr (coo->csr) failed");
	total[MKL].compress_tr += stop - start;
	fprintf(stderr, "-- MKL compress [COO'->CSR']: %.3fs\n", stop - start);


	start = spasm_wtime();
	status = mkl_sparse_convert_csr(mkl_C, SPARSE_OPERATION_TRANSPOSE, &mkl_D);
	stop = spasm_wtime();
	if (status != SPARSE_STATUS_SUCCESS)
		errx(1, "MKL convert_csr (csr->csr) failed");
	total[MKL].transpose_tr += stop - start;
	fprintf(stderr, "-- MKL transpose [CSR'->CSR]: %.3fs\n", stop - start);

	status = mkl_sparse_destroy(mkl_T);
	if (status != SPARSE_STATUS_SUCCESS)
		errx(1, "MKL destroy (T) failed");
	status = mkl_sparse_destroy(mkl_A);
	if (status != SPARSE_STATUS_SUCCESS)
		errx(1, "MKL destroy (A) failed");
	status = mkl_sparse_destroy(mkl_B);
	if (status != SPARSE_STATUS_SUCCESS)
		errx(1, "MKL destroy (B) failed");
	status = mkl_sparse_destroy(mkl_C);
	if (status != SPARSE_STATUS_SUCCESS)
		errx(1, "MKL destroy (C) failed");
	status = mkl_sparse_destroy(mkl_D);
	if (status != SPARSE_STATUS_SUCCESS)
		errx(1, "MKL destroy (D) failed");
}
#endif

void run_test(char *filename)
{
	FILE *f = fopen(filename, "r");
	if (f == NULL)
		err(1, "impossible to open %s", filename);

	spasm_triplet *T = spasm_load_mm(f);
	fclose(f);
		
	// transpose of T.
	spasm_triplet R;
	R.i = T->j;
	R.j = T->i;
	R.x = T->x;
	R.x = T->x;
	R.n = T->m;
	R.m = T->n;
	R.nz = T->nz;
	R.nzmax = T->nzmax;

	run_test_classical(T, &R);
	run_test_stdsort(T, &R);

#ifdef HAVE_TBB
	run_test_tbbsort(T, &R);
#endif

#ifdef HAVE_MKL
	run_test_MKL(T);
#endif

	fflush(stdout);
	spasm_triplet_free(T);
}

void show_grand_totals()
{
	fprintf(stderr, "\n\nGRAND TOTALS:\n");
	fprintf(stderr, "  Gustavson:\n");
	fprintf(stderr, "    compress:      %.3fs\n", total[GUSTAVSON].compress);
	fprintf(stderr, "    compress_tr:   %.3fs\n", total[GUSTAVSON].compress_tr);
	fprintf(stderr, "    transpsose:    %.3fs\n", total[GUSTAVSON].transpose);
	fprintf(stderr, "    transpsose_tr: %.3fs\n", total[GUSTAVSON].transpose_tr);
	
	fprintf(stderr, "  tbb::parallel_sort:\n");
	fprintf(stderr, "    compress:      %.3fs\n", total[TBBSORT].compress);
	fprintf(stderr, "    compress_tr:   %.3fs\n", total[TBBSORT].compress_tr);
	fprintf(stderr, "    transpsose:    %.3fs\n", total[TBBSORT].transpose);
	fprintf(stderr, "    transpsose_tr: %.3fs\n", total[TBBSORT].transpose_tr);
	
	fprintf(stderr, "  std::sort:\n");
	fprintf(stderr, "    compress:      %.3fs\n", total[STDSORT].compress);
	fprintf(stderr, "    compress_tr:   %.3fs\n", total[STDSORT].compress_tr);
	fprintf(stderr, "    transpsose:    %.3fs\n", total[STDSORT].transpose);
	fprintf(stderr, "    transpsose_tr: %.3fs\n", total[STDSORT].transpose_tr);

	fprintf(stderr, "  MKL:\n");
	fprintf(stderr, "    compress:      %.3fs\n", total[MKL].compress);
	fprintf(stderr, "    compress_tr:   %.3fs\n", total[MKL].compress_tr);
	fprintf(stderr, "    transpsose:    %.3fs\n", total[MKL].transpose);
	fprintf(stderr, "    transpsose_tr: %.3fs\n", total[MKL].transpose_tr);
}

int main()
{
#ifdef HAVE_TBB
	tbb_version();
#endif
#ifdef HAVE_MKL
	mkl_version();
#endif
#ifdef BENCHMARK_SMALL_MATRICES
	for (int i = 0; i < N; i++) {
	 	char filename[128];
	 	sprintf(filename, "%s/%s.mtx", MATRIX_PATH, matrices[i]);
		
		printf("#---------------------------------------- %s\n", filename);
		run_test(filename);

		fprintf(stderr, "\n\n");
	}
	
	show_grand_totals();
#endif


#ifdef BENCHMARK_LARGE_MATRICES
	for (int i = 1; i <= 58; i++) {  // just this one is enough to exhibit the crash
	 	char filename[128];
	 	sprintf(filename, "%s/RSA/pre_transpose%d.mtx", MATRIX_PATH, i);
		
		printf("#---------------------------------------- %s\n", filename);
		run_test(filename);

		fprintf(stderr, "\n\n");
		
	}
	show_grand_totals();
#endif

	return EXIT_SUCCESS;
}
