#include <algorithm>
#include <cstdint>

/* simple compress/transpose algorithm that directly rely on explicit sorting, 
   using off-the-shelf sort functions (std::sort from the STL, or Intel's parallel quicksort). 
*/

extern "C" {
	#include "mini_spasm.h"
}

extern "C" void stdsort_compress(const spasm_triplet * T, spasm * A, struct matrix_entry_t * Te);
extern "C" void stdsort_transpose(const spasm * A, spasm * R, struct matrix_entry_t * Te);


#ifdef HAVE_TBB
#include <tbb/parallel_sort.h>
#include <tbb/task_scheduler_init.h>
extern "C" void tbbsort_compress(const spasm_triplet * T, spasm * A, struct matrix_entry_t * Te, int num_threads);
extern "C" void tbbsort_transpose(const spasm * A, spasm * R, struct matrix_entry_t * Te, int num_threads);
#endif


inline bool operator < (const struct matrix_entry_t &a, const struct matrix_entry_t &b)
{
	return a.i < b.i;
}

static void finalize(int n, int nnz, const matrix_entry_t * Te, spasm * A)
{
	int *Ap = A->p;
	int *Aj = A->j;
	double *Ax = A->x;
	
	// scanning, copying and setting pointers
	// could be parallelized (bit this is not completely trivial)
	int i = -1;
	for (int k = 0; k < nnz; k++) {
		int next_i = Te[k].i;
		int next_j = Te[k].j;
		double next_x = Te[k].x;
		if (i < next_i) {
			i++;
			while (i < next_i) {
				Ap[i] = k;
				i++;
			}
			Ap[i] = k;
		}
		// append entry to current line
		Aj[k] = next_j;
		Ax[k] = next_x;
	}
	// finalization
	i += 1;
	while (i <= n) {
		Ap[i] = nnz;
		i += 1;
	}
}

void stdsort_compress(const spasm_triplet * T, spasm * A, struct matrix_entry_t * Te)
{
	int n = T->n;
	int nnz = T->nz;
	int *Ti = T->i;
	int *Tj = T->j;
	double *Tx = T->x;
	
	for (int k = 0; k < nnz; k++) {
		Te[k].i = Ti[k];
		Te[k].j = Tj[k];
		Te[k].x = Tx[k];
	}
	
	std::sort(Te, Te + nnz);
	finalize(n, nnz, Te, A);
}


void stdsort_transpose(const spasm * A, spasm * R, struct matrix_entry_t * Te)
{
	int n = A->n;
	int m = A->m;
	int nnz = spasm_nnz(A);
	const int *Ap = A->p;
	const int *Aj = A->j;
	const double *Ax = A->x;

	for (int i = 0; i < n; i++)
		for (int k = Ap[i]; k < Ap[i + 1]; k++) {
			int j = Aj[k];
			double x = Ax[k];
			Te[k].i = j;
			Te[k].j = i;
			Te[k].x = x;
		}
	
	std::sort(Te, Te + nnz);
	finalize(m, nnz, Te, R);
}

#ifdef HAVE_TBB
void tbbsort_compress(const spasm_triplet * T, spasm * A, struct matrix_entry_t * Te, int num_threads)
{
	tbb::task_scheduler_init tsi(num_threads);

	int n = T->n;
	int nnz = T->nz;
	int *Ti = T->i;
	int *Tj = T->j;
	double *Tx = T->x;

	double a = spasm_wtime();

	#pragma omp parallel for
	for (int k = 0; k < nnz; k++) {
		Te[k].i = Ti[k];
		Te[k].j = Tj[k];
		Te[k].x = Tx[k];
	}

	double b = spasm_wtime();
	
	tbb::parallel_sort(Te, Te + nnz);
	
	double c = spasm_wtime();

	finalize(n, nnz, Te, A);

	/* 
	// extra timing
	double d = spasm_wtime();
	printf("      Subtimes:\n");
	printf("        Fill_Te: %.3f\n", b-a);
	printf("        tbb::sort: %.3f\n", c-b);
	printf("        finalize: %.3f\n", d-c);
	*/
}

void tbbsort_transpose(const spasm * A, spasm * R, struct matrix_entry_t * Te, int num_threads)
{
	tbb::task_scheduler_init tsi(num_threads);

	int n = A->n;
	int m = A->m;
	int nnz = spasm_nnz(A);
	const int *Ap = A->p;
	const int *Aj = A->j;
	const double *Ax = A->x;

	double a = spasm_wtime();
	#pragma omp parallel for
	for (int i = 0; i < n; i++)
		for (int k = Ap[i]; k < Ap[i + 1]; k++) {
			int j = Aj[k];
			double x = Ax[k];
			Te[k].j = i;
			Te[k].i = j;
			Te[k].x = x;
		}
	
	double b = spasm_wtime();
	
	tbb::parallel_sort(Te, Te + nnz);
	
	double c = spasm_wtime();

	finalize(m, nnz, Te, R);
	/*
	// extra timing
	double d = spasm_wtime();
	printf("      Subtimes:\n");
	printf("        Fill_Te: %.3f\n", b-a);
	printf("        tbb::sort: %.3f\n", c-b);
	printf("        finalize: %.3f\n", d-c);
	*/
}
#endif