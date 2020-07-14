///
/// \file simple_sort.cpp
/// \author Charles Bouillaguet & Jérôme Bonacchi
/// \brief simple compress/transpose algorithms that directly rely on explicit sorting,
///  using off-the-shelf sort functions (std::sort from the STL or Intel's parallel sort).
/// \date 2020-07-09
///
/// @copyright Copyright (c) 2020
///

#include <algorithm>
#include <cstdint>
#include <omp.h>

extern "C" {
	#include "mini_spasm.h"
}

extern "C" void stdsort_compress(const spasm_triplet * T, spasm * A, struct matrix_entry_t * Te);
extern "C" void stdsort_transpose(const spasm * A, spasm * R, struct matrix_entry_t * Te);

#ifdef HAVE_TBB
#include <iostream>
#include <tbb/parallel_sort.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tbb_stddef.h>

extern "C" void tbbsort_compress(const spasm_triplet * T, spasm * A, struct matrix_entry_t * Te, int num_threads);
extern "C" void tbbsort_transpose(const spasm * A, spasm * R, struct matrix_entry_t * Te, int num_threads);
extern "C" void tbb_version();

#endif // HAVE_TBB

///
/// \brief Compare the row of the matrix_entry_t data type.
///
/// \param a Left-hand side object
/// \param b Right-hand side object
/// \return true If the row index of a is less than the row index of b
/// \return false If the row index of a is greater than, or equal to, the row index of b
///
inline bool operator < (const struct matrix_entry_t &a, const struct matrix_entry_t &b)
{
	return a.i < b.i;
}

///
/// \brief Convert a sparse matrix in COO format into CSR format.
///
/// \param n Number of rows
/// \param nnz Number of non-zero entries
/// \param Te Sparse matrix in COO format to convert
/// \param A Sparse matrix in CSR format
///
static void finalize(int n, int nnz, const matrix_entry_t * Te, spasm * A)
{
	int *Ap = A->p;
	int *Aj = A->j;
	double *Ax = A->x;
	
	// scanning, copying and setting pointers
	// could be parallelized (but this is not completely trivial)
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

///
/// \brief Convert the sparse matrix T in COO format into A in CSR format by
/// using std::sort.
/// 
/// \param T Sparse matrix in COO format to convert
/// \param A T converted into CSR format
/// \param Te T converted into another COO format
///
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

///
/// \brief Convert the sparse matrix A in CSR format into Te in COO format.
/// Then, use std::sort to transpose Te. Finally, convert Te into R (CSR 
/// format)
/// 
/// \param A Sparse matrix in CSR format to transpose
/// \param R Sparse matrix in CSR format that is the transpose of A
/// \param Te Sparse matrix in COO format that is the transpose of A
/// \param num_threads Number of threads used in TBB's parallel sort
///
///
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

///
/// \brief Convert the sparse matrix T in COO format into A in CSR format by
/// using TBB's parallel_sort.
/// 
/// \param T Sparse matrix in COO format to convert
/// \param A T converted into CSR format
/// \param Te T converted into another COO format
/// \param num_threads Number of threads used in TBB's parallel sort
///
void tbbsort_compress(const spasm_triplet * T, spasm * A, struct matrix_entry_t * Te, int num_threads)
{
	tbb::task_scheduler_init tsi(num_threads);

	int n = T->n;
	int nnz = T->nz;
	int *Ti = T->i;
	int *Tj = T->j;
	double *Tx = T->x;

	// double a = spasm_wtime();
	omp_set_num_threads(num_threads);
	#pragma omp parallel for
	for (int k = 0; k < nnz; k++) {
		Te[k].i = Ti[k];
		Te[k].j = Tj[k];
		Te[k].x = Tx[k];
	}

	// double b = spasm_wtime();
	
	tbb::parallel_sort(Te, Te + nnz);
	
	// double c = spasm_wtime();

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

///
/// \brief Convert the sparse matrix A in CSR format into Te in COO format.
/// Then, use TBB's parallel_sort to transpose Te. Finally, convert Te into R
/// (CSR format)
/// 
/// \param A Sparse matrix in CSR format to transpose
/// \param R Sparse matrix in CSR format that is the transpose of A
/// \param Te Sparse matrix in COO format that is the transpose of A
/// \param num_threads Number of threads used in TBB's parallel sort
///
void tbbsort_transpose(const spasm * A, spasm * R, struct matrix_entry_t * Te, int num_threads)
{
	tbb::task_scheduler_init tsi(num_threads);

	int n = A->n;
	int m = A->m;
	int nnz = spasm_nnz(A);
	const int *Ap = A->p;
	const int *Aj = A->j;
	const double *Ax = A->x;

	// double a = spasm_wtime();
	omp_set_num_threads(num_threads);
	#pragma omp parallel for
	for (int i = 0; i < n; i++)
		for (int k = Ap[i]; k < Ap[i + 1]; k++) {
			int j = Aj[k];
			double x = Ax[k];
			Te[k].j = i;
			Te[k].i = j;
			Te[k].x = x;
		}
	
	// double b = spasm_wtime();
	
	tbb::parallel_sort(Te, Te + nnz);
	
	// double c = spasm_wtime();

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

///
/// \brief Print the version of TBB use for compilation and at runtime.
///
void tbb_version()
{
	std::cout << "TBB version: compiled=" << TBB_INTERFACE_VERSION
	<< ", runtime=" << tbb::TBB_runtime_interface_version() << '\n';
}

#endif // HAVE_TBB
