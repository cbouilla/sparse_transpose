///
/// \file classical_sort.c
/// \author Charles Bouillaguet & Jérôme Bonacchi
/// \brief Implementation of the "classical" sort algorithm from Gustavson.
/// \date 2020-07-09
/// 
/// @copyright Copyright (c) 2020
///

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <err.h>

#include "classical_sort.h"
#include "mini_spasm.h"

// assemble the CSR representation of T
// W = scratch space, size == #rows + 1
void classical_compress(const spasm_triplet * T, spasm * A, int *W)
{
	const int *Ti = T->i;
	const int *Tj = T->j;
	const double *Tx = T->x;
	int *Ap = A->p;
	int *Aj = A->j;
	double *Ax = A->x;
	
	int n = T->n;
	int nnz = T->nz;

	// setup
	for (int i = 0; i < n; i++)  /* gcc replaces this with AVX2-optimized memset */
		W[i] = 0;
	
	// counting entries on each row
	for (int k = 0; k < nnz; k++) {
		int i = Ti[k];
		int w = W[i];
		w += 1;
		W[i] = w;
	}

	// prefix-sum + setup Ap
	int s = 0;
	for (int i = 0; i < n; i++) {
		int w = W[i];
		Ap[i] = s;
		W[i] = s;
		s += w;
	}
	Ap[n] = s;

	// dispatch
	for (int k = 0; k < nnz; k++) {
		int i = Ti[k];
		int j = Tj[k];
		double x = Tx[k];
		int l = W[i];
		Aj[l] = j;
		Ax[l] = x;
		l += 1;
		W[i] = l;
	}
}

// W = scratch space, size == #columns + 1
void classical_transpose(const spasm * A, spasm * R, int *W)
{
	const int *Ap = A->p;
	const int *Aj = A->j;
	const double *Ax = A->x;
	int *Rp = R->p;
	int *Rj = R->j;
	double *Rx = R->x;
	
	int n = A->n;
	int m = A->m;

	// setup
	for (int i = 0; i < m; i++)  /* gcc replaces this with AVX2-optimized memset */
		W[i] = 0;
	
	// counting entries on each column
	for (int i = 0; i < n; i++)
		for (int k = Ap[i]; k < Ap[i + 1]; k++) {
			int j = Aj[k];
			int w = W[j];
			w += 1;
			W[j] = w;
	}

	// prefix-sum + setup Rp
	int s = 0;
	for (int j = 0; j < m; j++) {
		int w = W[j];
		Rp[j] = s;
		W[j] = s;
		s += w;
	}
	Rp[m] = s;

	// dispatch
	for (int i = 0; i < n; i++)
		for (int k = Ap[i]; k < Ap[i + 1]; k++) {
			int j = Aj[k];
			double x = Ax[k];
			int l = W[j];
			Rj[l] = i;
			Rx[l] = x;
			l += 1;
			W[j] = l;
		}
}

/* Wang et. al variant using the extra array. Faster than the classical variant...

void wang_transpose(const spasm_triplet * A, spasm * T, int *W, int *Z)
{
	(void) W;
	const int *Ai = A->i;
	const int *Aj = A->j;
	int *Tp = T->p;
	int *Tj = T->j;
	int m = A->m;
	int nnz = A->nz;

	for (int j = 0; j < m; j++)
		Tp[j] = 0;
	for (int i = 0; i < n; i++)
		for (int k = Ap[i]; k < Ap[i + 1]; k++) {
			int j = Aj[k];
			int w = Tp[j];
			Z[k] = w;
			w += 1;
			Tp[j] = w;
		}

	int s = 0;
	for (int j = 0; j < m; j++) {
		int w = Tp[j];
		Tp[j] = s;
		s += w;
	}
	Tp[m] = s;

	for (int k = 0; k < nnz; k++) {
		int j = Aj[k];
		int s = Z[k];
		int r = Tp[j];
		int i = Ai[k];
		int l = r + s;
		Tj[l] = i;
	}
}
*/
