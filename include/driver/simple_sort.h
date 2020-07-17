///
/// \file simple_sort.h
/// \author Charles Bouillaguet & Jérôme Bonacchi
/// \date 2020-07-09
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_DRIVER_SIMPLE_SORT_H
#define INCLUDE_DRIVER_SIMPLE_SORT_H

#include "mini_spasm.h"

void stdsort_compress(const spasm_triplet *T, spasm *A,
                      struct matrix_entry_t *Te);
void stdsort_transpose(const spasm *A, spasm *R, struct matrix_entry_t *Te);

void tbbsort_compress(const spasm_triplet *T, spasm *A,
                      struct matrix_entry_t *Te, int num_threads);
void tbbsort_transpose(const spasm *A, spasm *R, struct matrix_entry_t *Te,
                       int num_threads);
void tbb_version();

#endif /* INCLUDE_DRIVER_SIMPLE_SORT_H */
