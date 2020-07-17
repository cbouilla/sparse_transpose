///
/// \file classical_sort.h
/// \author Charles Bouillaguet & Jérôme Bonacchi
/// \date 2020-07-09
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_DRIVER_CLASSICAL_SORT_H
#define INCLUDE_DRIVER_CLASSICAL_SORT_H

#include "mini_spasm.h"

void classical_compress(const spasm_triplet *T, spasm *A, int *W);
void classical_transpose(const spasm *A, spasm *R, int *W);

#endif /* INCLUDE_DRIVER_CLASSICAL_SORT_H */
