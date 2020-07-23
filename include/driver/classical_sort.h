///
/// \file classical_sort.h
/// \author Charles Bouillaguet & Jérôme Bonacchi
/// \brief Classical algorithms to convert et tranpose matrices
/// \date 2020-07-09
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_DRIVER_CLASSICAL_SORT_H
#define INCLUDE_DRIVER_CLASSICAL_SORT_H

#include "mini_spasm.h"

///
/// \brief Converts COO format matrix into CSR format matrix using Gustavon's
/// algorithm.
///
/// \param T Input matrix in COO format
/// \param A Output matrix in CSR format
/// \param W Scratch space, size == #columns + 1
///
void classical_compress(const spasm_triplet *T, spasm *A, int *W);

///
/// \brief Transposes CSR format matrix using Gustavon's algorithm.
///
/// \param A Input matrix in CSR format
/// \param R Output matrix in CSR format
/// \param W Working vector
///
void classical_transpose(const spasm *A, spasm *R, int *W);

#endif /* INCLUDE_DRIVER_CLASSICAL_SORT_H */
