///
/// \file classical_sort.h
/// \author Charles Bouillaguet and Jérôme Bonacchi
/// \brief This file implements the classical algorithms to convert and to
/// tranpose matrices.
/// \date 2020-07-09
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_DRIVER_CLASSICAL_SORT_H
#define INCLUDE_DRIVER_CLASSICAL_SORT_H

#include <stdint.h>

///
/// \brief Type used for dimensions and entries.
///
typedef uint32_t u32;

#include "mini_spasm.h"

///
/// \brief Converts a matrix in triplet format into a matrix in CSR format by
/// using the classical algorithm.
///
/// \param[in] T the input matrix in triplet format
/// \param[out] A the output matrix in CSR format
/// \param[out] W a scratch space (size == // TODO #columns + 1)
///
void classical_compress(const spasm_triplet *T, spasm *A, u32 *W);

///
/// \brief Transposes a matrix in CSR format by using the classical algorithm.
///
/// \param[in] A the input matrix in CSR format
/// \param[out] R the output matrix in CSR format
/// \param[out] W a scratch space (size == // TODO #columns + 1)
///
void classical_transpose(const spasm *A, spasm *R, u32 *W);

#endif /* INCLUDE_DRIVER_CLASSICAL_SORT_H */
