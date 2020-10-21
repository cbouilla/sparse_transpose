///
/// \file classical_sort.h
/// \author Charles Bouillaguet (Github: cbouilla) and Jérôme Bonacchi (Github:
/// MarsParallax)
/// \brief This file contains the classical (Gustavson's) algorithm for
/// matrix conversion from COO format to CSR format and matrix transposition
/// in CSR format.
/// \date 2020
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_DRIVER_CLASSICAL_SORT_H
#define INCLUDE_DRIVER_CLASSICAL_SORT_H

#include "sparse.h"

///
/// \brief Converts a matrix in COO format into a matrix in CSR format by
/// using the classical algorithm.
///
/// \param[in] T the input matrix in COO format
/// \param[out] A the output matrix in CSR format
/// \param[out] W a scratch space (size == // TODO #columns + 1)
///
void classical_compress(const mtx_COO *T, mtx_CSR *A, u32 *W);

///
/// \brief Transposes a matrix in CSR format by using the classical
/// (Gustavson's) algorithm.
///
/// \param[in] A the input matrix in CSR format
/// \param[out] R the output matrix in CSR format
/// \param[out] W a scratch space (size == // TODO #columns + 1)
///
void classical_transpose(const mtx_CSR *A, mtx_CSR *R, u32 *W);

///
/// \brief Transposes a matrix in CSR format by using the altered version of the
/// classical algorithm: Gustavon's algorithm variant using an extra array.
/// Faster than the classical algorithm.
///
/// \param[in] A the input matrix in CSR format
/// \param[out] R the output matrix in CSR format
/// \param[out] Z a scratch space (size == nnz)
///
void wang_transpose(const mtx_COO *T, mtx_CSR *A, u32 *Z);

#endif /* INCLUDE_DRIVER_CLASSICAL_SORT_H */
