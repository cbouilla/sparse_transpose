///
/// \file simple_sort.h
/// \author Charles Bouillaguet and Jérôme Bonacchi
/// \brief This files uses std::sort and tbb::parallel_sort algorithms to
/// convert and to tranpose matrices.
/// \date 2020-07-09
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_DRIVER_SIMPLE_SORT_H
#define INCLUDE_DRIVER_SIMPLE_SORT_H

#include "mini_spasm.h"
#include "tools.h"

///
/// \brief Converts a sparse matrix in triplet format into a matrix in CSR
/// format.
///
/// \param[in] n the number of rows
/// \param[in] nnz the number of nonzero entries
/// \param[in] Te the input matrix in triplet format to convert
/// \param[out] A the output matrix in CSR format
///
void finalize(const u32 n, const u32 nnz, const struct matrix_entry_t *Te, spasm *A);

///
/// \brief Converts a matrix in triplet format into a matrix in CSR
/// format by using std::sort.
///
/// \param[in] T the input matrix in triplet format to convert
/// \param[out] A the output matrix in CSR format
/// \param[out] Te the output matrix in another triplet format
///
void stdsort_compress(const spasm_triplet *T, spasm *A,
                      struct matrix_entry_t *Te);

///
/// \brief Transposes a mtrix in CSR format by using std::sort.
/// Fist, converts a matrix in CSR format into a matrix in triplet
/// format. Then, uses std::sort to transpose it. Finally, converts it into a
/// matrix in CSR format.
///
/// \param[in] A the input matrix in CSR format to transpose
/// \param[out] R the output matrix in CSR format
/// \param[out] Te the output matrix in another triplet format
///
void stdsort_transpose(const spasm *A, spasm *R, struct matrix_entry_t *Te);

///
/// \brief Converts a matrix in triplet format into a matrix in CSR
/// format by using tbb::parallel_sort.
///
/// \param[in] T the input matrix in triplet format to convert
/// \param[out] A the output matrix in CSR format
/// \param[out] Te the output matrix in another triplet format
/// \param[in] num_threads the number of threads used in TBB's parallel sort
///
void tbbsort_compress(const spasm_triplet *T, spasm *A,
                      struct matrix_entry_t *Te, const u32 num_threads);

///
/// \brief Transposes a mtrix in CSR format by using tbb::parallel_sort.
/// Fist, converts a matrix in CSR format into a matrix in triplet
/// format. Then, uses tbb::parallel_sort to transpose it. Finally, converts it
/// into a matrix in CSR format.
///
/// \param[in] A the input matrix in CSR format to transpose
/// \param[out] R the output matrix in CSR format
/// \param[out] Te the output matrix in another triplet format
/// \param[in] num_threads the number of threads used in TBB's parallel sort
///
void tbbsort_transpose(const spasm *A, spasm *R, struct matrix_entry_t *Te,
                       const u32 num_threads);

///
/// \brief Prints the version of TBB used for compilation and used at runtime.
///
void tbb_version();

#endif /* INCLUDE_DRIVER_SIMPLE_SORT_H */
