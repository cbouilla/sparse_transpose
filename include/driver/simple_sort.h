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

#include "sparse.h"

typedef uint32_t u32;

///
/// \brief Converts a sparse matrix in COO format into a matrix in CSR
/// format.
///
/// \param[in] n the number of rows
/// \param[in] nnz the number of nonzero entries
/// \param[in] Te the input matrix in COO format to convert
/// \param[out] A the output matrix in CSR format
///
void finalize(const u32 n, const u32 nnz, const mtx_entry *Te,
              mtx_CSR *A);

///
/// \brief Converts a matrix in COO format into a matrix in CSR
/// format by using std::sort.
///
/// \param[in] T the input matrix in COO format to convert
/// \param[out] A the output matrix in CSR format
/// \param[out] Te the output matrix in another COO format
///
void stdsort_compress(const mtx_COO *T, mtx_CSR *A, mtx_entry *Te);

///
/// \brief Transposes a mtrix in CSR format by using std::sort.
/// Fist, converts a matrix in CSR format into a matrix in COO
/// format. Then, uses std::sort to transpose it. Finally, converts it into a
/// matrix in CSR format.
///
/// \param[in] A the input matrix in CSR format to transpose
/// \param[out] R the output matrix in CSR format
/// \param[out] Te the output matrix in another COO format
///
void stdsort_transpose(const mtx_CSR *A, mtx_CSR *R, mtx_entry *Te);

///
/// \brief Converts a matrix in COO format into a matrix in CSR
/// format by using tbb::parallel_sort.
///
/// \param[in] T the input matrix in COO format to convert
/// \param[out] A the output matrix in CSR format
/// \param[out] Te the output matrix in another COO format
/// \param[in] num_threads the number of threads used in TBB's parallel sort
///
void tbbsort_compress(const mtx_COO *T, mtx_CSR *A, mtx_entry *Te,
                      const u32 num_threads);

///
/// \brief Transposes a mtrix in CSR format by using tbb::parallel_sort.
/// Fist, converts a matrix in CSR format into a matrix in COO
/// format. Then, uses tbb::parallel_sort to transpose it. Finally, converts it
/// into a matrix in CSR format.
///
/// \param[in] A the input matrix in CSR format to transpose
/// \param[out] R the output matrix in CSR format
/// \param[out] Te the output matrix in another COO format
/// \param[in] num_threads the number of threads used in TBB's parallel sort
///
void tbbsort_transpose(const mtx_CSR *A, mtx_CSR *R, mtx_entry *Te,
                       const u32 num_threads);

///
/// \brief Prints the version of TBB used for compilation and used at runtime.
///
void tbb_version();

#endif /* INCLUDE_DRIVER_SIMPLE_SORT_H */
