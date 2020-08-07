///
/// \file tools.h
/// \author Charles Bouillaguet and Jérôme Bonacchi
/// \brief This files containes `#define` and variables useful for many other
/// files.
/// \date 2020-07-21
///
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_TOOLS_H
#define INCLUDE_TOOLS_H

#include <stdint.h>

#define N_REPEAT 7
#define BENCHMARK_SMALL_MATRICES
#define BENCHMARK_LARGE_MATRICES

/// The offsets in the "total" array below
#define GUSTAVSON 0
#define MKL 1
#define STDSORT 2
#define TBBSORT 3
#define SCANTRANS 4
#define MERGETRANS 5
#define BB 6
#define N_METHOD 7

#define OUTPUT_FILENAME "benchmarks.csv"

/// name, cflags, cxxflags, threads, matrix, compress, transpose, compress_tr,
/// transpose_tr, step, total_step
#define OUTPUT_FORMAT "%s,%s,%s,%d, %s,%.15f,%.15f,%.15f,%.15f,%d,%d\n"

#define MATRIX_PATH "/Infos/lmd/2019/master/ue/MU4IN903-2020fev"

#define N_SMALL_MATRICES 21

/// The matrices in "Parallel Transposition of Sparse Data Structures" by Wang,
/// Liu, Hou and Feng
static const char *matrices[N_SMALL_MATRICES] = {
	"language",
	"ASIC_680k",
	"circuit5M",
	"flickr",
	"memchip",
	"rajat21",
	"sme3Dc",
	"stomach",
	"transient",
	"webbase-1M",
	"wiki-Talk",
	"cage14",
	"eu-2005",
	"FullChip",
	"mac_econ_fwd500",
	"para-4",
	"rajat29",
	"Stanford_Berkeley",
	"torso1",
	"venkat01",
	"web-Google"
};

#define N_LARGE_MATRICES 58

typedef uint8_t u8;
typedef uint32_t u32;

///
/// \brief Type used to time algorithms.
///
typedef struct {
	double compress;     ///< duration to "compress" the matrix    (convert   COO  -> CSR) 
	double compress_tr;  ///< duration to "compress" the transpose (convert   COO' -> CSR') 
	double transpose;    ///< duration to transpose the matrix     (transpose CSR  -> CSR')
	double transpose_tr; ///< duration to transpose the transpose  (transpose CSR' -> CSR)
} algorithm_times;

static algorithm_times total[N_METHOD];

///
/// \brief Reset to 0 the values of `duration`.
/// 
/// \param[in, out] duration the data to reset
///
inline void clear_times(algorithm_times *duration)
{
  duration->compress = 0.0;
  duration->compress_tr = 0.0;
  duration->transpose = 0.0;
  duration->transpose_tr = 0.0;
}

#endif /* INCLUDE_TOOLS_H */
