# Sparse Matrix Transposition

- [Sparse Matrix Transposition](#sparse-matrix-transposition)
  - [Installing](#installing)
  - [Matrices](#matrices)
  - [Use](#use)
    - [Compiler](#compiler)
    - [Intel Math Kernel Library (MKL)](#intel-math-kernel-library-mkl)
    - [Intel Threading Building Block (TBB)](#intel-threading-building-block-tbb)
  - [Version](#version)
  - [Documentation](#documentation)
  - [TODO](#todo)

## Installing

Tested with:

- Python 3.8.3: [PSF licence agreement](https://docs.python.org/3/license.html)
- Matplotlib 3.2.1: [The Matplotlib Development Team](https://matplotlib.org/3.2.1/users/license.html)
- NumPy 1.18.5: [Copyright (c) 2005, NumPy Developers](https://numpy.org/doc/stable/license.html)
- pandas 1.0.3: [BSD 3-Clause License](https://pandas.pydata.org/pandas-docs/stable/getting_started/overview.html)
- spTrans: [Wang, Hao and Liu, Weifeng and Hou, Kaixi and Feng, Wu-chun, LGPL-2.1](https://github.com/vtsynergy/sptrans)
- SpaSM: [The SpaSM group, GPL-3.0](https://github.com/cbouilla/spasm/)
- Intel Math Kernel Library (MKL) 2020.0.166: [Intel Simplified Software License](https://software.intel.com/content/www/us/en/develop/articles/end-user-license-agreement.html#inpage-nav-3)
- Intel Threading Building Block (TBB) 2020.0.166: [Intel Simplified Software License](https://software.intel.com/content/www/us/en/develop/articles/end-user-license-agreement.html#inpage-nav-3)
- ICC
- C 11
- C 11 standard librairies
- C++ 11
- C++ 11 STL
- GCC 6.3.0 20170516 (Debian 6.3.0-18+deb9u1) and 10.1.0

## Matrices

Information can be found [here](matrices.md).

## Use

One can run `init.sh` to load the MKL, TBB, ICC/ICPC:

```sh
source init.sh
```

Paths must be changed in order to your own configuration.

### Compiler

When running `make` one can use :

```sh
make [CC=[gcc | icc]] [CXX=[g++ | icpc]]
```

### Intel Math Kernel Library (MKL)

In order to use the MKL, one should run:

```sh
source /path/to/mklvars.sh intel64 [ilp64]
```

The last option activates 64-bits integers (not tested, may not work).

Then, when running `make` one can use :

```sh
make [USE_MKL=[tbb | gomp | iomp | serial]]
```

For more information about how to use the MKL, one can refer to [Intel® Math Kernel Library Link Line Advisor](https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html).

### Intel Threading Building Block (TBB)

In order to use the TBB library, one should run:

```sh
source /path/to/tbbvars.sh intel64
```

Then, when running `make` one can use :

```sh
make [USE_TBB=?]
```

## Version

1. count + prefix-sum dans transpose()
2. finalize dans transpose()
3. lo_bucket -> hi_bucket dans transpose_bucket()
4. WCB et OUT en struct entry_t *
5. 2 lignes de cache + OMP sur Rp avec n = 1
6. 3. et OMP sur Rp avec n = 1
7. 6. et maximisation du radix
8. 4. et écriture de Rj et Rx dans la boucle de Rp
9. 8. et test/print
10. 3. avec variation du radix selon nnz

## Documentation

Documentation can be generated thanks to Doxygen by running:

```sh
make docs
```

## TODO

- [ ] Remove warnings
- [ ] cite/licence SpaSM/CADO-NFS ?
- [ ] rename versions, keep 3 and 10
- [ ] comment code and make the doc (comments may be incorrect)
- [ ] enhance prints
- [ ] check todos (classical_sort.h, sparse.h, sparse.c, transpose*.h, transpose.c)
- [ ] parallelise and vectorise with AVX128 AVX256 AVX512
- [ ] rename files with "::" in their name
- [ ] refactor code
- [ ] compute memory use of Wang, scanTrans with 1 thread and Gustavson
- [ ] parallel LSD radix sort
- [ ] parallel MSD radix sort
- [ ] radix sort:  use another sort when their only a few elements in buckets, use half flush/purge, reduce conditionnal branch
