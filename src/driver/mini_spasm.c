#include <assert.h>
#include <err.h>
#include <math.h>
#include <sys/time.h>

#include "mini_spasm.h"
#include "mmio.h"

double spasm_wtime()
{
  struct timeval ts;
  gettimeofday(&ts, NULL);
  return (double)ts.tv_sec + ts.tv_usec / 1E6;
}

int spasm_nnz(const spasm *A) { return A->p[A->n]; }

void *spasm_malloc(size_t size)
{
  void *x = aligned_alloc(64, size);
  if (x == NULL)
    err(1, "malloc failed");
  return x;
}

/* return a string representing n in 4 bytes */
static void spasm_human_format(int64_t n, char *target)
{
  if (n < 1000)
  {
    sprintf(target, "%d", (int)n);
    return;
  }
  if (n < 1000000)
  {
    sprintf(target, "%.1fk", n / 1e3);
    return;
  }
  if (n < 1000000000)
  {
    sprintf(target, "%.1fM", n / 1e6);
    return;
  }
  if (n < 1000000000000ll)
  {
    sprintf(target, "%.1fG", n / 1e9);
    return;
  }
  if (n < 1000000000000000ll)
  {
    sprintf(target, "%.1fT", n / 1e12);
    return;
  }
}

/* free a sparse matrix */
void spasm_csr_free(spasm *A)
{
  if (A == NULL)
    return;
  free(A->p);
  free(A->j);
  free(A->x); /* trick : free does nothing on NULL pointer */
  free(A);
}

void spasm_triplet_free(spasm_triplet *T)
{
  free(T->i);
  free(T->j);
  free(T->x); /* trick : free does nothing on NULL pointer */
  free(T);
}

/* allocate a sparse matrix (triplet form) */
spasm_triplet *spasm_triplet_alloc(int nzmax)
{
  spasm_triplet *T;

  T = spasm_malloc(sizeof(spasm_triplet));
  T->m = 0;
  T->n = 0;
  T->nzmax = nzmax;
  T->nz = 0;
  T->i = spasm_malloc(nzmax * sizeof(int));
  T->j = spasm_malloc(nzmax * sizeof(int));
  T->x = spasm_malloc(nzmax * sizeof(double));
  return T;
}

/* allocate a sparse matrix (compressed-row form) */
spasm *spasm_csr_alloc(int n, int m, int nzmax)
{
  spasm *A;

  A = spasm_malloc(sizeof(spasm)); /* allocate the cs struct */
  A->m = m;                        /* define dimensions and nzmax */
  A->n = n;
  A->nzmax = nzmax;
  A->p = spasm_malloc((n + 1) * sizeof(A->p));
  A->j = spasm_malloc(nzmax * sizeof(A->j));
  A->x = spasm_malloc(nzmax * sizeof(A->x));
  return A;
}

/* add an entry to a triplet matrix; enlarge it if necessary */
static inline void spasm_add_entry(spasm_triplet *T, int i, int j, double x)
{
  assert((i >= 0) && (i < T->n) && (j >= 0) && (j < T->m));
  int nz = T->nz;
  assert(nz < T->nzmax);

  T->i[nz] = i;
  T->j[nz] = j;
  T->x[nz] = x;

  T->nz = nz + 1;
  // T->n = spasm_max(T->n, i + 1);
  // T->m = spasm_max(T->m, j + 1);
}

/*
 * Load a matrix in MatrixMarket sparse format.
 * Heavily inspired by the example program:
 *     http://math.nist.gov/MatrixMarket/mmio/c/example_read.c
 */
spasm_triplet *spasm_load_mm(FILE *f)
{
  MM_typecode matcode;
  int n, m, nnz;

  double start = spasm_wtime();
  if (mm_read_banner(f, &matcode) != 0)
    errx(1, "Could not process Matrix Market banner.\n");
  char *typecode = mm_typecode_to_str(matcode);

  if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode))
    errx(1, "Matrix Market type: [%s] not supported", typecode);

  if (!mm_is_general(matcode))
    errx(1, "Matrix market type [%s] not supported", typecode);

  if (mm_is_symmetric(matcode) || mm_is_skew(matcode))
    errx(1, "Matrix market type [%s] not supported", typecode);

  int real = mm_is_real(matcode) || mm_is_integer(matcode);
  int pattern = mm_is_pattern(matcode);

  if (!real && !pattern)
    errx(1, "Matrix market type [%s] not supported", typecode);

  if (mm_read_mtx_crd_size(f, &n, &m, &nnz) != 0)
    errx(1, "Cannot read matrix size");

  char s_nnz[16];
  spasm_human_format(nnz, s_nnz);
  fprintf(stderr, "[IO] loading %d x %d MTX [%s] %s NNZ ...", n, m, typecode,
          s_nnz);
  fflush(stderr);
  free(typecode); // mm_typecode_to_str return a malloc'd char*

  spasm_triplet *T = spasm_triplet_alloc(nnz);
  T->n = n;
  T->m = m;
  for (int i = 0; i < nnz; i++)
  {
    int u, v;
    double x;
    if (real)
    {
      if (3 != fscanf(f, "%d %d %lg\n", &u, &v, &x))
        errx(1, "parse error entry %d\n", i);
      spasm_add_entry(T, u - 1, v - 1, x);
    }
    else
    {
      if (2 != fscanf(f, "%d %d\n", &u, &v))
        errx(1, "parse error entry %d\n", i);
      spasm_add_entry(T, u - 1, v - 1, 1.0);
    }
  }
  fprintf(stderr, " Actually %d x %d [%.1fs]\n", T->n, T->m,
          spasm_wtime() - start);
  return T;
}

void spasm_triplet_gemv(const spasm_triplet *T, const double *x, double *y)
{
  const int *Ti = T->i;
  const int *Tj = T->j;
  const double *Tx = T->x;
  int m = T->m;
  int nnz = T->nz;

  for (int j = 0; j < m; j++)
    y[j] = 0;
  for (int k = 0; k < nnz; k++)
  {
    int i = Ti[k];
    int j = Tj[k];
    y[j] += Tx[k] * x[i];
  }
}

void spasm_csr_gemv(const spasm *A, const double *x, double *y)
{
  const int *Ap = A->p;
  const int *Aj = A->j;
  const double *Ax = A->x;
  int n = A->n;
  int m = A->m;

  for (int j = 0; j < m; j++)
    y[j] = 0;
  for (int i = 0; i < n; i++)
    for (int k = Ap[i]; k < Ap[i + 1]; k++)
    {
      int j = Aj[k];
      y[j] += Ax[k] * x[i];
    }
}

void spasm_triplet_transpose(const spasm_triplet *T, spasm_triplet *R)
{
  R->i = T->j;
  R->j = T->i;
  R->x = T->x;
  R->n = T->m;
  R->m = T->n;
  R->nz = T->nz;
  R->nzmax = T->nzmax;
}
