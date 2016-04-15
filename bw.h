#include "m4ri.h"
#include "flint.h"
#include "nmod_sparse_mat.h"

void
nmod_sparse_mat_mul_m4ri_mat(mzd_t* C, const nmod_sparse_mat_t A, const mzd_t *B);

void
_bw(mzd_t *K, const nmod_sparse_mat_t M, const int skip, const int epsilon, const int m_shift, const int n_shift);

void
bw(mzd_t *K, const nmod_sparse_mat_t M);
