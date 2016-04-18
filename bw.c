#include <stdio.h>
#include <gmp.h>
#include "m4ri.h"
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_sparse_mat.h"
#define MZD_MUL_CUTOFF 0

void
nmod_sparse_mat_mul_m4ri_mat_w1(mzd_t* C, const nmod_sparse_mat_t A, const mzd_t *B)
{
    int i, k;
    for (i = 0; i < A->r; i++)
    {
        for (k = 0; k < A->row_supports[i]; k++)
            C->rows[i][0] ^= B->rows[A->rows[i][k].pos][0];
    }
}
void
clever(word *restrict a, word *restrict b)
{
    a[0] ^= b[0];
    a[1] ^= b[1];
}

void
nmod_sparse_mat_mul_m4ri_mat_w2(mzd_t* C, const nmod_sparse_mat_t A, const mzd_t *B)
{
    int i, k;
    for (i = 0; i < A->r; i++)
    {
        for (k = 0; k < A->row_supports[i]; k++)
        {
            C->rows[i][0] ^= B->rows[A->rows[i][k].pos][0];
            C->rows[i][1] ^= B->rows[A->rows[i][k].pos][1];
        }
    }
}

void
nmod_sparse_mat_mul_m4ri_mat_w4(mzd_t* C, const nmod_sparse_mat_t A, const mzd_t *B)
{
    int i, k;
    for (i = 0; i < A->r; i++)
    {
        for (k = 0; k < A->row_supports[i]; k++)
        {
            C->rows[i][0] ^= B->rows[A->rows[i][k].pos][0];
            C->rows[i][1] ^= B->rows[A->rows[i][k].pos][1];
            C->rows[i][2] ^= B->rows[A->rows[i][k].pos][2];
            C->rows[i][3] ^= B->rows[A->rows[i][k].pos][3];
        }
    }
}

void
nmod_sparse_mat_mul_m4ri_mat(mzd_t* C, const nmod_sparse_mat_t A, const mzd_t *B)
{
    int i, k;
    mzd_set_ui(C, 0);
    if (B->width == 1)
    {
        nmod_sparse_mat_mul_m4ri_mat_w1(C, A, B);
        return;
    }
    else if (B->width == 2)
    {
        nmod_sparse_mat_mul_m4ri_mat_w2(C, A, B);
        return;
    }
    else if (B->width == 4)
    {
        nmod_sparse_mat_mul_m4ri_mat_w4(C, A, B);
        return;
    }

    for (i = 0; i < A->r; i++)
    {
        for (k = 0; k < A->row_supports[i]; k++)
        {
            mzd_combine_even_in_place(C, i, 0, B, A->rows[i][k].pos, 0);
        }
    }
}

struct pair
{
    int key;
    int value;
};

int compare(const void* a, const void* b)
{
     int ka = ( (struct pair *) a)->key;
     int kb = ( (struct pair *) b)->key;

     return (ka - kb);
}

void
ALGO1(mzd_t ** P, mzd_t * Xe, int * delta)
{
    mzd_t * XeT, * PT, *tmp;
    int i, j;
    int m = Xe->nrows;
    int n = Xe->ncols - Xe->nrows;

    tmp = mzd_transpose(NULL, Xe);
    XeT = mzd_init(Xe->ncols, Xe->nrows);
    PT = mzd_init(m + n, m + n);

    struct pair Delta[m + n];
    for (i = 0; i < m + n; i++)
    {
        Delta[i].key = delta[i];
        Delta[i].value = i;
    }

    qsort(Delta, m + n, sizeof(struct pair), compare);

    mzd_set_ui(PT, 0);
    for (i = 0; i < m + n; i++)
        mzd_write_bit(PT, i, Delta[i].value, 1);
    mzd_mul(XeT, PT, tmp, MZD_MUL_CUTOFF);                       /*!!!!!!! MISSING IN PAPER */
    mzd_free(tmp);

    for (i = 0; i < m + n; i++)
        delta[i] = Delta[i].key;

    int busy[m + n];
    for (i = 0; i < m + n; i++)
        busy[i] = 0;

    for (i = 0; i < m; i++)
    {
        int j0;
        for (j0 = 0; j0 < m + n; j0++) /* find a pivot */
            if (mzd_read_bit(XeT, j0, i) && !busy[j0])
                break;
        if (j0 < m + n)                                                 //MISSING IN PAPER??
        {
            busy[j0] = 1;
            for (j = j0 + 1; j < m + n; j++)
            {
                if (mzd_read_bit(XeT, j, i))
                {
                    mzd_row_add(XeT, j0, j);
                    mzd_row_add(PT, j0, j);
                }
            }
        }
        else {
            assert(0);
        }
    }

    mzd_transpose(Xe, XeT);
    mzd_transpose(P[0], PT);
    mzd_set_ui(P[1], 0);

    for (j = 0; j < m + n; j++)
    {
        if (busy[j])
        {
            /* TODO do this better, with mul and neg? */
            for (i = 0; i < P[0]->nrows; i++)
            {
                mzd_write_bit(P[1], i, j, mzd_read_bit(P[0], i, j));
                mzd_write_bit(P[0], i, j, 0);
            }
        }
    }
    mzd_free(XeT);
    mzd_free(PT);
}

/*
slong
berlekamp_massey(mzd_t ** qi, mzd_t **mty, const slong size)
{
    return 0;
}
*/

/*
 * This is an implementation of Coppersmith's block Wiedemann algorithm.
 * It follows the article "Fast computation of linear generators for matrix sequences
 * and application to the block Wiedemann algorithm" by Emmanuel Thome.
 */
void
_bw(mzd_t *K, const nmod_sparse_mat_t M, const int skip, const int epsilon, const int m_shift, const int n_shift)
{
    mzd_t *e, *x, *id, *zero, **mtz, **mty, **xmty, **P, **f;
    const slong N = M->r, m = 2 << m_shift, n = 2 << n_shift;
    slong i, j, count = 0;
    int done = 0, ftries, t = (m + n - 1)/n + skip, delta[m + n];
    int L = (N+m-1)/m + N/n + epsilon, max_diff;
    /* TODO check the correct things are const */

    /* set up dense matrices */
    x = mzd_init(m, N);
    mtz = (mzd_t **)malloc((L + 1) * sizeof(mzd_t *));
    mty = &mtz[1];
    xmty = (mzd_t **)malloc(L * sizeof(mzd_t *));
    P = (mzd_t **)malloc(2 * sizeof(mzd_t *)); /* a pair of matrices */
    f = (mzd_t **)malloc((t + 1) * sizeof(mzd_t *));
    P[0] = mzd_init(m + n, m + n);
    P[1] = mzd_init(m + n, m + n);
    mtz[0] = mzd_init(N, n);

    for (i = 0; i < L; i++)
    {
        mty[i] = mzd_init(N, n);
        xmty[i] = mzd_init(m, n);
    }

    for (i = 0; i <= t; i++)
        f[i] = mzd_init(n, m + n);

    for (i = 0; i < n + m; i++)
        delta[i] = t;

    zero = mzd_init(n, m);
    id = mzd_init(n, n);
    mzd_set_ui(id, 1);
    mzd_concat(f[t], zero, id);
    mzd_free(id);
    mzd_free(zero);

    e = mzd_init(m, m + n);
    while (!done)
    {
        /* random initial y = Mz = mty[0] and x */
        mzd_randomize(mtz[0]);
        mzd_randomize(x);

        for (i = 0; i < L; i++)
        {
            nmod_sparse_mat_mul_m4ri_mat(mty[i], M, mtz[i]);
            mzd_mul(xmty[i], x, mty[i], MZD_MUL_CUTOFF);
        }

        for (ftries = 0; ftries <= 10 && !done; ftries++)
        {
            mzd_mul(e, xmty[0], f[t], MZD_MUL_CUTOFF);
            mzd_t * m_cols;
            for (i = 0; i < t; i++)
            {
                /*ithCoefficientOfVectorPolynomialProduct(B.getField(),m,m+n,AX,f,t);*/
                mzd_randomize(f[i]);
                for (j = 0; j < f[i]->nrows; j++)
                    mzd_row_clear_offset(f[i], j, m);
                mzd_addmul(e, xmty[t - i], f[i], MZD_MUL_CUTOFF);
            }


            m_cols = mzd_submatrix(NULL, e, 0, 0, m, m);
            done = (mzd_echelonize(m_cols, 0) == m);
            mzd_free(m_cols);
        }
    }

    for (max_diff = 0; max_diff <= N/m; t++, max_diff++)
    {
        mzd_t * f_tmp;

        /*mzd_t * e = ithCoefficientOfVectorPolynomialProduct(B.getField(),m,m+n,AX,f,t); */

        if (done) /* first run only */
        {
            done = 0;
        }
        else
        {
            mzd_set_ui(e, 0);
            for (i = 0; i <= t; i++)
                mzd_addmul(e, xmty[t - i], f[i], MZD_MUL_CUTOFF);
        }

        ALGO1(P, e, delta);

        /* f = productWithLinear(B.getField(), n, m + n, f, P);*/
        f = (mzd_t **)realloc(f, (t + 2) * sizeof(mzd_t *));
        f[t + 1] = mzd_init(n, m + n);
        f_tmp = mzd_init(n, m + n);
        mzd_mul(f[t + 1], f[t], P[1], MZD_MUL_CUTOFF);
        for (i = t; i >= 1; i--)
        {
            mzd_mul(f_tmp, f[i], P[0], MZD_MUL_CUTOFF);
            mzd_copy(f[i], f_tmp);
            mzd_addmul(f[i], f[i - 1], P[1], MZD_MUL_CUTOFF);
        }
        mzd_mul(f_tmp, f[0], P[0], MZD_MUL_CUTOFF);
        mzd_copy(f[0], f_tmp);

        max_diff = 0;
        for (i = 0; i < m + n; i++)
        {
            for (j = 0; j < m + n; j++)
            {
                if (mzd_read_bit(P[1], j, i))
                {
                    delta[i]++;
                    break;
                }
            }

            if (t - delta[i] > max_diff)
                max_diff = t - delta[i];
        }

        mzd_free(f_tmp);
    }

    /* extract kernel vector */
    printf("fin\n");

    mzd_t *w, *wT, *Mw, *f_w, ** mtzT, **fT;
    w = mzd_init(N, 1);
    wT = mzd_init(1, N);
    Mw = mzd_init(N, 1);
    mtzT = (mzd_t **)malloc((L + 1) * sizeof(mzd_t *));
    fT = (mzd_t **)malloc((t + 1) * sizeof(mzd_t *));
    for (j = 0; j < L + 1; j++)
        mtzT[j] = mzd_transpose(NULL, mtz[j]);
    for (j = 0; j < t + 1; j++)
        fT[j] = mzd_transpose(NULL, f[j]);
    for (j = 0; j < m + n; j++)
    {
        /*printf("%ld, %d, %d\n", j, delta[j], t- delta[j]);*/
        mzd_set_ui(wT, 0);
        for (i = 0; i <= delta[j]; i++)
        {
            f_w = mzd_init_window(fT[i], j, 0, j + 1, f[i]->nrows);
            /*printf("%ld, %d\n", delta[j] - i, L);*/
            _mzd_mul_va(wT, f_w , mtzT[delta[j] - i], 0);
            mzd_free_window(f_w);
        }

        /*if (!w.transposed()[0].isZero() && (M*w).transposed()[0].isZero())*/
        if (!mzd_is_zero(wT))
        {
            mzd_transpose(w, wT);
            nmod_sparse_mat_mul_m4ri_mat(Mw, M, w); /* TODO can we use existing knowledge here */
            if (mzd_is_zero(Mw))
            {
                count++;
                for (i = 0; i < w->nrows; i++)
                {
                    mzd_write_bit(K, i, 0, mzd_read_bit(w, i, 0));
                }
            }
        }
    }
    printf("found %ld kernel vecs:\n", count);
    for (j = 0; j < t + 1; j++)
        mzd_free(fT[j]);
    for (j = 0; j < L + 1; j++)
        mzd_free(mtzT[j]);
    free(mtzT);
    free(fT);
    mzd_free(w);
    mzd_free(Mw);

    /* cleanup */
    for (i = 0; i <= t; i++)
        mzd_free(f[i]);

    mzd_free(mtz[0]);
    for (i = 0; i < L; i++)
    {
        mzd_free(mty[i]);
        mzd_free(xmty[i]);
    }

    mzd_free(e);
    mzd_free(x);
    mzd_free(P[0]);
    mzd_free(P[1]);
    free(mtz);
    free(xmty);
    free(P);
    free(f);
}

void
bw(mzd_t *K, const nmod_sparse_mat_t M)
{
    _bw(K, M, 1, 1, 8, 8);
}
