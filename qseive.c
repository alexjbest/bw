#include <stdio.h>
#include <gmp.h>
#include "bw.h"
#include "fmpz.h"
#include "ulong_extras.h"

void
_qseive(const mp_limb_t n, const mp_limb_t B)
{
    nmod_sparse_mat_t M;
    mp_limb_t quad, *quads, *xs, x, i = 0, j, piB = n_prime_pi(B);
    const mp_limb_t * ps = n_primes_arr_readonly(piB + 2);
    const double * pinvs = n_prime_inverses_arr_readonly(piB + 2);
    mzd_t *K;

    /* init */
    quads = (mp_limb_t *)malloc((piB + 1)*sizeof(mp_limb_t *));
    xs = (mp_limb_t *)malloc((piB + 1)*sizeof(mp_limb_t *));
    K = mzd_init(piB + 1, 1);
    nmod_sparse_mat_init(M, piB + 1, piB + 1, 2);

    printf("init done\n");
    printf("using %ld primes\n", piB);
    /* seive */

    for (x = n_sqrt(n), i = 0; i <= piB; x++)
    {
        quad = x*x - n;
        if (quad == 0)
            continue;
        for (j = 0; j < piB; j++)
            n_remove2_precomp(&quad, ps[j], pinvs[j]);
        if (quad == 1) /* was B-smooth */
        {
            quads[i] = x*x - n;
            quad = x*x - n;
            for (j = 0; j < piB; j++)
            {
                if (n_remove2_precomp(&quad, ps[j], pinvs[j]) % 2)
                    _nmod_sparse_mat_set_entry(M, j, i, M->row_supports[j], 1);
            }
            xs[i] = x;
            i++;
        }
    }
    printf("data collection done\n");

    n_cleanup_primes();

    _bw(K, M, 1, 2, 7, 7);

    printf("procesing complete\n");
    mzd_print(K);

    int done = 0;
    for (j = 0; !done; j++)
    {
        fmpz_t a, b, diff, N;
        fmpz_init_set_ui(a, 1);
        fmpz_init_set_ui(b, 1);
        fmpz_init_set_ui(N, n);
        fmpz_init(diff);
        for (i = 0; i < piB; i++)
        {
            if (mzd_read_bit(K, i, j))
            {
                fmpz_mul_ui(a, a, xs[i]);
                fmpz_mul_ui(b, b, quads[i]);
            }
        }
        assert(fmpz_is_square(b));
        fmpz_sqrt(b, b);
        if (fmpz_mod_ui(a, a, n) != fmpz_mod_ui(b, b, n) && fmpz_mod_ui(a, a, n) != n - fmpz_mod_ui(b, b, n))
        {
            done = 1;

            fmpz_print(a);
            printf("\n");
            fmpz_print(b);
            printf("\n");
            fmpz_sub(diff, a, b);
            fmpz_gcd(a, diff, N);
            fmpz_divexact(b, N, a);
            fmpz_print(a);
            printf("\n");
            fmpz_print(b);
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(N);
        fmpz_clear(diff);
    }
    /* cleanup */
    free(quads);
    free(xs);

    mzd_free(K);

    nmod_sparse_mat_clear(M);

    return;
}

void
qseive(mp_limb_t n)
{
    double logn = log(n), loglogn = log(logn);
    mp_limb_t B = (mp_limb_t)(exp((0.5f + 0.38f)*sqrt(logn*loglogn)));

    _qseive(n, B);
}

int
main(int argc, char * argv[])
{
    mp_limb_t n = 10142789312725007;
    n = n_nextprime(1000000000,0)*n_nextprime(1200000000,0);
    qseive(n);
}
