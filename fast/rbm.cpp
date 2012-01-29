#include "rbm_mc.hpp"

#include <armadillo>
#include <assert.h>
#include <boost/container/vector.hpp>
#include <boost/tuple/tuple.hpp>
#include <gsl/gsl_qrng.h>
#include <iostream>

using namespace arma;
using namespace boost;
using namespace std;

const double PRNG_RANGE = 3;

vec param_to_simplex(int n, int k, mat w, vec b, vec c) {
    mat v = binary(n);
    mat h = binary(k).t();

    mat energy = 
          v * w * h
        + v * b * ones<mat>(1, 1 << k)
        + ones<mat>(1 << n, 1) * c.t() * h;
    mat psi = exp(energy);

    vec uprob = sum(psi, 1);  // Sum (i.e. marginalize) across hidden states.
    vec prob = uprob / sum(uprob);

    return prob;
}

/* Ex. for n = 3:   0 0 0
 *                  0 0 1
 *                  0 1 0
 *                   ...
 *                  1 1 1
 */
mat binary(int n) {
    int nrows = 1 << n;
    mat ret(nrows, n);

    for (int col=0; col < n; ++col) {
        int nsteps = 1 << (n - col - 1);
        int step = 0;
        int val = 0;
        for (int row=0; row < nrows; ++row) {
            ret(row, col) = val;
            if (++step == nsteps) {
                val = val ? 0 : 1;
                step = 0;
            }
        }
    }
    return ret;
}

/* Use quasi-random sequences from GSL */
vector<vec> sample(int n, int k, int ntrials) {
    // The gsl_qrng_halton algorithm supports sampling from a space of at most
    // 1229 dimensions.  It is also deterministic, I think.
    int nparams = n*k + n + k;
    assert(nparams < 1229);

    vector<vec> res;

    gsl_qrng* q = gsl_qrng_alloc(gsl_qrng_halton, nparams);
    double* params = new double[nparams];


    for (int t=0; t < ntrials; ++t) {
        gsl_qrng_get(q, params);
        int pidx = 0;

        mat w(n, k);
        for (int i=0; i < n; ++i)
            for (int j=0; j < k; ++j)
                w(i, j) = params[pidx++];


        vec b(n);
        for (int i=0; i < n; ++i)
            b(i) = params[pidx++];

        vec c(k);
        for (int j=0; j < k; ++j)
            c(j) = params[pidx++];

        res.push_back(param_to_simplex(n, k, w, b, c));
    }

    delete[] params;
    gsl_qrng_free(q);

    return res;
}

/* Below here are unused things. */

void test() {
}


tuple<mat, vec, vec> rand_param(int n, int k) {
    mat w = (2 * PRNG_RANGE * randu<mat>(n, k)) + PRNG_RANGE;
    vec b = (2 * PRNG_RANGE * randu<vec>(n)) + PRNG_RANGE;
    vec c = (2 * PRNG_RANGE * randu<vec>(k)) + PRNG_RANGE;
}
