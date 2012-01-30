#include "rbm.hpp"

#include <armadillo>
#include <assert.h>
#include <boost/container/vector.hpp>
#include <boost/tuple/tuple.hpp>
#include <gsl/gsl_qrng.h>
#include <iostream>
#include <float.h>
#include <math.h>

using namespace arma;
using namespace boost;
using namespace std;

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

vec param_to_simplex(int n, int k, mat w, vec b, vec c) {
    mat v = binary(n);
    mat h = binary(k).t();

    mat energy = 
          v * w * h
        + v * b * ones<mat>(1, 1 << k)
        + ones<mat>(1 << n, 1) * c.t() * h;
    mat psi = exp(energy);

    vec uprob = sum(psi, 1);
    vec prob = uprob / sum(uprob);

    return prob;
}


const double PRNG_RANGE = 2;

inline double uniform(double x) {
    return x * 2 * PRNG_RANGE - PRNG_RANGE;
}

/* Use quasi-random sequences from GSL */
vector<vec> sample(int n, int k, int ntrials) {
    // The gsl_qrng_halton algorithm supports sampling from a space of at most
    // 1229 dimensions.  It is also deterministic, I think.
    int nparams = n*k + n + k;
    assert(nparams < 1229);

    vector<vec> res;

    for (int t=0; t < ntrials; ++t) {
        mat w = randu<mat>(n, k) * 2 * PRNG_RANGE - PRNG_RANGE;
        mat b = randu<mat>(n) * 2 * PRNG_RANGE - PRNG_RANGE;
        mat c = randu<mat>(k) * 2 * PRNG_RANGE - PRNG_RANGE;
        res.push_back(param_to_simplex(n, k, w, b, c));
    }

    /*
    gsl_qrng* q = gsl_qrng_alloc(gsl_qrng_halton, nparams);
    double* params = new double[nparams];

    for (int t=0; t < ntrials; ++t) {
        gsl_qrng_get(q, params);
        int pidx = 0;

        mat w(n, k);
        for (int i=0; i < n; ++i)
            for (int j=0; j < k; ++j)
                w(i, j) = uniform(params[pidx++]);


        vec b(n);
        for (int i=0; i < n; ++i)
            b(i) = uniform(params[pidx++]);

        vec c(k);
        for (int j=0; j < k; ++j)
            c(j) = uniform(params[pidx++]);

        res.push_back(param_to_simplex(n, k, w, b, c));
    }

    delete[] params;
    gsl_qrng_free(q);
    */

    return res;
}

double total_variation_distance(vec x, vec y) {
    assert(x.n_elem == y.n_elem);

    double acc = 0;
    for (int i=0; i < x.n_elem; ++i)
        acc += abs(x[i] - y[i]);

    return 0.5 * acc;
}

double hausdorff(vector<vec> xs, vector<vec> ys) {

    int n = xs.size();
    int m = ys.size();

    // maximize
    double a = DBL_MIN;
    for (int i=0; i < n; ++i) {
        // minimize
        double c = DBL_MAX;
        for (int j=0; j < m; ++j) {
            double d = total_variation_distance(xs[i], ys[j]);
            c = d < c ? d : c;
        }
        a = c > a ? c : a;
    }

    // maximize
    double b = DBL_MIN;
    for (int j=0; j < m; ++j) {
        // minimize
        double c = DBL_MAX;
        for (int i=0; i < n; ++i) {
            double d = total_variation_distance(xs[i], ys[j]);
            c = d < c ? d : c;
        }
        b = c > b ? c : b;
    }

    return (a < b) ? b : a;
}
