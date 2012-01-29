#ifndef RBM_MC_H
#define RBM_MC_H

#include <armadillo>
#include <boost/container/vector.hpp>
#include <gsl/gsl_qrng.h>

arma::vec param_to_simplex(int n, int k, arma::mat w, arma::vec b, arma::vec c);
std::vector<arma::vec> sample(int n, int k, int ntrials);
arma::mat binary(int);

#endif
