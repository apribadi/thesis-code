#ifndef RBM_H
#define RBM_H

#include <armadillo>
#include <boost/container/vector.hpp>
#include <gsl/gsl_qrng.h>

arma::vec param_to_simplex(int n, int k, arma::mat w, arma::vec b, arma::vec c);
std::vector<arma::vec> sample(int n, int k, int ntrials);
std::vector<arma::vec> sample_lattice(int n, int k);
arma::mat binary(int);
double total_variation_distance(arma::vec, arma::vec);
double hausdorff(std::vector<arma::vec>, std::vector<arma::vec>);

#endif
