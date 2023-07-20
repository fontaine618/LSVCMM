#include "RcppArmadillo.h"
#include "Kernel.hpp"

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef GaussianKernel_hpp
#define GaussianKernel_hpp

class GaussianKernel : public Kernel{

public:

  GaussianKernel();
  arma::mat eval(const arma::rowvec &t0, const arma::colvec &t1);

  void add_to_results(Rcpp::List& results);
};

#endif
