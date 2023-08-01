#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef GaussianKernel_hpp
#define GaussianKernel_hpp

// class GaussianKernel : public Kernel{
class GaussianKernel{

public:

  double scale, eval0;
  GaussianKernel();
  // using Kernel::eval;
  arma::mat eval(const arma::rowvec &t0, const arma::colvec &t1);
  std::vector<arma::mat> eval(const std::vector<arma::rowvec> &t0, const arma::colvec &t1);
  std::vector<arma::mat> eval(const arma::rowvec &t0, const std::vector<arma::colvec> &t1);

  void add_to_results(Rcpp::List& results);
};

#endif
