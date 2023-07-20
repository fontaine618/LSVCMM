#include "RcppArmadillo.h"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef Kernel_hpp
#define Kernel_hpp

class Kernel{

public:
  double scale, eval0;

  Kernel();

  arma::mat eval(const arma::rowvec &t0, const arma::colvec &t1);

  std::vector<arma::mat> eval(const std::vector<arma::rowvec> &t0, const arma::colvec &t1);
  std::vector<arma::mat> eval(const arma::rowvec &t0, const std::vector<arma::colvec> &t1);

  void add_to_results(Rcpp::List& results);

};

class UniformKernel : public Kernel{};

#endif
