#include "RcppArmadillo.h"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef Interpolator_hpp
#define Interpolator_hpp

class Interpolator{

public:

  arma::rowvec time;

  Interpolator(const arma::rowvec &time);

  arma::mat interpolator_matrix(const arma::colvec &new_time);
  std::vector<arma::mat> interpolator_matrix(const std::vector<arma::colvec> &new_time);

};

#endif
