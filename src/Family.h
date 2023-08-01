#include "RcppArmadillo.h"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef Family_hpp
#define Family_hpp

class Gaussian{

public:

  double dispersion;

  Gaussian();

  arma::colvec unit_variance(const arma::colvec &mean);

  std::vector<arma::colvec> unit_variance(const std::vector<arma::colvec> &mean);

  double log_likelihood(const arma::colvec &residuals, const arma::mat &Precision);

  void add_to_results(Rcpp::List& results);

};

#endif
