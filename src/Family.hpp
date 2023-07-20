#include "RcppArmadillo.h"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef Family_hpp
#define Family_hpp

class Family{

public:

  double dispersion;

  Family();

  arma::colvec unit_variance(const arma::colvec &mean);

  std::vector<arma::colvec> unit_variance(const std::vector<arma::colvec> &mean);

  double log_likelihood(const arma::colvec &residuals, const arma::mat &Precision);

  void add_to_results(Rcpp::List& results);

};

// Default is Gaussian
class Gaussian : public Family{

  void add_to_results(Rcpp::List& results);
};

#endif
