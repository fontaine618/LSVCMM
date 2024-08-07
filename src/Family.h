#include "RcppArmadillo.h"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef Family_hpp
#define Family_hpp

class Family{
public:
  double dispersion, power;
  std::vector<arma::colvec> unit_variance(const std::vector<arma::colvec> &mean);
  static Family* Create(std::string name, double power);
  virtual ~Family(){};

  virtual arma::colvec unit_variance(const arma::colvec &mean) = 0;
  virtual void add_to_results(Rcpp::List& results) = 0;
};

class Gaussian : public Family{
public:
  Gaussian();
  using Family::unit_variance;
  arma::colvec unit_variance(const arma::colvec &mean);
  void add_to_results(Rcpp::List& results);
};

class Tweedie : public Family{
public:
  Tweedie(double power);
  using Family::unit_variance;
  arma::colvec unit_variance(const arma::colvec &mean);
  void add_to_results(Rcpp::List& results);
};

#endif
