#include "RcppArmadillo.h"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef LinkFunction_hpp
#define LinkFunction_hpp

class LinkFunction{
public:
  std::vector<arma::colvec> eval(const std::vector<arma::colvec> &linear_predictor);
  std::vector<arma::colvec> derivative(const std::vector<arma::colvec> &linear_predictor);
  static LinkFunction* Create(std::string name);

  virtual ~LinkFunction(){};
  virtual arma::colvec eval(const arma::colvec &linear_predictor) = 0;
  virtual arma::colvec derivative(const arma::colvec &linear_predictor) = 0;
  virtual void add_to_results(Rcpp::List& results) = 0;
};

class IdentityLink : public LinkFunction{
public:
  IdentityLink();
  using LinkFunction::eval;
  using LinkFunction::derivative;
  arma::colvec eval(const arma::colvec &linear_predictor);
  arma::colvec derivative(const arma::colvec &linear_predictor);
  void add_to_results(Rcpp::List& results);
};

class LogLink : public LinkFunction{
public:
  LogLink();
  using LinkFunction::eval;
  using LinkFunction::derivative;
  arma::colvec eval(const arma::colvec &linear_predictor);
  arma::colvec derivative(const arma::colvec &linear_predictor);
  void add_to_results(Rcpp::List& results);
};

#endif
