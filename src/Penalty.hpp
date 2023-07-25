#include "RcppArmadillo.h"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef Penalty_hpp
#define Penalty_hpp

class Penalty{

public:

  double lambda = 0.;
  double alpha = 1.;
  double power = 0.;
  bool penalize_intercept = FALSE;

  arma::mat W1;
  arma::colvec w2;

  Penalty(
    double lambda,
    double alpha,
    double power,
    bool penalize_intercept
  );

  Penalty();

  arma::mat proximal(const arma::mat &B, double stepsize);

  arma::mat proximal_l1(const arma::mat &B, arma::mat M);

  arma::mat proximal_l2(const arma::mat &B, arma::colvec m);

  void update_weights(const arma::mat &B);

  void unit_weights(const arma::mat &B);

  void add_to_results(Rcpp::List& results);

  double lambda_max(arma::mat B, const arma::mat &gB, const double stepsize);

};

#endif
