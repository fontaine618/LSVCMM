#include "RcppArmadillo.h"
#include "Logger.cpp"
#include "Control.cpp"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef WorkingCovariance_hpp
#define WorkingCovariance_hpp


class CompoundSymmetry{

public:

  double variance_ratio = 1.0;
  bool estimate_parameters = TRUE;

  // TODO: perhaps store the Z[i] in data and define a function here that computes it

  CompoundSymmetry();
  CompoundSymmetry(double variance_ratio, bool estimate_parameters);

  arma::mat compute_covariance(const arma::colvec &time);
  std::vector<arma::mat> compute_covariance(const std::vector<arma::colvec> &time);

  arma::mat compute_precision(const arma::colvec &time);
  std::vector<arma::mat> compute_precision(const std::vector<arma::colvec> &time);

  uint update_parameters(
      const std::vector<arma::colvec> &sr,
      const std::vector<arma::colvec> &t,
      const std::vector<arma::mat> &P,
      const double dispersion,
      Logger *logger,
      uint round,
      Control *control
  );

  double profile_likelihood(
      const std::vector<arma::colvec> &sr,
      const std::vector<arma::mat> &P
  );

  std::vector<double> derivatives(
      const std::vector<arma::colvec> &sr,
      const std::vector<arma::mat> &P,
      const double dispersion
  );

  void add_to_results(Rcpp::List& results);

};

#endif
