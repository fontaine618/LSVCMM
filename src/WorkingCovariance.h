#include "RcppArmadillo.h"
#include "Logger.cpp"
#include "Control.cpp"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef WorkingCovariance_hpp
#define WorkingCovariance_hpp


class WorkingCovariance{
public:
  double correlation = 0.5;
  double variance_ratio = 1.0;
  bool estimate_parameters = TRUE;
  virtual ~WorkingCovariance(){};
  virtual arma::mat compute_precision(const arma::colvec &time) = 0;
  std::vector<arma::mat> compute_precision(const std::vector<arma::colvec> &time);
  virtual uint update_parameters(
      const std::vector<arma::colvec> &sr,
      const std::vector<arma::colvec> &t,
      const std::vector<arma::mat> &P,
      const double dispersion,
      Logger *logger,
      uint round,
      Control *control
  ) = 0;
  double profile_likelihood(
      const std::vector<arma::colvec> &sr,
      const std::vector<arma::mat> &P
  );
  virtual std::vector<double> derivatives(
      const std::vector<arma::colvec> &sr,
      const std::vector<arma::colvec> &t,
      const std::vector<arma::mat> &P,
      const double dispersion
  ) = 0;
  virtual void add_to_results(Rcpp::List& results) = 0;
  static WorkingCovariance* Create(std::string name, double variance_ratio, double correlation, bool estimate_parameters);
};

class CompoundSymmetry : public WorkingCovariance{
public:
  CompoundSymmetry();
  CompoundSymmetry(double variance_ratio, double correlation, bool estimate_parameters);
  using WorkingCovariance::compute_precision;
  arma::mat compute_precision(const arma::colvec &time);
  uint update_parameters(
      const std::vector<arma::colvec> &sr,
      const std::vector<arma::colvec> &t,
      const std::vector<arma::mat> &P,
      const double dispersion,
      Logger *logger,
      uint round,
      Control *control
  );
  std::vector<double> derivatives(
      const std::vector<arma::colvec> &sr,
      const std::vector<arma::colvec> &t,
      const std::vector<arma::mat> &P,
      const double dispersion
  );
  void add_to_results(Rcpp::List& results);
};

#endif
