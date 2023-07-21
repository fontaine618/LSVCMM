#include "RcppArmadillo.h"
#include "Data.hpp"
#include "Kernel.hpp"
#include "LinkFunction.cpp"
#include "Family.hpp"
#include "WorkingCovariance.hpp"
#include "Penalty.hpp"
#include "Control.cpp"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef Model_hpp
#define Model_hpp

class Model{

private:


public:

  uint px, pu;

  arma::mat B;
  arma::mat gB;
  double LB = 0.1;

  arma::colvec a;
  arma::colvec ga;
  double La = 0.01;

  Rcpp::List results = Rcpp::List::create();

  Kernel kernel;
  Penalty penalty;
  WorkingCovariance workingCovariance;
  LinkFunction linkFunction;
  Family family;
  Control control;
  arma::rowvec estimated_time;
  uint nt;

  Model();
  Model(
    uint px,
    uint pu,
    arma::rowvec estimated_time,
    Penalty penalty,
    WorkingCovariance workingCovariance,
    LinkFunction linkFunction,
    Family family,
    Control control,
    Kernel kernel
  );

  std::vector<arma::colvec> linear_predictor(Data &data); // linear predictor for each observation
  void update_mean(Data &data); // should be called whenever B or a is updated
  void update_precision(Data &data); // should be called whenever workingCovariance is updated

  double quadratic_term(const Data &data); // objective for mean only
  double logdet_precision(const Data &data); // update only when P[i] are updated
  double logdet_scaling(const Data &data); // update whenever mean parms were updated
  double logdet_dispersion(const Data &data); // update whenever dispersion is updated
  double quasi_log_likelihood(double quad_term, double logdet_term); // NB logdet should be scaling+dispersion-precision

  void update_gradients(const Data &data);
  arma::mat hessian(const Data &data);
  void proximal_gradient_step();
  void fit(Data &data);
  void initialize(Data &data);
  void prepare_stepsize(Data &data);
  double lambda_max(Data &data);

  void prepare_results(const Data &data);
  void prepare_ics(const Data &data);
  Rcpp::List save();


};

#endif
