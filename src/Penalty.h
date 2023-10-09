#include "RcppArmadillo.h"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]

#ifndef Penalty_hpp
#define Penalty_hpp

class Penalty{
  // Base class assumes a weighted sparse group lasso structure
  //    sum_s alpha * w_j *  p_lambda(|b_j|) + (1-alpha) * w * p_{lambda*sqrt(T)}(||b_j||)
  // Variants use different weights
  // - p'_lambda(|B|) for the entrywise penalty
  // - p'_{lambda*sqrt(T)}(||B||) for the group penalty
  // Usage:
  // - Create a penalty object
  // - Provide unpenalized solution in update_B0
  // - Whenever a new lambda is set, need to update the weights (not necessary for ASGL, but it's cheap)
  // - Call proximal to get the update given the gradient and stepsize

public:
  static Penalty* Create(std::string name, double lambda, double alpha, double power, double a, bool penalize_intercept);
  double lambda = 0.;
  double alpha = 1.;
  bool compute_penalty_weights = true;
  bool penalize_intercept = FALSE;
  arma::mat W1;
  arma::colvec w2;
  arma::mat AbsB0; // Value of |B| without penalty
  arma::colvec normB0; // row norm of B0
  arma::mat proximal(const arma::mat &B, double stepsize);
  double eval(const arma::mat &B);
  void update_B0(const arma::mat &B); // should be called only once, once we get an unpenalized estimate
  double lambda_max(arma::mat B, const arma::mat &gB, const double stepsize);
  void update_weights(); // should be called any time lambda or B0 is updated
  void unit_weights(const arma::mat &B);
  void large_weights(const arma::mat &B);

  virtual ~Penalty(){};
  virtual void add_to_results(Rcpp::List& results) = 0;
  virtual arma::mat derivative(const arma::mat &B, double scaling) = 0;
};

class AdaptiveSparseGroupLasso: public Penalty{
  // Weights given by p'_lambda(|B|)=lambda*|B|^-power
  // and lambda*sqrt(T)||B||^-power

public:

  double power = 0.5;
  AdaptiveSparseGroupLasso(
    double lambda,
    double alpha,
    double power,
    bool penalize_intercept
  );
  void add_to_results(Rcpp::List& results);
  arma::mat derivative(const arma::mat &B, double scaling);
};

class SparseGroupSCAD: public Penalty{
  // Weights given by p'_lambda(|B|) for the SCAD penalty
public:
  double a = 3.7;
  SparseGroupSCAD(
    double lambda,
    double alpha,
    double a,
    bool penalize_intercept
  );
  void add_to_results(Rcpp::List& results);
  arma::mat derivative(const arma::mat &B, double scaling);
};

class SparseGroupMCP: public Penalty{
  // Weights given by p'_lambda(|B|) for the SCAD penalty
public:
  double a = 3.7;
  SparseGroupMCP(
    double lambda,
    double alpha,
    double a,
    bool penalize_intercept
  );
  void add_to_results(Rcpp::List& results);
  arma::mat derivative(const arma::mat &B, double scaling);
};

#endif
