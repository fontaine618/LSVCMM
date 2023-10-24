#include "RcppArmadillo.h"
#include "uint.h"

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef Control_hpp
#define Control_hpp

class Control {

public:

  uint max_rounds = 20;
  uint max_iter = 100;
  double rel_tol= 1e-6;
  uint verbose = 0;
  std::string update_method = "PGD";
  double backtracking_fraction = 0.9;
  bool two_step_estimation = true;
  double stepsize_factor = 0.1;

  Control(
    uint max_rounds,
    uint max_iter,
    double rel_tol,
    uint verbose,
    std::string update_method,
    double backtracking_fraction,
    bool two_step_estimation,
    double stepsize_factor
  ){
    this->max_iter = max_iter;
    this->max_rounds = max_rounds;
    this->rel_tol = rel_tol;
    this->verbose = verbose;
    this->update_method = update_method;
    this->backtracking_fraction = backtracking_fraction;
    this->two_step_estimation = two_step_estimation;
    this->stepsize_factor = stepsize_factor;
  }

  Control(){}

  void add_to_results(Rcpp::List& results){
    results["control.max_iter"] = this->max_iter;
    results["control.max_rounds"] = this->max_rounds;
    results["control.rel_tol"] = this->rel_tol;
    results["control.verbose"] = this->verbose;
    results["control.update_method"] = this->update_method;
    results["control.backtracking_fraction"] = this->backtracking_fraction;
    results["control.two_step_estimation"] = this->two_step_estimation;
    results["control.stepsize_factor"] = this->stepsize_factor;
  }

};

#endif
