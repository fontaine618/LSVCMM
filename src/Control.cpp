#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef Control_hpp
#define Control_hpp

class Control {

public:

  uint max_rounds = 20;
  uint max_iter = 100;
  double rel_tol= 1e-6;
  uint verbose = 0;

  Control(uint max_rounds, uint max_iter, double rel_tol, uint verbose){
    this->max_iter = max_iter;
    this->max_rounds = max_rounds;
    this->rel_tol = rel_tol;
    this->verbose = verbose;
  }

  Control(){}

};

#endif
