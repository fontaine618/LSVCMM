#include "RcppArmadillo.h"
#include "uint.h"

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef LinkFunction_hpp
#define LinkFunction_hpp


class Identity {

public:
  Identity(){}

  arma::colvec eval(const arma::colvec &linear_predictor){
    return linear_predictor;
  }

  arma::colvec derivative(const arma::colvec &linear_predictor){
    return arma::ones(linear_predictor.size());
  }

  void add_to_results(Rcpp::List& results){
    results["link_function.name"] = "identity";
  }

  // overload for lists
  std::vector<arma::colvec> eval(const std::vector<arma::colvec> &linear_predictor){
    std::vector<arma::colvec> out(linear_predictor.size());
    for(uint i=0; i<linear_predictor.size(); i++) out[i] = this->eval(linear_predictor[i]);
    return out;
  }

  std::vector<arma::colvec> derivative(const std::vector<arma::colvec> &linear_predictor){
    std::vector<arma::colvec> out(linear_predictor.size());
    for(uint i=0; i<linear_predictor.size(); i++) out[i] = this->derivative(linear_predictor[i]);
    return out;
  }

};

#endif
