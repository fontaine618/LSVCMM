#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef LinkFunction_hpp
#define LinkFunction_hpp

class LinkFunction {

public:

  LinkFunction(){}

  arma::colvec eval(const arma::colvec &linear_predictor){
    return linear_predictor;
  }

  // maps linear predictor to mean
  std::vector<arma::colvec> eval(const std::vector<arma::colvec> &linear_predictor){
    std::vector<arma::colvec> out(linear_predictor.size());
    for(uint i=0; i<linear_predictor.size(); i++) out[i] = this->eval(linear_predictor[i]);
    return out;
  }

  arma::colvec derivative(const arma::colvec &linear_predictor){
    return arma::ones(linear_predictor.size());
  }

  std::vector<arma::colvec> derivative(const std::vector<arma::colvec> &linear_predictor){
    std::vector<arma::colvec> out(linear_predictor.size());
    for(uint i=0; i<linear_predictor.size(); i++) out[i] = this->derivative(linear_predictor[i]);
    return out;
  }

  void add_to_results(Rcpp::List& results){}

};

class Identity : public LinkFunction{

  void add_to_results(Rcpp::List& results){
    results["link_function.name"] = "identity";
  }};

#endif
