#include "RcppArmadillo.h"
#include "LinkFunction.h"
#include "uint.h"

// [[Rcpp::depends(RcppArmadillo)]]


// =============================================================================
// Interface + Factory
LinkFunction* LinkFunction::Create(std::string name){
  if(name=="identity") return new IdentityLink();
  if(name=="log") return new LogLink();
  Rcpp::stop("Unknown family: %s", name);
}

std::vector<arma::colvec> LinkFunction::eval(const std::vector<arma::colvec> &linear_predictor){
  std::vector<arma::colvec> out(linear_predictor.size());
  for(uint i=0; i<linear_predictor.size(); i++) out[i] = this->eval(linear_predictor[i]);
  return out;
}

std::vector<arma::colvec> LinkFunction::derivative(const std::vector<arma::colvec> &linear_predictor){
  std::vector<arma::colvec> out(linear_predictor.size());
  for(uint i=0; i<linear_predictor.size(); i++) out[i] = this->derivative(linear_predictor[i]);
  return out;
}

// =============================================================================
// Identity
IdentityLink::IdentityLink(){}

arma::colvec IdentityLink::eval(const arma::colvec &linear_predictor){
  return linear_predictor;
}

arma::colvec IdentityLink::derivative(const arma::colvec &linear_predictor){
  return arma::ones(linear_predictor.size());
}

void IdentityLink::add_to_results(Rcpp::List& results){
  results["link_function.name"] = "identity";
}

// =============================================================================
// Log
LogLink::LogLink(){}

arma::colvec LogLink::eval(const arma::colvec &linear_predictor){
  return arma::exp(linear_predictor);
}

arma::colvec LogLink::derivative(const arma::colvec &linear_predictor){
  return arma::exp(linear_predictor);
}

void LogLink::add_to_results(Rcpp::List& results){
  results["link_function.name"] = "log";
}


