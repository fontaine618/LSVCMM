#include "RcppArmadillo.h"
#include "Family.hpp"

//[[Rcpp::depends(RcppArmadillo)]]


Family::Family(){
  this->dispersion = 1.;
}

std::vector<arma::colvec> Family::unit_variance(const std::vector<arma::colvec> &mean){
  std::vector<arma::colvec> out(mean.size());
  for(uint i=0; i<mean.size(); i++) out[i] = this->unit_variance(mean[i]);
  return out;
}

arma::colvec Family::unit_variance(const arma::colvec &mean){
  return arma::ones(mean.size());
}

void Family::add_to_results(Rcpp::List& results){}

void Gaussian::add_to_results(Rcpp::List& results){
  results["family.dispersion"] = this->dispersion;
  results["family.name"] = "gaussian";
}
