#include "RcppArmadillo.h"
#include "Family.hpp"

//[[Rcpp::depends(RcppArmadillo)]]


Gaussian::Gaussian(){
  this->dispersion = 1.;
}

std::vector<arma::colvec> Gaussian::unit_variance(const std::vector<arma::colvec> &mean){
  std::vector<arma::colvec> out(mean.size());
  for(uint i=0; i<mean.size(); i++) out[i] = this->unit_variance(mean[i]);
  return out;
}

arma::colvec Gaussian::unit_variance(const arma::colvec &mean){
  return arma::ones(mean.size());
}

void Gaussian::add_to_results(Rcpp::List& results){
  results["family.dispersion"] = this->dispersion;
  results["family.name"] = "gaussian";
}
