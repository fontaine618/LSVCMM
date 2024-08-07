#include "RcppArmadillo.h"
#include "Family.h"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]


// =============================================================================
// Interface + Factory
Family* Family::Create(std::string name, double power){
  if(name=="gaussian") return new Gaussian();
  if(name=="tweedie") return new Tweedie(power);
  Rcpp::stop("Unknown family: %s", name);
}

std::vector<arma::colvec> Family::unit_variance(const std::vector<arma::colvec> &mean){
  std::vector<arma::colvec> out(mean.size());
  for(uint i=0; i<mean.size(); i++) out[i] = this->unit_variance(mean[i]);
  return out;
}


// =============================================================================
Gaussian::Gaussian(){
  this->dispersion = 1.;
}

arma::colvec Gaussian::unit_variance(const arma::colvec &mean){
  return arma::ones(mean.size());
}

void Gaussian::add_to_results(Rcpp::List& results){
  results["family.dispersion"] = this->dispersion;
  results["family.name"] = "gaussian";
}



// =============================================================================
Tweedie::Tweedie(double power){
  this->dispersion = 1.;
  this->power = power;
}

arma::colvec Tweedie::unit_variance(const arma::colvec &mean){
  return arma::pow(mean, this->power);
}

void Tweedie::add_to_results(Rcpp::List& results){
  results["family.dispersion"] = this->dispersion;
  results["family.name"] = "tweedie";
  results["family.power"] = this->power;
}
