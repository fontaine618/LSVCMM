#include "RcppArmadillo.h"
#include "Kernel.h"
#include "uint.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Kernel
std::vector<arma::mat> Kernel::eval(const std::vector<arma::rowvec> &t0, const arma::colvec &t1){
  std::vector<arma::mat> out(t0.size());
  for(uint i=0; i<t0.size(); i++) out[i] = this->eval(t0[i], t1);
  return out;
}

std::vector<arma::mat> Kernel::eval(const arma::rowvec &t0, const std::vector<arma::colvec> &t1){
  std::vector<arma::mat> out(t1.size());
  for(uint i=0; i<t1.size(); i++) out[i] = this->eval(t0, t1[i]);
  return out;
}

Kernel* Kernel::Create(std::string name){
  if(name=="gaussian") return new GaussianKernel();
  if(name=="epa") return new EpaKernel();
  Rcpp::stop("unrecognized kernel name: " + name);
}

// Gaussian Kernel
GaussianKernel::GaussianKernel(){
  this->scale = 1.;
  this->eval0 = 1. / sqrt(arma::datum::pi);
}

arma::mat GaussianKernel::eval(const arma::rowvec &t0, const arma::colvec &t1){
  arma::mat od(t1.n_rows, t0.n_cols);
  od.each_col() = t1;
  od.each_row() -= t0;
  od = arma::exp(- arma::square(od / this->scale) ) / (sqrt(arma::datum::pi) * this->scale);
  return od;
}

void GaussianKernel::add_to_results(Rcpp::List& results){
  results["kernel.name"] = "gaussian";
  results["kernel.scale"] = this->scale;
}

// Epanechnikov Kernel
EpaKernel::EpaKernel(){
  this->scale = 1.;
  this->eval0 = 0.75;
}

arma::mat EpaKernel::eval(const arma::rowvec &t0, const arma::colvec &t1){
  arma::mat od(t1.n_rows, t0.n_cols);
  od.each_col() = t1;
  od.each_row() -= t0;
  od = 0.75 * (1. - arma::square(od / this->scale) ) / this->scale;
  od.elem(arma::find(od<0)).zeros();
  return od;
}

void EpaKernel::add_to_results(Rcpp::List& results){
  results["kernel.name"] = "epa";
  results["kernel.scale"] = this->scale;
}
