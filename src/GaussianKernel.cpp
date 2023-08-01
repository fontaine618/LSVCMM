#include "RcppArmadillo.h"
#include "GaussianKernel.h"
#include "uint.h"

// [[Rcpp::depends(RcppArmadillo)]]

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

std::vector<arma::mat> GaussianKernel::eval(const std::vector<arma::rowvec> &t0, const arma::colvec &t1){
  std::vector<arma::mat> out(t0.size());
  for(uint i=0; i<t0.size(); i++) out[i] = this->eval(t0[i], t1);
  return out;
}

std::vector<arma::mat> GaussianKernel::eval(const arma::rowvec &t0, const std::vector<arma::colvec> &t1){
  std::vector<arma::mat> out(t1.size());
  for(uint i=0; i<t1.size(); i++) out[i] = this->eval(t0, t1[i]);
  return out;
}

void GaussianKernel::add_to_results(Rcpp::List& results){
  results["kernel.name"] = "gaussian";
  results["kernel.scale"] = this->scale;
}
