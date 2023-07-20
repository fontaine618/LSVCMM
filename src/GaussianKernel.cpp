#include "RcppArmadillo.h"
#include "GaussianKernel.hpp"

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

void GaussianKernel::add_to_results(Rcpp::List& results){
  results["kernel.name"] = "gaussian";
  results["kernel.scale"] = this->scale;
}
