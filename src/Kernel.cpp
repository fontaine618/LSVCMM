#include "RcppArmadillo.h"
#include "Kernel.hpp"

//[[Rcpp::depends(RcppArmadillo)]]


Kernel::Kernel(){
  this->scale = 1.;
  this->eval0 = 0.5;
}

arma::mat Kernel::eval(const arma::rowvec &t0, const arma::colvec &t1){
  arma::mat od(t1.n_rows, t0.n_cols);
  arma::mat out(t1.n_rows, t0.n_cols, arma::fill::zeros);

  od.each_col() = t1;
  od.each_row() -= t0;

  out.elem(arma::find(arma::abs(od) <= this->scale)).ones();

  return out / (2. * this->scale);
}

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

void Kernel::add_to_results(Rcpp::List& results){
  results["kernel.name"] = "uniform";
  results["kernel.scale"] = this->scale;
}

