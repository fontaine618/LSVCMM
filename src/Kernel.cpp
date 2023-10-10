#include "RcppArmadillo.h"
#include "Kernel.h"
#include "uint.h"

// [[Rcpp::depends(RcppArmadillo)]]

// =============================================================================
// Kernel
// std::vector<arma::mat> Kernel::eval(const std::vector<arma::rowvec> &t0, const arma::colvec &t1){
//   std::vector<arma::mat> out(t0.size());
//   for(uint i=0; i<t0.size(); i++) out[i] = this->eval(t0[i], t1);
//   return out;
// }

std::vector<arma::mat> Kernel::eval(const std::vector<arma::colvec> &t1){
  std::vector<arma::mat> out(t1.size());
  for(uint i=0; i<t1.size(); i++) out[i] = this->eval(t1[i]);
  return out;
}

Kernel* Kernel::Create(std::string name, bool rescale_boundary, const arma::rowvec &t0){
  if(name=="gaussian") return new GaussianKernel(rescale_boundary, t0);
  if(name=="epa") return new EpaKernel(rescale_boundary, t0);
  if(name=="linear") return new LinearKernel(rescale_boundary, t0);
  if(name=="triweight") return new TriweightKernel(rescale_boundary, t0);
  Rcpp::stop("unrecognized kernel name: " + name);
}

arma::mat Kernel::eval(const arma::colvec &t1){
  arma::mat od(t1.n_rows, t0.n_cols);
  od.each_col() = t1;
  od.each_row() -= t0;
  od = arma::abs(od);
  arma::mat out = this->function(od / this->scale) / this->scale;
  if(this->bounded) out.elem(arma::find(od>this->scale)).zeros();
  if(this->rescale_boundary){
    // recall that first argument is the estimated time, second is the observed time
    arma::rowvec factor = this->rescale_factor();
    // factor.print("rescaling factor for kernel boundary fix");
    out.each_row() %= factor;
  }
  return out;
}

void Kernel::update_scale(double scale){
  this->scale = scale;
  this->factor = this->rescale_factor();
}

// =============================================================================
// Gaussian Kernel
GaussianKernel::GaussianKernel(bool rescale_boundary, const arma::rowvec &t0){
  this->rescale_boundary = rescale_boundary;
  this->t0 = t0;
  this->bounded = false;
  this->scale = 1.;
  this->eval0 = 1. / sqrt(arma::datum::pi);
}

arma::mat GaussianKernel::function(const arma::mat &od){
  return arma::exp(- arma::square(od) ) / sqrt(arma::datum::pi);
}

void GaussianKernel::add_to_results(Rcpp::List& results){
  results["kernel.name"] = "gaussian";
  results["kernel.scale"] = this->scale;
}

arma::rowvec GaussianKernel::rescale_factor(){
  arma::rowvec U = sqrt(arma::datum::pi) * (1. - t0) / this->scale;
  arma::rowvec L = sqrt(arma::datum::pi) * (0. - t0) / this->scale;
  return 1. / (arma::normcdf(U) - arma::normcdf(L));
}

// =============================================================================
// Epanechnikov Kernel
EpaKernel::EpaKernel(bool rescale_boundary, const arma::rowvec &t0){
  this->rescale_boundary = rescale_boundary;
  this->t0 = t0;
  this->scale = 1.;
  this->eval0 = 0.75;
}

arma::mat EpaKernel::function(const arma::mat &od){
  return 0.75 * (1. - arma::square(od) );
}

void EpaKernel::add_to_results(Rcpp::List& results){
  results["kernel.name"] = "epa";
  results["kernel.scale"] = this->scale;
}

arma::rowvec EpaKernel::rescale_factor(){
  return t0*0. + 1.;
}

// =============================================================================
// Linear Kernel (Triangular)
LinearKernel::LinearKernel(bool rescale_boundary, const arma::rowvec &t0){
  this->rescale_boundary = rescale_boundary;
  this->t0 = t0;
  this->scale = 1.;
  this->eval0 = 1.;
}

arma::mat LinearKernel::function(const arma::mat &od){
  return (1. - od);
}

void LinearKernel::add_to_results(Rcpp::List& results){
  results["kernel.name"] = "linear";
  results["kernel.scale"] = this->scale;
}

arma::rowvec LinearKernel::rescale_factor(){
  return t0*0. + 1.;
}


// =============================================================================
// Triweight Kernel
TriweightKernel::TriweightKernel(bool rescale_boundary, const arma::rowvec &t0){
  this->rescale_boundary = rescale_boundary;
  this->t0 = t0;
  this->scale = 1.;
  this->eval0 = 35./32.;
}

arma::mat TriweightKernel::function(const arma::mat &od){
  return (35./32.) * arma::pow(1. - arma::square(od),  3);
}

void TriweightKernel::add_to_results(Rcpp::List& results){
  results["kernel.name"] = "triweight";
  results["kernel.scale"] = this->scale;
}

arma::rowvec TriweightKernel::rescale_factor(){
  return t0*0. + 1.;
}
