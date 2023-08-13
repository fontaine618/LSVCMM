#include "RcppArmadillo.h"
#include "uint.h"

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef Kernel_hpp
#define Kernel_hpp

// Abstract Kernel class (interface)

class Kernel{
public:
  double scale, eval0;
  virtual ~Kernel(){};
  virtual arma::mat eval(const arma::rowvec &t0, const arma::colvec &t1) = 0;
  std::vector<arma::mat> eval(const std::vector<arma::rowvec> &t0, const arma::colvec &t1);
  std::vector<arma::mat> eval(const arma::rowvec &t0, const std::vector<arma::colvec> &t1);
  virtual void add_to_results(Rcpp::List& results) = 0;
  static Kernel* Create(std::string name);
};

class GaussianKernel : public Kernel{
public:
  GaussianKernel();
  using Kernel::eval;
  arma::mat eval(const arma::rowvec &t0, const arma::colvec &t1);
  void add_to_results(Rcpp::List& results);
};

class EpaKernel : public Kernel{
public:
  EpaKernel();
  using Kernel::eval;
  arma::mat eval(const arma::rowvec &t0, const arma::colvec &t1);
  void add_to_results(Rcpp::List& results);
};

#endif
