#include "RcppArmadillo.h"
#include "uint.h"

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef Kernel_hpp
#define Kernel_hpp

// Abstract Kernel class (interface)

class Kernel{
public:
  double scale, eval0;
  bool bounded = true;
  bool rescale_boundary = true;
  arma::rowvec t0;
  arma::rowvec factor;
  void update_scale(double scale);
  arma::mat eval(const arma::colvec &t1);
  // std::vector<arma::mat> eval(const std::vector<arma::rowvec> &t0, const arma::colvec &t1);
  std::vector<arma::mat> eval(const std::vector<arma::colvec> &t1);

  static Kernel* Create(std::string name, bool rescale_boundary, const arma::rowvec &t0);
  virtual ~Kernel(){};
  virtual void add_to_results(Rcpp::List& results) = 0;
  virtual arma::mat function(const arma::mat &od) = 0;
  virtual arma::rowvec rescale_factor() = 0;
};

class GaussianKernel : public Kernel{
public:
  GaussianKernel(bool rescale_boundary, const arma::rowvec &t0);
  arma::mat function(const arma::mat &od);
  void add_to_results(Rcpp::List& results);
  arma::rowvec rescale_factor();
};

class EpaKernel : public Kernel{
public:
  EpaKernel(bool rescale_boundary, const arma::rowvec &t0);
  arma::mat function(const arma::mat &od);
  void add_to_results(Rcpp::List& results);
  arma::rowvec rescale_factor();
};

class LinearKernel : public Kernel{
public:
  LinearKernel(bool rescale_boundary, const arma::rowvec &t0);
  arma::mat function(const arma::mat &od);
  void add_to_results(Rcpp::List& results);
  arma::rowvec rescale_factor();
};

class TriweightKernel : public Kernel{
public:
  TriweightKernel(bool rescale_boundary, const arma::rowvec &t0);
  arma::mat function(const arma::mat &od);
  void add_to_results(Rcpp::List& results);
  arma::rowvec rescale_factor();
};

#endif
