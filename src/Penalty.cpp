#include "RcppArmadillo.h"
#include "Penalty.hpp"

//[[Rcpp::depends(RcppArmadillo)]]

arma::rowvec proximal_L1(
    const arma::rowvec &b,
    const arma::rowvec &m
){
  uint nt = b.n_elem;
  arma::rowvec s(nt);
  s = arma::abs(b) - m;
  s = arma::sign(b) % arma::clamp(s, 0., arma::datum::inf);
  return s;
}

arma::rowvec proximal_L2(
    const arma::rowvec &s,
    const double m
){
  double sn = arma::norm(s);
  if(sn == 0.) sn = 1.; // to avoid dividing by 0, the value doesn't matter since s=0 in that case
  return fmax(1. - m/sn, 0.) * s;
}

arma::rowvec proximal_L1L2_row(
    const arma::rowvec &b,
    const arma::rowvec &m1,
    const double m2
){
  return proximal_L2(proximal_L1(b, m1), m2);
}

arma::mat proximal_L1L2(
    const arma::mat &B,
    const arma::mat &M1,
    const arma::colvec &m2
){
  arma::mat out = B * 0.;
  for(uint j=0; j<B.n_rows; j++){
    out.row(j) = proximal_L1L2_row(B.row(j), M1.row(j), m2(j));
  }
  return out;
}

Penalty::Penalty(
  double lambda,
  double alpha,
  double power,
  bool penalize_intercept
){
  this->lambda = lambda;
  this->alpha = alpha;
  this->power = power;
  this->penalize_intercept = penalize_intercept;
}

Penalty::Penalty(){}

arma::mat Penalty::proximal(const arma::mat &B, double stepsize){
  arma::mat M1 = this->W1 * this->alpha * this->lambda * stepsize;
  arma::colvec m2 = this->w2 * (1 - this->alpha) * this->lambda * stepsize * sqrt(B.n_cols);
  return proximal_L1L2(B, M1, m2);
}

void Penalty::update_weights(const arma::mat &B){
  this->W1 = arma::pow(arma::abs(B), -this->power);
  this->w2 = arma::pow(arma::sum(arma::pow(B, 2), 1), -this->power/2);
  if(!this->penalize_intercept){
    this->W1.row(0) = 0.;
    this->w2(0) = 0.;
  }
}

void Penalty::unit_weights(const arma::mat &B){
  this->W1 = arma::mat(B.n_rows, B.n_cols, arma::fill::zeros);
  this->w2 = arma::colvec(B.n_rows, arma::fill::zeros);
}

void Penalty::add_to_results(Rcpp::List& results){
  results["penalty.alpha"] = this->alpha;
  results["penalty.lambda"] = this->lambda;
  results["penalty.adaptive"] = this->power;
  results["penalty.name"] = "adaptive_sparse_group_lasso";
}

double Penalty::lambda_max(const arma::mat &gB){
  double lambda_max = 0.;
  double num, denum, gnorm, wnorm;
  for(uint j=0; j<gB.n_rows; j++){
    // L1 bound
    if(this->alpha > 0.){
      for(uint t=0; t<gB.n_cols; t++){
        if(this->W1(j, t) > 0.){
          num = fabs(gB(j, t));
          denum = this->alpha * this->W1(j, t) + (1-this->alpha) * this->w2(j);
          lambda_max = fmax(lambda_max, num / denum);
        }
      }
    }
    // L2 bound
    if(this->alpha < 1.){
      gnorm = arma::norm(gB.row(j), "fro");
      wnorm = arma::norm(this->W1.row(j), "fro");
      lambda_max = fmax(lambda_max, gnorm / (this->alpha * wnorm + (1-this->alpha) * this->w2(j)));
    }
  }
  return lambda_max;
}
