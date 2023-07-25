#include "RcppArmadillo.h"
#include "WorkingCovariance.hpp"

//[[Rcpp::depends(RcppArmadillo)]]


std::vector<arma::mat> CompoundSymmetry::compute_covariance(const std::vector<arma::colvec> &time){
  std::vector<arma::mat> out(time.size());
  for(uint i=0; i<time.size(); i++) out[i] = this->compute_covariance(time[i]);
  return out;
}
std::vector<arma::mat> CompoundSymmetry::compute_precision(const std::vector<arma::colvec> &time){
  std::vector<arma::mat> out(time.size());
  for(uint i=0; i<time.size(); i++) out[i] = this->compute_precision(time[i]);
  return out;
}

CompoundSymmetry::CompoundSymmetry(){}

CompoundSymmetry::CompoundSymmetry(double variance_ratio, bool estimate_parameters){
  this->variance_ratio = variance_ratio;
  this->estimate_parameters = estimate_parameters;
}

arma::mat CompoundSymmetry::compute_covariance(const arma::colvec &time){
  arma::mat out = arma::eye(time.n_elem, time.n_elem);
  out += this->variance_ratio;
  return out;
}

arma::mat CompoundSymmetry::compute_precision(const arma::colvec &time){
  uint ni = time.n_elem;
  arma::mat out = arma::eye(ni, ni);
  out -= 1.0/(ni + 1.0/this->variance_ratio);
  return out;
}

double CompoundSymmetry::profile_likelihood(
    const std::vector<arma::colvec> &sr,
    const std::vector<arma::mat> &P
){
  double pllk = 0.;
  for(uint i=0; i<sr.size(); i++){
    pllk += arma::dot(sr[i], P[i] * sr[i]);
    pllk -= arma::log_det_sympd(P[i]);
  }
  return -0.5 * pllk;
}

std::vector<double> CompoundSymmetry::derivatives(
    const std::vector<arma::colvec> &sr,
    const std::vector<arma::mat> &P,
    const double dispersion
){
  double d1 = 0.;
  double d2 = 0.;
  for(uint i=0; i<sr.size(); i++){
    uint ni = sr[i].n_elem;
    // arma::mat dV1 = arma::mat(ni, ni, arma::fill::ones);
    // arma::mat dV2 = arma::mat(ni, ni, arma::fill::zeros);
    // arma::mat dP1 = -P[i] * dV1 * P[i];
    // arma::mat PdV1 = P[i] * dV1;
    // arma::mat PdV2 = P[i] * dV2;
    // arma::mat dP1dV1 = dP1 * dV1;
    // arma::mat dP2 = - dP1dV1 * P[i] - PdV2 * P[i] - PdV1 * dP1;
    // d1 += arma::trace(PdV1);
    // d1 += - arma::dot(sr[i], dP1 * sr[i]);
    // d2 += arma::trace(dP1dV1);
    // d2 += arma::trace(PdV2);
    // d2 += arma::dot(sr[i], dP2 * sr[i]);
    arma::mat dV = arma::mat(ni, ni, arma::fill::ones);
    arma::mat PdV = P[i] * dV;
    arma::mat PdVPsr = PdV * P[i] * sr[i];
    double d1tr = arma::trace(PdV);
    double d1ip = arma::dot(sr[i], PdVPsr);
    double d2tr = arma::trace(PdV * PdV);
    double d2ip = arma::dot(sr[i], PdV * PdVPsr);
    d1 += d1tr - d1ip / dispersion;
    d2 += -d2tr + 2. * d2ip / dispersion;
  }
  return std::vector<double>{-0.5 * d1, -0.5 * d2};
}

void CompoundSymmetry::update_parameters(
    const std::vector<arma::colvec> &sr,
    const std::vector<arma::colvec> &t,
    const std::vector<arma::mat> &P,
    const double dispersion
){
  if(!this->estimate_parameters) return;
  double pllk = this->profile_likelihood(sr, P);
  double pllk_prev = pllk;
  std::vector<arma::mat> Ptmp = P;
  for(uint iter=0; iter<20; iter++){
    std::vector<double> d = this->derivatives(sr, Ptmp, dispersion);
    if(d[1] > 0.){
      this->variance_ratio *= (d[0] > 0.) ? 2. : 0.5; // jump somewhere else
    }else{
      double step = -d[0] / d[1];
      this->variance_ratio += step;
    }
    this->variance_ratio = fmax(this->variance_ratio, 1e-6);
    Ptmp = this->compute_precision(t);
    pllk = this->profile_likelihood(sr, Ptmp);
    if(fabs(pllk - pllk_prev) / fabs(pllk) < 1e-6) {
      Rcpp::Rcout << "          .V." << iter << " pllk=" << pllk << " (re_ratio=" << this->variance_ratio <<
        ", d1=" << d[0] << ", d2=" << d[1] << ")\n";
      return;
    }
    pllk_prev = pllk;
  }
  // NB: precision update is not stored, need to redo it outside
}

void CompoundSymmetry::add_to_results(Rcpp::List& results){
  results["working_covariance.variance_ratio"] = this->variance_ratio;
  results["working_covariance.name"] = "compound_symmetry";
}
