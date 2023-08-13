#include "RcppArmadillo.h"
#include "WorkingCovariance.h"
#include "Logger.cpp"
#include "Control.cpp"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]

// =============================================================================
// Interface
WorkingCovariance* WorkingCovariance::Create(
    std::string name,
    double variance_ratio,
    double correlation,
    bool estimate_parameters
  ){
  if(name=="compound_symmetry") return new CompoundSymmetry(variance_ratio, 1., estimate_parameters);
  Rcpp::stop("unrecognized working covariance name: " + name);
}

std::vector<arma::mat> WorkingCovariance::compute_precision(const std::vector<arma::colvec> &time){
  std::vector<arma::mat> out(time.size());
  for(uint i=0; i<time.size(); i++) out[i] = this->compute_precision(time[i]);
  return out;
}

// TODO: this is a bit wrong since it does not use the dispersion
double WorkingCovariance::profile_likelihood(
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

// =============================================================================
// Compound Symmetry
CompoundSymmetry::CompoundSymmetry(){
  this->variance_ratio = 1.0;
  this->correlation = 1.0;
  this->estimate_parameters = TRUE;
}
CompoundSymmetry::CompoundSymmetry(double variance_ratio, double correlation, bool estimate_parameters){
  this->variance_ratio = variance_ratio;
  this->correlation = correlation;
  this->estimate_parameters = estimate_parameters;
}

arma::mat CompoundSymmetry::compute_precision(const arma::colvec &time){
  uint ni = time.n_elem;
  arma::mat out = arma::eye(ni, ni);
  out -= 1.0/(ni + 1.0/this->variance_ratio);
  return out;
}

std::vector<double> CompoundSymmetry::derivatives(
    const std::vector<arma::colvec> &sr,
    const std::vector<arma::colvec> &t,
    const std::vector<arma::mat> &P,
    const double dispersion
){
  double d1 = 0.;
  double d2 = 0.;
  for(uint i=0; i<sr.size(); i++){
    uint ni = sr[i].n_elem;
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

uint CompoundSymmetry::update_parameters(
    const std::vector<arma::colvec> &sr,
    const std::vector<arma::colvec> &t,
    const std::vector<arma::mat> &P,
    const double dispersion,
    Logger* logger,
    uint round,
    Control *control
){
  if(!this->estimate_parameters) return 0;
  double pllk = this->profile_likelihood(sr, P);
  double pllk_prev = pllk;
  std::vector<arma::mat> Ptmp = P;
  uint iter;
  for(iter=0; iter<control->max_iter; iter++){
    std::vector<double> d = this->derivatives(sr, t, Ptmp, dispersion);
    if(d[1] > 0.){
      this->variance_ratio *= (d[0] > 0.) ? 2. : 0.5; // jump somewhere else
    }else{
      double step = -d[0] / d[1];
      this->variance_ratio += step;
    }
    this->variance_ratio = fmax(this->variance_ratio, 1e-6);
    Ptmp = this->compute_precision(t);
    pllk = this->profile_likelihood(sr, Ptmp);
    logger->add_variance_iteration_results(round, iter, pllk);
    if(control->verbose > 2) Rcpp::Rcout << "          " << round <<".V." << iter <<
      " obj=" << pllk << "\n";
    if(fabs(pllk - pllk_prev) / fabs(pllk) < control->rel_tol) {
      break;
    }
    pllk_prev = pllk;
  }
  return iter;
}

void CompoundSymmetry::add_to_results(Rcpp::List& results){
  results["working_covariance.estimate"] = this->estimate_parameters;
  results["working_covariance.ratio"] = this->variance_ratio;
  results["working_covariance.name"] = "compound_symmetry";
}
