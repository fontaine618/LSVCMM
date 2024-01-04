#include "RcppArmadillo.h"
#include "Data.h"
#include "uint.h"

// [[Rcpp::depends(RcppArmadillo)]]

std::vector<arma::mat> to_list_by_subject(
    const arma::uvec & subject,
    const arma::mat & array
){
  arma::uvec ids = unique(subject);
  std::vector<arma::mat> out(ids.n_elem);
  arma::ucolvec is;
  arma::mat rows;
  for(uint i : ids){
    is = arma::find(subject==i);
    rows = array.rows(is);
    out[i] = rows;
  }
  return out;
}

std::vector<arma::colvec> to_list_by_subject(
    const arma::uvec & subject,
    const arma::colvec & array
){
  arma::uvec ids = unique(subject);
  std::vector<arma::colvec> out(ids.n_elem);
  arma::ucolvec is;
  arma::colvec rows;
  for(uint i : ids){
    is = arma::find(subject==i);
    rows = array.rows(is);
    out[i] = rows;
  }
  return out;
}

Data::Data(
  const arma::colvec & response, // n x 1 in R
  const arma::ucolvec & subject, // n x 1 in {0, 1, ..., N-1}
  const arma::colvec & response_time, // n x 1 in R, but normally in [0,1] is scale_time=T from outside
  const arma::mat & vcm_covariates, // n x px in R, but typically 0/1 (vc intercept has to be added outside)
  const arma::mat & fixed_covariates, // n x pu in R, need to copy outside if constant, but allows changing covariates
  const arma::colvec & offset, // n x 1 in R, but normally in [0,1] is scale_time=T from outside
  const arma::colvec & weight // n x 1 in R, but normally in [0,1] is scale_time=T from outside
){
  // to list format and initialize object
  this->y = to_list_by_subject(subject, response);
  this->t = to_list_by_subject(subject, response_time);
  this->o = to_list_by_subject(subject, offset);
  this->w = to_list_by_subject(subject, weight);
  this->U = to_list_by_subject(subject, fixed_covariates);
  this->X = to_list_by_subject(subject, vcm_covariates);

  // store dimensions for quick access
  this->px = vcm_covariates.n_cols;
  this->pu = fixed_covariates.n_cols;
  this->n = response.n_elem;
  this->N = y.size();

  // initialize
  this->P = std::vector<arma::mat>(this->N);
  this->W = std::vector<arma::mat>(this->N);
  this->I = std::vector<arma::mat>(this->N);
  this->r = std::vector<arma::colvec>(this->N);
  this->m = std::vector<arma::colvec>(this->N);
  this->s = std::vector<arma::colvec>(this->N);
  this->lp = std::vector<arma::colvec>(this->N);
  this->sr = std::vector<arma::colvec>(this->N);
  this->sPsr = std::vector<arma::colvec>(this->N);
  this->foldid = arma::uvec(N, arma::fill::zeros);
}

Data::Data(
  const std::vector<arma::colvec> & y,
  const std::vector<arma::colvec> & t,
  const std::vector<arma::colvec> & o,
  const std::vector<arma::colvec> & w,
  const std::vector<arma::mat> & X,
  const std::vector<arma::mat> & U,
  const std::vector<arma::mat> & P,
  const std::vector<arma::mat> & I,
  const std::vector<arma::mat> & W
){
  // this constructor is to ease the get method for CV/Bootstrap
  this->y = y;
  this->t = t;
  this->o = o;
  this->w = w;
  this->U = U;
  this->X = X;
  this->P = P;
  this->W = W;
  this->I = I;
  this->px = X[0].n_cols;
  this->pu = U[0].n_cols;
  this->N = y.size();
  uint n = 0;
  for(uint i=0; i<y.size(); i++) n += y[i].n_elem;
  this->n = n;

  // initialize
  this->r = std::vector<arma::colvec>(this->N);
  this->m = std::vector<arma::colvec>(this->N);
  this->s = std::vector<arma::colvec>(this->N);
  this->lp = std::vector<arma::colvec>(this->N);
  this->sr = std::vector<arma::colvec>(this->N);
  this->sPsr = std::vector<arma::colvec>(this->N);
  this->foldid = arma::uvec(N, arma::fill::zeros);
}

void Data::prepare_folds(uint nfolds){
  arma::uvec order = arma::randperm(this->N);
  uint N_per_fold = this->N / nfolds;
  arma::uvec foldid(this->N);
  uint N_in_current_fold = 0;
  uint current_fold = 0;
  for(uint i=0; i<this->N; i++){
    // NB last fold may be bigger/smaller than firsts
    if(N_in_current_fold>=N_per_fold && current_fold<nfolds-1){
      current_fold ++;
      N_in_current_fold = 0;
    }
    foldid[order[i]] = current_fold;
    N_in_current_fold++;
  }
  this->foldid = foldid;
}

Data Data::get(arma::uvec ids){
  N = ids.n_elem;
  std::vector<arma::colvec> y(N);
  std::vector<arma::colvec> t(N);
  std::vector<arma::colvec> o(N);
  std::vector<arma::colvec> w(N);
  std::vector<arma::mat> X(N);
  std::vector<arma::mat> U(N);
  std::vector<arma::mat> P(N);
  std::vector<arma::mat> I(N);
  std::vector<arma::mat> W(N);
  for(uint i=0; i<N; i++){
    y[i] = this->y[ids[i]];
    t[i] = this->t[ids[i]];
    o[i] = this->o[ids[i]];
    w[i] = this->w[ids[i]];
    X[i] = this->X[ids[i]];
    U[i] = this->U[ids[i]];
    P[i] = this->P[ids[i]];
    I[i] = this->I[ids[i]];
    W[i] = this->W[ids[i]];
  }
  return Data(y, t, o, w, X, U, P, I, W);
}

Data Data::get_fold(uint fold){
  arma::uvec ids = arma::find(this->foldid == fold);
  return this->get(ids);
}

Data Data::get_other_folds(uint fold){
  arma::uvec ids = arma::find(this->foldid != fold);
  return this->get(ids);
}

Data Data::resample(bool resample_within_subject){
  uint N = this->y.size();
  arma::ivec ids = arma::randi(N, arma::distr_param(0, N-1));
  arma::uvec ids2 = arma::conv_to<arma::uvec>::from(ids);
  Data out = this->get(ids2);
  if(resample_within_subject) return out.resample_within_subject();
  return out;
}

Data Data::resample_within_subject(){
  std::vector<arma::colvec> y(N);
  std::vector<arma::colvec> t(N);
  std::vector<arma::colvec> o(N);
  std::vector<arma::colvec> w(N);
  std::vector<arma::mat> X(N);
  std::vector<arma::mat> U(N);
  std::vector<arma::mat> P(N);
  std::vector<arma::mat> I(N);
  std::vector<arma::mat> W(N);
  for(uint i=0; i<N; i++){
    uint n = this->y[i].n_elem;
    arma::ivec ids_ = arma::randi(n, arma::distr_param(0, n-1));
    arma::uvec ids = arma::conv_to<arma::uvec>::from(ids_);
    y[i] = this->y[i].elem(ids);
    t[i] = this->t[i].elem(ids);
    o[i] = this->o[i].elem(ids);
    w[i] = this->w[i].elem(ids);
    X[i] = this->X[i].rows(ids);
    U[i] = this->U[i].rows(ids);
    P[i] = this->P[i].rows(ids);
    I[i] = this->I[i].rows(ids);
    W[i] = this->W[i].rows(ids);
  }
  return Data(y, t, o, w, X, U, P, I, W);
}
