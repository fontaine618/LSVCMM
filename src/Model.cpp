#include "RcppArmadillo.h"
#include "Data.hpp"
#include "Kernel.hpp"
#include "LinkFunction.cpp"
#include "Family.hpp"
#include "WorkingCovariance.hpp"
#include "Penalty.hpp"
#include "Model.hpp"

//[[Rcpp::depends(RcppArmadillo)]]

Model::Model(
  arma::rowvec estimated_time,
  Penalty penalty,
  WorkingCovariance workingCovariance,
  LinkFunction linkFunction,
  Family family,
  Control control,
  Kernel kernel
){
  this->estimated_time = estimated_time;
  this->nt = estimated_time.n_elem;
  this->penalty = penalty;
  this->workingCovariance = workingCovariance;
  this->linkFunction = linkFunction;
  this->family = family;
  this->control = control;
  this->kernel = kernel;
}

Model::Model(){}

std::vector<arma::colvec> Model::_linear_predictor(Data data){
  std::vector<arma::colvec> lp(data.N);
  for(uint i=0; i<data.N; i++){
    lp[i] = data.o[i];
    lp[i] += arma::sum((data.X[i] * this->B) % data.I[i], 1);
    lp[i] += data.U[i] * this->a;
  }
  return lp;
}

void Model::update_mean(Data data){
  data.lp = this->_linear_predictor(data);
  data.m = this->linkFunction.eval(data.lp);
  for(uint i=0; i<data.N; i++){
    data.r[i] = data.y[i] - data.m[i];
    data.s[i] = arma::sqrt(this->family.unit_variance(data.m[i]));
    data.sr[i] = data.r[i] / data.s[i];
    data.sPsr[i] = (data.P[i] * data.sr[i]) / data.s[i];
  }
}

void Model::update_precision(Data data){
  data.P = this->workingCovariance.compute_precision(data.t);
}

double Model::quadratic_term(const Data data){
  double q = 0;
  for(uint i=0; i<data.N; i++){
    q += arma::dot(data.r[i], data.sPsr[i]);
  }
  return q;
}

double Model::logdet_precision(const Data data){
  double ld = 0;
  for(uint i=0; i<data.N; i++){
    ld += arma::log_det_sympd(data.P[i]);
  }
  return ld;
}

double Model::logdet_scaling(const Data data){
  double ld = 0;
  for(uint i=0; i<data.N; i++){
    ld += arma::accu(arma::log(data.s[i]));
  }
  return ld;
}

double Model::logdet_dispersion(const Data data){
  double ld = 0;
  for(uint i=0; i<data.N; i++){
    ld += data.y[i].n_elem * log(2*arma::datum::pi*this->family.dispersion);
  }
  return ld;
}

double Model::quasi_log_likelihood(double quad_term, double logdet_term){
  return -0.5 * (quad_term + logdet_term);
}

void Model::update_gradients(const Data data){
  this->gB.zeros();
  this->ga.zeros();
  uint nt = this->B.n_cols;
  std::vector<arma::colvec> d = this->linkFunction.derivative(data.lp);
  for(uint i=0; i<data.N; i++){
    arma::colvec tmp = d[i] % data.sPsr[i];
    for(uint j=0; j<nt; j++){
      this->gB += data.X[i].t() * (data.W[i].col(j) % tmp);
    }
    this->ga += data.U[i].t() * data.sPsr[i];
  }
}

void Model::proximal_gradient_step(){
  this->a += this->ga / this->La;
  this->B += this->gB / this->LB;
  this->B = this->penalty.proximal(this->B, this->LB);
}

void Model::initialize(Data data){}

arma::mat Model::hessian(const Data data){
  // assumes mean and precision are updated

  // FIXME: not sure if this is correct for non-identity link?
  // should there be a second term?

  // compute all blocks
  arma::mat Haa = arma::mat(data.pu, data.pu, arma::fill::zeros);
  std::vector<arma::mat> HaB(this->nt);
  std::vector<arma::mat> HBB(this->nt);
  for(uint t=0; t<this->nt; t++){
    HaB[t] = arma::mat(data.pu, data.px, arma::fill::zeros);
    HBB[t] = arma::mat(data.px, data.px, arma::fill::zeros);
  }
  std::vector<arma::colvec> d = this->linkFunction.derivative(data.lp);
  for(uint i=0; i<data.N; i++){
    arma::colvec ds = d[i] / data.s[i];
    for(uint t=0; t<this->nt; t++){
      arma::colvec wds = ds % data.W[i].col(t);
      arma::mat wdsPsdw = data.P[i];
      wdsPsdw.each_row() %= wds;
      wdsPsdw.each_col() %= wds;
      arma::mat wdsPsdwX = wdsPsdw * data.X[i];
      arma::mat wdsPsdwU = wdsPsdw * data.U[i];
      Haa += data.U[i].t() * wdsPsdwU;
      HaB[t] += data.U[i].t() * wdsPsdwX;
      HBB[t] += data.X[i].t() * wdsPsdwX;
    }
  }

  // construct full Hessian
  uint dim = data.pu + data.px*this->nt;
  arma::mat hessian = arma::mat(dim, dim, arma::fill::zeros);

  if(data.pu > 0) hessian.submat(0, 0, data.pu - 1, data.pu - 1) = Haa;
  for(uint k=0; k<this->nt; k++){
    if(data.pu > 0){
      hessian.submat(0,
                     data.pu + k*data.px,
                     data.pu - 1,
                     data.pu + (k+1)*data.px - 1) = HaB[k];
      hessian.submat(data.pu + k*data.px,
                     0,
                     data.pu + (k+1)*data.px - 1,
                     data.pu - 1) = HaB[k].t();
    }
    hessian.submat(data.pu + k*data.px,
                   data.pu + k*data.px,
                   data.pu + (k+1)*data.px - 1,
                   data.pu + (k+1)*data.px - 1) = HBB[k];
  }
  return hessian;
}

void Model::prepare_stepsize(Data data){
  arma::mat hessian = this->hessian(data);
  uint dim = data.pu + data.px*this->nt;
  arma::mat Haa = hessian.submat(0, 0, data.pu - 1, data.pu - 1);
  arma::mat HBB = hessian.submat(data.pu - 1, data.pu - 1, dim - 1, dim - 1);

  // Find heuristic values for La, Lb
  arma::vec eigval;

  double La = 0.;
  if(data.pu > 0) La = arma::eig_sym(Haa).max();

  double Lb = arma::eig_sym(HBB).max();

  // Heuristic might be a little small, increase until fine
  arma::mat L = arma::eye(dim, dim) * LB;
  for(uint j=0; j<data.pu; j++) L(j, j) = La;

  arma::mat LmH = L - hessian;
  double factor = 1.;
  while(!LmH.is_sympd()){
    factor *= 1.01;
    L *= 1.01;
    LmH = L - hessian;
    Rcpp::Rcout << "         factor=" << factor << "\n";
  }

  this->La = La*factor;
  this->LB = Lb*factor;
}

double Model::lambda_max(Data data){
  this->penalty.lambda = 1e10;
  this->fit(data);
  return this->penalty.lambda_max(this->gB);
}

void Model::fit(Data data){
  this->update_precision(data);
  this->update_mean(data);
  double quad_term = this->quadratic_term(data);
  double ld_scaling = this->logdet_scaling(data);
  double ld_precision = this->logdet_precision(data);
  double ld_dispersion = this->logdet_dispersion(data);
  double llk = this->quasi_log_likelihood(quad_term, ld_scaling + ld_dispersion - ld_precision);
  double llk_old = llk;
  for(uint round=0; round<this->control.max_rounds; round++){
    // mean update inner loop
    double quad_term_old = quad_term;
    for(uint iter=0; iter<this->control.max_iter; iter++){
      this->update_gradients(data);
      this->proximal_gradient_step();
      this->update_mean(data);
      quad_term = this->quadratic_term(data);
      Rcpp::Rcout << "         " << round << "." << "M" << "." << iter << ": obj=" << quad_term << "\n";
      if(fabs(quad_term - quad_term_old) / fabs(quad_term_old)< this->control.rel_tol) break;
      quad_term_old = quad_term;
    }
    // variance update
    this->workingCovariance.update_parameters(data.sr);
    this->update_precision(data);
    quad_term = this->quadratic_term(data);
    this->family.dispersion = quad_term / data.n;
    ld_scaling = this->logdet_scaling(data);
    ld_precision = this->logdet_precision(data);
    ld_dispersion = this->logdet_dispersion(data);
    llk = this->quasi_log_likelihood(quad_term, ld_scaling + ld_dispersion - ld_precision);
    Rcpp::Rcout << "         " << round << " llk=" << llk << "\n";
    if(fabs(llk - llk_old) / fabs(llk_old) < this->control.rel_tol) break;
    llk_old = llk;
  }
  this->results["llk"] = llk;
  this->results["rss"] = quad_term;
  this->prepare_results(data);
}

void Model::prepare_results(const Data data){
  this->family.add_to_results(this->results);
  this->penalty.add_to_results(this->results);
  this->linkFunction.add_to_results(this->results);
  this->workingCovariance.add_to_results(this->results);
  this->kernel.add_to_results(this->results);
  this->prepare_ics(data);
}

void Model::prepare_ics(const Data data){
  arma::rowvec nhat = arma::mat(1, this->B.n_cols, arma::fill::zeros);
  for(uint i=0; i<data.N; i++){
    nhat += arma::sum(data.W[i], 0);
  }
  arma::mat active = this->B *0.;
  active.elem(arma::find(arma::abs(this->B) > 0)).ones();
  arma::rowvec df = arma::sum(active, 0);
  arma::rowvec df_kernel = df * kernel.eval0 / (kernel.scale * df.n_elem);

  results["df"] = arma::accu(df);
  results["df_kernel"] = arma::accu(df_kernel);
  results["df_logn"] = arma::dot(df, arma::log(nhat));
  results["df_logn_kernel"] = arma::dot(df_kernel, arma::log(nhat));

  results["aic"] = -2 * (double)results["llk"] + 2 * (double)results["df"];
  results["bic"] = -2 * (double)results["llk"] + (double)results["df_logn"];
  results["aich"] = -2 * (double)results["llk"] + 2 * (double)results["df_kernel"];
  results["bich"] = -2 * (double)results["llk"] + (double)results["df_logn_kernel"];
}

Rcpp::List Model::save(){
  return Rcpp::List::create(
    Rcpp::Named("a", this->a),
    Rcpp::Named("B", this->B),
    Rcpp::Named("results", this->results)
  );
}
