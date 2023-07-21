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
  uint px,
  uint pu,
  arma::rowvec estimated_time,
  Penalty penalty,
  WorkingCovariance workingCovariance,
  LinkFunction linkFunction,
  Family family,
  Control control,
  Kernel kernel
){
  this->px = px;
  this->pu = pu;
  this->estimated_time = estimated_time;
  this->nt = estimated_time.n_elem;
  this->penalty = penalty;
  this->workingCovariance = workingCovariance;
  this->linkFunction = linkFunction;
  this->family = family;
  this->control = control;
  this->kernel = kernel;

  this->B = arma::mat(px, estimated_time.n_elem, arma::fill::zeros);
  this->a = arma::colvec(pu, arma::fill::zeros);
  this->gB = arma::mat(px, estimated_time.n_elem, arma::fill::zeros);
  this->ga = arma::colvec(pu, arma::fill::zeros);
}

Model::Model(){}

std::vector<arma::colvec> Model::linear_predictor(Data &data){
  std::vector<arma::colvec> lp(data.N);
  for(uint i=0; i<data.N; i++){
    lp[i] = data.o[i];
    lp[i] += arma::sum((data.X[i] * this->B) % data.I[i], 1);
    if(this->pu) lp[i] += data.U[i] * this->a;
  }
  return lp;
}

void Model::update_mean(Data &data){
  data.lp = this->linear_predictor(data);
  data.m = this->linkFunction.eval(data.lp);
  for(uint i=0; i<data.N; i++){
    data.r[i] = data.y[i] - data.m[i];
    data.s[i] = arma::sqrt(this->family.unit_variance(data.m[i]));
    data.sr[i] = data.r[i] / data.s[i];
    data.sPsr[i] = (data.P[i] * data.sr[i]) / data.s[i];
  }
}

void Model::update_precision(Data &data){
  data.P = this->workingCovariance.compute_precision(data.t);
}

double Model::quadratic_term(const Data &data){
  double q = 0;
  for(uint i=0; i<data.N; i++){
    q += arma::dot(data.r[i], data.sPsr[i]);
  }
  return q;
}

double Model::logdet_precision(const Data &data){
  double ld = 0;
  for(uint i=0; i<data.N; i++){
    ld += arma::log_det_sympd(data.P[i]);
  }
  return ld;
}

double Model::logdet_scaling(const Data &data){
  double ld = 0;
  for(uint i=0; i<data.N; i++){
    ld += arma::accu(arma::log(data.s[i]));
  }
  return ld;
}

double Model::logdet_dispersion(const Data &data){
  double ld = 0;
  for(uint i=0; i<data.N; i++){
    ld += data.y[i].n_elem * log(2*arma::datum::pi*this->family.dispersion);
  }
  return ld;
}

double Model::quasi_log_likelihood(double quad_term, double logdet_term){
  return -0.5 * (quad_term + logdet_term);
}

void Model::update_gradients(const Data &data){
  this->gB.zeros();
  this->ga.zeros();
  uint nt = this->B.n_cols;
  std::vector<arma::colvec> d = this->linkFunction.derivative(data.lp);
  for(uint i=0; i<data.N; i++){
    arma::colvec dsPsr = d[i] % data.sPsr[i];
    for(uint j=0; j<nt; j++){
      arma::vec wdsPsr = data.W[i].col(j) % dsPsr;
      this->gB.col(j) += data.X[i].t() * wdsPsr;
    }
    if(data.pu) this->ga += data.U[i].t() * dsPsr;
  }
}

void Model::proximal_gradient_step(){
  this->a -= this->ga / this->La;
  this->B -= this->gB / this->LB;
  this->B = this->penalty.proximal(this->B, this->LB);
}

void Model::initialize(Data &data){}

arma::mat Model::hessian(const Data &data){
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
    arma::mat dsPsd = data.P[i];
    dsPsd.each_row() %= ds.t();
    dsPsd.each_col() %= ds;
    arma::mat dsPsdU = dsPsd * data.U[i];
    Haa = data.U[i].t() * dsPsdU;
    for(uint t=0; t<this->nt; t++){
      arma::colvec w = data.W[i].col(t);
      arma::mat wX = data.X[i];
      wX.each_col() %= w;
      arma::mat dsPsdwX = dsPsd * wX;
      HaB[t] += data.U[i].t() * dsPsdwX;
      HBB[t] += wX.t() * dsPsdwX;
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

void Model::prepare_stepsize(Data &data){
  arma::mat hessian = this->hessian(data);
  uint dim = data.pu + data.px*this->nt;
  arma::mat Haa;
  if(data.pu > 0) Haa = hessian.submat(0, 0, data.pu - 1, data.pu - 1);
  arma::mat HBB = hessian.submat(data.pu, data.pu, dim - 1, dim - 1);

  // Find heuristic values for La, Lb
  arma::vec eigval;

  double La = 0.;
  if(data.pu > 0) La = arma::eig_sym(Haa).max();

  double LB = arma::eig_sym(HBB).max();

  // Heuristic might be a little small, increase until fine
  arma::mat L = arma::eye(dim, dim) * LB;
  for(uint j=0; j<data.pu; j++) L(j, j) = La;

  arma::mat LmH = L - hessian;
  double factor = 1.;
  while(!LmH.is_sympd()){
    factor *= 1.01;
    L *= 1.01;
    LmH = L - hessian;
  }
  Rcpp::Rcout << "         factor=" << factor << ", min eval=" << arma::eig_sym(LmH).min() << "\n";

  this->La = La*factor;
  this->LB = LB*factor;
}

double Model::lambda_max(Data &data){
  this->penalty.lambda = 1e10;
  this->fit(data);
  return this->penalty.lambda_max(this->gB);
}

void Model::fit(Data &data){
  // Rcpp::Rcout << "1\n";
  this->update_precision(data);
  // Rcpp::Rcout << "2\n";
  this->update_mean(data);
  // Rcpp::Rcout << "3\n";
  double quad_term = this->quadratic_term(data);
  // Rcpp::Rcout << "4\n";
  double ld_scaling = this->logdet_scaling(data);
  // Rcpp::Rcout << "5\n";
  double ld_precision = this->logdet_precision(data);
  // Rcpp::Rcout << "6\n";
  double ld_dispersion = this->logdet_dispersion(data);
  // Rcpp::Rcout << "7\n";
  double llk = this->quasi_log_likelihood(quad_term, ld_scaling + ld_dispersion - ld_precision);
  // Rcpp::Rcout << "8\n";
  double llk_old = llk;
  for(uint round=0; round<this->control.max_rounds; round++){
    // mean update inner loop
    double quad_term_old = quad_term;
    for(uint iter=0; iter<this->control.max_iter; iter++){
      // Rcpp::Rcout << "9\n";
      this->update_gradients(data);
      // Rcpp::Rcout << "10\n";
      this->proximal_gradient_step();
      // Rcpp::Rcout << "11\n";
      this->update_mean(data);
      // Rcpp::Rcout << "12\n";
      quad_term = this->quadratic_term(data);
      Rcpp::Rcout << "         " << round << "." << "M" << "." << iter << ": obj=" << quad_term << "\n";
      if(fabs(quad_term - quad_term_old) / fabs(quad_term_old)< this->control.rel_tol) break;
      quad_term_old = quad_term;
    }
    // variance update
    // Rcpp::Rcout << "13\n";
    this->workingCovariance.update_parameters(data.sr);
    // Rcpp::Rcout << "14\n";
    this->update_precision(data);
    // Rcpp::Rcout << "15\n";
    quad_term = this->quadratic_term(data);
    // Rcpp::Rcout << "16\n";
    this->family.dispersion = quad_term / data.n;
    // Rcpp::Rcout << "17\n";
    ld_scaling = this->logdet_scaling(data);
    // Rcpp::Rcout << "18\n";
    ld_precision = this->logdet_precision(data);
    // Rcpp::Rcout << "19\n";
    ld_dispersion = this->logdet_dispersion(data);
    // Rcpp::Rcout << "20\n";
    llk = this->quasi_log_likelihood(quad_term, ld_scaling + ld_dispersion - ld_precision);
    // Rcpp::Rcout << "21\n";
    Rcpp::Rcout << "         " << round << " llk=" << llk << "\n";
    if(fabs(llk - llk_old) / fabs(llk_old) < this->control.rel_tol) break;
    llk_old = llk;
  }
  this->results["llk"] = llk;
  this->results["rss"] = quad_term;
  this->prepare_results(data);
}

void Model::prepare_results(const Data &data){
  this->family.add_to_results(this->results);
  this->penalty.add_to_results(this->results);
  this->linkFunction.add_to_results(this->results);
  this->workingCovariance.add_to_results(this->results);
  this->kernel.add_to_results(this->results);
  this->prepare_ics(data);
}

void Model::prepare_ics(const Data &data){
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
