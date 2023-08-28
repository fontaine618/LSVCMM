#include "RcppArmadillo.h"
#include "Data.h"
#include "Kernel.h"
#include "LinkFunction.cpp"
#include "Family.h"
#include "WorkingCovariance.h"
#include "Penalty.h"
#include "Model.h"
#include "Logger.cpp"
#include "Control.cpp"
#include "uint.h"

//[[Rcpp::depends(RcppArmadillo)]]

Model::Model(
  uint px,
  uint pu,
  arma::rowvec estimated_time,
  Penalty* penalty,
  WorkingCovariance* workingCovariance,
  Identity* linkFunction,
  Gaussian* family,
  Kernel* kernel,
  Control* control
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

  this->logger = new Logger();
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
  data.m = this->linkFunction->eval(data.lp);
  for(uint i=0; i<data.N; i++){
    data.r[i] = data.y[i] - data.m[i];
    data.s[i] = arma::sqrt(this->family->unit_variance(data.m[i]));
    data.sr[i] = data.r[i] / data.s[i];
    data.sPsr[i] = (data.P[i] * data.sr[i]) / data.s[i];
  }
}

void Model::update_precision(Data &data){
  data.P = this->workingCovariance->compute_precision(data.t);
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
  return data.n * log(2*arma::datum::pi*this->family->dispersion);
}

double Model::quasi_log_likelihood(double quad_term, double logdet_term){
  return -0.5 * (quad_term + logdet_term);
}

void Model::update_gradients(const Data &data){
  this->gB.zeros();
  this->ga.zeros();
  uint nt = this->B.n_cols;
  std::vector<arma::colvec> d = this->linkFunction->derivative(data.lp);
  for(uint i=0; i<data.N; i++){
    arma::colvec dsPsr = d[i] % data.sPsr[i];
    for(uint j=0; j<nt; j++){
      arma::vec wdsPsr = data.W[i].col(j) % dsPsr;
      this->gB.col(j) += data.X[i].t() * wdsPsr;
    }
    if(data.pu) this->ga += data.U[i].t() * dsPsr;
  }
  this->gB *= -1;
  this->ga *= -1;
}

void Model::proximal_gradient_step(Data &data){
  this->update_mean(data);
  this->update_gradients(data);
  this->a -= this->ga * this->ssa;
  this->B -= this->gB * this->ssB;
  this->B = this->penalty->proximal(this->B, this->ssB);
}

void Model::backtracking_proximal_gradient_step(Data &data){
  double beta = this->control->backtracking_fraction;

  arma::colvec current_a = this->a;
  arma::mat current_B = this->B;
  this->update_mean(data);
  this->update_gradients(data);
  double current_quad = this->quadratic_term(data);

  this->ssa /= beta; // 1st iteration will undo this
  this->ssB /= beta;
  arma::colvec proposed_a;
  arma::mat proposed_B;
  double LHS;
  double RHS;

  for(uint iter=0; iter<100; iter++){
    this->ssa *= beta;
    this->ssB *= beta;
    proposed_a = current_a - this->ga * this->ssa;
    proposed_B = current_B - this->gB * this->ssB;
    proposed_B = this->penalty->proximal(proposed_B, this->ssB);
    this->a = proposed_a;
    this->B = proposed_B;
    // current_B.print("current_B");
    // proposed_B.print("proposed_B");
    this->update_mean(data);
    LHS = this->quadratic_term(data);
    RHS = current_quad;
    RHS -= arma::accu(this->ga % (current_a - proposed_a));
    RHS -= arma::accu(this->gB % (current_B - proposed_B));
    RHS += arma::accu(arma::square(current_a - proposed_a)) / (2. * this->ssa);
    RHS += arma::accu(arma::square(current_B - proposed_B)) / (2. * this->ssB);
    if(this->control->verbose > 2) Rcpp::Rcout << "         Backtracking iteration: " << iter <<
                                    " LHS=" << LHS << " RHS=" << RHS <<
                                    " ssB=" << this->ssB << std::endl;
    if(LHS <= RHS + this->control->rel_tol*fabs(RHS)) break;
  }

  // NB: selected proposal already stored in a and B
}

void Model::accelerated_proximal_gradient_step(){
  Rcpp::stop("APGD not working yet");
  double pt = this->momentum;
  double newt = 0.5 * (1. + sqrt(1 + 4.*pt*pt));
  arma::colvec tmpa = this->a - this->ga / this->La;
  arma::mat tmpB = this->B - this->gB / this->LB;
  tmpB = this->penalty->proximal(tmpB, 1./this->LB);
  this->a = tmpa + (pt - 1.) * (tmpa - this->aprev) / newt;
  this->B = tmpB + (pt - 1.) * (tmpB - this->Bprev) / newt;
  this->aprev = tmpa;
  this->Bprev = tmpB;
  this->momentum = newt;
}

void Model::backtracking_accelerated_proximal_gradient_step(const Data &data){
  Rcpp::stop("BAPGD not working yet");
}

void Model::initialize(Data &data){}

arma::mat Model::hessian(const Data &data){
  // assumes mean and precision are updated

  // FIXME: not sure if this is correct for non-identity link?
  // should there be a second term?

  uint dim = data.pu + data.px * this->nt;
  arma::mat H = arma::mat(dim, dim, arma::fill::zeros);

  std::vector<arma::colvec> d = this->linkFunction->derivative(data.lp);
  for(uint i=0; i<data.N; i++){
    uint ni = data.W[i].n_rows;
    arma::colvec ds = d[i] / data.s[i];
    arma::mat dsPsd = data.P[i];
    dsPsd.each_row() %= ds.t();
    dsPsd.each_col() %= ds;
    arma::mat UwX = arma::mat(ni, dim, arma::fill::zeros);
    if(data.pu) UwX.cols(0, data.pu-1) = data.U[i];
    for(uint t=0; t<this->nt; t++){
      arma::mat wX = data.X[i];
      wX.each_col() %= data.W[i].col(t);
      UwX.cols(data.pu + t*data.px, data.pu + (t+1)*data.px - 1) = wX;
    }
    H += UwX.t() * dsPsd * UwX;
  }
  return H;
}

void Model::prepare_stepsize(Data &data){
  // NB: this does not find proper Lipshitx bounds, just heuristics,
  // but it find decent values to start from with backtracking.
  // Perhaps this is unnecessary, but it is not too expensive, I believe.
  this->update_mean(data);
  this->update_precision(data);
  arma::mat hessian = this->hessian(data);
  uint dim = data.pu + data.px*this->nt;
  arma::mat Haa;
  if(data.pu > 0) Haa = hessian.submat(0, 0, data.pu - 1, data.pu - 1);
  arma::mat HBB = hessian.submat(data.pu, data.pu, dim - 1, dim - 1);

  // Find heuristic values for La, Lb

  double La = 0.;
  if(data.pu > 0) La = arma::eig_sym(Haa).max();

  // arma::eig_sym(HBB).print("Evals HBB");
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
  // if(this->control.verbose > 1) Rcpp::Rcout << "         factor=" << factor <<
  //   ", min eval=" << arma::eig_sym(LmH).min() << "\n";

  this->La = La*factor;
  this->LB = LB*factor;
}

void Model::reset_stepsize(){
  // current hack to get PGD to work
  this->ssa = 0.1/this->La;
  this->ssB = 0.1/this->LB;
}

double Model::lambda_max(Data &data){
  this->penalty->lambda = 1e10;
  this->fit(data);
  this->logger->reset();
  return this->penalty->lambda_max(this->B, this->gB, 1./this->LB);
}

void Model::fit(Data &data){
  // // hack for the moment
  // if(this->kernel->scale > 0.5) this->control->update_method = "BPGD";
  // else this->control->update_method = "PGD";

  this->logger->reset();
  this->momentum = 1.;
  this->aprev = this->a;
  this->Bprev = this->B;
  this->update_precision(data);
  this->update_mean(data);
  double quad_term = this->quadratic_term(data);
  double penalty_term = this->penalty->eval(this->B);
  double objective = 0.5*quad_term + penalty_term;
  double objective_old = objective;
  double ld_scaling = this->logdet_scaling(data);
  double ld_precision = this->logdet_precision(data);
  double ld_dispersion = this->logdet_dispersion(data);
  double llk = this->quasi_log_likelihood(
    quad_term/this->family->dispersion,
    ld_scaling + ld_dispersion - ld_precision
  );
  double llk_old = llk;
  for(uint round=0; round<this->control->max_rounds; round++){
    this->reset_stepsize();
    double quad_term_old = quad_term;
    uint mean_iter;
    for(mean_iter=0; mean_iter<this->control->max_iter; mean_iter++){
      if(this->control->update_method == "PGD") this->proximal_gradient_step(data);
      else if (this->control->update_method == "BPGD") this->backtracking_proximal_gradient_step(data);
      else if (this->control->update_method == "APGD") this->accelerated_proximal_gradient_step();
      else if (this->control->update_method == "BAPGD") this->backtracking_accelerated_proximal_gradient_step(data);
      else Rcpp::stop("Unknown update method.");
      this->update_mean(data); // needs to be done after updating
      quad_term = this->quadratic_term(data);
      penalty_term = this->penalty->eval(this->B);
      objective = 0.5*quad_term + penalty_term;
      this->logger->add_mean_iteration_results(round, mean_iter, objective);
      if(fabs(quad_term - quad_term_old) / fabs(quad_term_old)< this->control->rel_tol){
        if(this->control->verbose > 2) Rcpp::Rcout << "         " << round << "." << "M" <<
          "." << mean_iter << ": obj=" << quad_term << "\n";
        break;
      }
      quad_term_old = quad_term;
      if(objective > objective_old){
        // Rcpp::Rcout << "         Objective increased\n";
        this->ssa *= this->control->backtracking_fraction;
        this->ssB *= this->control->backtracking_fraction;
      };
      objective_old = objective;
    }
    // variance update: mean should already be up to date
    uint variance_iter = this->workingCovariance->update_parameters(
      data.sr, data.t, data.P, this->family->dispersion, this->logger, round, this->control
    );
    // need to recompute precision matrices with new parameters
    if(this->workingCovariance->estimate_parameters) this->update_precision(data);
    quad_term = this->quadratic_term(data);
    this->family->dispersion = quad_term / data.n;
    ld_scaling = this->logdet_scaling(data);
    ld_precision = this->logdet_precision(data);
    ld_dispersion = this->logdet_dispersion(data);
    llk = this->quasi_log_likelihood(
      quad_term/this->family->dispersion,
      ld_scaling + ld_dispersion - ld_precision
    );
    if(this->control->verbose > 1) Rcpp::Rcout << "         " << round << " llk=" << llk << "\n";
    this->logger->add_round_results(round+1, mean_iter+1, variance_iter+1, llk);
    if(fabs(llk - llk_old) / fabs(llk_old) < this->control->rel_tol) break;
    llk_old = llk;
  }
  this->results["llk"] = llk;
  this->results["rss"] = quad_term;
  this->prepare_results(data);
}

void Model::prepare_results(const Data &data){
  this->family->add_to_results(this->results);
  this->penalty->add_to_results(this->results);
  this->linkFunction->add_to_results(this->results);
  this->workingCovariance->add_to_results(this->results);
  this->kernel->add_to_results(this->results);
  this->control->add_to_results(this->results);
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
  double mult = this->kernel->eval0 / (this->kernel->scale * df.n_elem);
  mult = fmin(1., mult); // when h is small, we need to avoid going above 1
  mult = fmax(1. / df.n_elem, mult); // when h is large, we need to keep 1 df
  arma::rowvec df_kernel = df * mult;

  results["df"] = arma::accu(df) + data.pu;
  results["df_kernel"] = arma::accu(df_kernel) + data.pu;
  results["df_logn"] = arma::dot(df, arma::log(nhat)) + data.pu * log(data.n);
  results["df_logn_kernel"] = arma::dot(df_kernel, arma::log(nhat)) + data.pu * log(data.n);
  results["df_max"] = data.px * this->nt + data.pu;

  results["aic"] = -2 * (double)results["llk"] + 2 * (double)results["df"];
  results["aich"] = -2 * (double)results["llk"] + 2 * (double)results["df_kernel"];
  results["bic"] = -2 * (double)results["llk"] + (double)results["df_logn"];
  results["bich"] = -2 * (double)results["llk"] + (double)results["df_logn_kernel"];
  results["ebic"] = -2 * (double)results["llk"] + (double)results["df_logn"] +
    (double)results["df"] * log((double)results["df_max"]);
  results["ebich"] = -2 * (double)results["llk"] + (double)results["df_logn_kernel"] +
    (double)results["df_kernel"] * log((double)results["df_max"]);
}

Rcpp::List Model::save(){
  return Rcpp::List::create(
    Rcpp::Named("a", this->a),
    Rcpp::Named("B", this->B),
    Rcpp::Named("results", clone(this->results)),
    Rcpp::Named("log", clone(this->logger->to_rcpp()))
  );
}

void Model::set(Rcpp::List& model){
  // TODO
  // a, B, kscale, lambda, dispersion, variatince parameters
}
