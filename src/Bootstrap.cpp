#include "RcppArmadillo.h"
#include "Data.h"
#include "Kernel.h"
#include "Family.h"
#include "LinkFunction.h"
#include "WorkingCovariance.h"
#include "Penalty.h"
#include "Model.h"
#include "Interpolator.h"
#include "Path.cpp"
#include "Control.cpp"
#include "Logger.cpp"
#include "uint.h"


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export()]]
Rcpp::List LSVCMM_Bootstrap(
    arma::colvec &response, // DATA
    arma::ucolvec &subject,
    arma::colvec &response_time,
    arma::mat &vcm_covariates,
    arma::mat &fixed_covariates,
    arma::colvec &offset,
    arma::colvec &weight,
    std::string family_name, // FAMILY
    double family_power,
    std::string link,
    std::string kernel_name, // KERNEL
    arma::rowvec &estimated_time,
    double kernel_scale,
    bool rescale_boundary,
    std::string penalty_name, // PENALTY
    bool penalize_intercept,
    double alpha,
    double scad_a,
    double adaptive,
    double lambda,
    std::string working_covariance, // DEPENDENCY
    bool estimate_variance_components,
    double variance_ratio, // negative => proxy, i.e. log(n)
    double correlation,
    unsigned int max_rounds, // CONTROL
    unsigned int max_iter,
    double rel_tol,
    unsigned int verbose,
    std::string update_method,
    double backtracking_fraction,
    bool two_step_estimation,
    double stepsize_factor,
    uint n_samples,
    bool resample_within_subject
){
  Control* control = new Control(
    max_rounds,
    max_iter,
    rel_tol,
    verbose,
    update_method,
    backtracking_fraction,
    two_step_estimation,
    stepsize_factor
  );

  if(control->verbose) Rcpp::Rcout << "[LSVCMM] Initializing data \n";
  Data data = Data(
    response,
    subject,
    response_time,
    vcm_covariates,
    fixed_covariates,
    offset,
    weight
  );

  if(control->verbose) Rcpp::Rcout << "[LSVCMM] Initializing interpolator \n";
  estimated_time = arma::sort(estimated_time);
  Interpolator interpolator = Interpolator(estimated_time);
  data.I = interpolator.interpolator_matrix(data.t);

  if(control->verbose) Rcpp::Rcout << "[LSVCMM] Initializing kernel. kernel_name=" << kernel_name << "\n";
  Kernel* kernel = Kernel::Create(kernel_name, rescale_boundary, estimated_time);

  if(control->verbose) Rcpp::Rcout << "[LSVCMM] Initializing family. family_name=" << family_name << "\n";
  Family* family = Family::Create(family_name, family_power);

  if(control->verbose) Rcpp::Rcout << "[LSVCMM] Initializing link function. link_function=" << link << "\n";
  LinkFunction* link_function = LinkFunction::Create(link);

  if(control->verbose) Rcpp::Rcout << "[LSVCMM] Initializing penalty. penalty_name=" << penalty_name << "\n";
  Penalty* penalty = Penalty::Create(penalty_name, 0., alpha, adaptive, scad_a, penalize_intercept);

  if(control->verbose) Rcpp::Rcout << "[LSVCMM] Initializing working covariance. working_covariance=" << working_covariance << "\n";
  if(variance_ratio < 0.) variance_ratio = log(data.n);
  WorkingCovariance* working_cov = WorkingCovariance::Create(working_covariance, variance_ratio, correlation, estimate_variance_components);
  data.P = working_cov->compute_precision(data.t);

  if(control->verbose) Rcpp::Rcout << "[LSVCMM] Initializing model \n";
  Model* model = new Model(
    data.px,
    data.pu,
    estimated_time,
    penalty,
    working_cov,
    link_function,
    family,
    kernel,
    control
  );

  if(control->verbose) Rcpp::Rcout << "[LSVCMM] Estimating full model \n";
  model->a.zeros();
  model->B.zeros();
  model->kernel->update_scale(kernel_scale);
  data.W = model->kernel->eval(data.t);
  Rcpp::Rcout << "         Computing Lipschitz constants\n";
  model->prepare_stepsize(data);
  Rcpp::Rcout << "         Initialization to unpenalized independent model\n";
  // need to set to independent kernel with no penalty
  double old_vr = model->workingCovariance->variance_ratio;
  model->workingCovariance->variance_ratio = 0.;
  bool estimate_parameters = model->workingCovariance->estimate_parameters;
  model->workingCovariance->estimate_parameters = false;
  model->penalty->lambda = 0.;
  model->penalty->update_weights(); // to propagate
  model->fit(data);
  Rcpp::Rcout << "         Preparing penalty weights\n";
  model->penalty->update_B0(model->B);
  // reset covariance
  model->workingCovariance->variance_ratio = old_vr;
  model->workingCovariance->estimate_parameters = estimate_parameters;
  Rcpp::Rcout << "         Initialize covariance parameters\n";
  model->update_precision(data);
  model->workingCovariance->update_parameters(
      data.sr, data.t, data.P,
      model->family->dispersion, model->logger,
      0, model->control
  );
  if(model->control->two_step_estimation){
    Rcpp::Rcout << "         Two-step estimation: fix covariance parameters\n";
    model->workingCovariance->estimate_parameters = false;
  }
  model->penalty->lambda = lambda;
  model->penalty->update_weights();
  model->fit(data);
  Rcpp::List full_model = model->save();
  model->logger->reset();

  if(control->verbose) Rcpp::Rcout << "[LSVCMM] Starting bootstrap \n";
  arma::mat B = model->B;
  arma::colvec a = model->a;
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  Rcpp::List models(n_samples);
  Progress pbar(n_samples);
  for(uint b=0; b<n_samples; b++){
    model->B = B;
    model->a = a;
    // set seed
    set_seed_r(b);
    Data rdata = data.resample(resample_within_subject);
    model->fit(rdata);
    models[b] = model->save();
    model->logger->reset();
    pbar.increment();
  }
  if(control->verbose) Rcpp::Rcout << "[LSVCMM] DONE \n";

  return Rcpp::List::create(
    Rcpp::Named("boot", models),
    Rcpp::Named("model", full_model)
  );
}
