#include "RcppArmadillo.h"
#include "Data.hpp"
#include "Kernel.hpp"
#include "GaussianKernel.hpp"
#include "LinkFunction.cpp"
#include "Family.hpp"
#include "WorkingCovariance.hpp"
#include "Penalty.hpp"
#include "Model.hpp"
#include "Interpolator.hpp"
#include "Path.cpp"

// [[Rcpp::depends(RcppArmadillo)]]

// '@export
// [[Rcpp::export()]]
Rcpp::List LSVCMM(
    arma::colvec response, // DATA
    arma::ucolvec subject,
    arma::colvec response_time,
    arma::mat vcm_covariates,
    arma::mat fixed_covariates,
    arma::colvec offset,
    std::string family_name, // FAMILY
    std::string link,
    std::string kernel_name, // KERNEL
    arma::rowvec estimated_time,
    arma::vec kernel_scale, // if this has 0 elements, we use a grid
    double kernel_scale_factor,
    uint n_kernel_scale,
    bool penalize_intercept, // PENALTY
    double alpha,
    double adaptive,
    arma::vec lambda, // if this has 0 elements, we use a grid
    double lambda_factor,
    uint n_lambda,
    std::string working_covariance, // DEPENDENCY
    bool estimate_variance_components,
    double variance_ratio, // negative => proxy, i.e. log(n)
    uint max_rounds, // CONTROL
    uint max_iter,
    double rel_tol,
    uint verbose
){
  Control control = Control(max_rounds, max_iter, rel_tol, verbose);

  if(control.verbose) Rcpp::Rcout << "[LSVCMM] Initializing data \n";
  Data data = Data(
    response,
    subject,
    response_time,
    vcm_covariates,
    fixed_covariates,
    offset
  );
  estimated_time = arma::sort(estimated_time);
  Interpolator interpolator = Interpolator(estimated_time);
  data.I = interpolator.interpolator_matrix(data.t);

  if(control.verbose) Rcpp::Rcout << "[LSVCMM] Initializing kernel \n";
  Kernel kernel = Kernel(); // defaults to uniform
  if(kernel_name == "gaussian") kernel = GaussianKernel();
  data.W = kernel.eval(estimated_time, data.t);

  if(control.verbose) Rcpp::Rcout << "[LSVCMM] Initializing family \n";
  Family family = Family();
  if(family_name == "gaussian") family = Gaussian();
  LinkFunction link_function = LinkFunction();
  if(link == "identity") link_function = Identity();

  if(control.verbose) Rcpp::Rcout << "[LSVCMM] Initializing penalty \n";
  Penalty penalty = Penalty(0., alpha, adaptive, penalize_intercept);

  if(control.verbose) Rcpp::Rcout << "[LSVCMM] Initializing working covariance \n";
  WorkingCovariance working_cov = WorkingCovariance();
  if(working_covariance == "independence") working_cov = Independence();
  if(working_covariance == "compound_symmetry") {
    if(!estimate_variance_components and variance_ratio < 0.) variance_ratio = log(data.n);
    working_cov = CompoundSymmetry(variance_ratio, estimate_variance_components);
  }
  data.P = working_cov.compute_precision(data.t);

  if(control.verbose) Rcpp::Rcout << "[LSVCMM] Initializing model \n";
  Model model = Model(
    estimated_time,
    penalty,
    working_cov,
    link_function,
    family,
    control,
    kernel
  );

  if(control.verbose) Rcpp::Rcout << "[LSVCMM] Preparing grid search \n";
  if(kernel_scale.n_elem == 0){
    double range = estimated_time.max() - estimated_time.min();
    double min_gap = arma::diff(estimated_time).min();
    kernel_scale = arma::logspace<arma::rowvec>(log10(min_gap/2.), log10(range), n_kernel_scale);
    kernel_scale.print();
  }

  Path path;
  if(lambda.n_elem == 0){
    path = Path(model, kernel_scale, lambda_factor, n_lambda);
  }else{
    path = Path(model, kernel_scale, lambda);
  }

  if(control.verbose) Rcpp::Rcout << "[LSVCMM] RUN \n";
  Rcpp::List models = path.run(data);
  if(control.verbose) Rcpp::Rcout << "[LSVCMM] DONE \n";

  return models;
}
