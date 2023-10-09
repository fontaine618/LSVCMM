#include "RcppArmadillo.h"
#include "Model.h"
#include "Data.h"
#include "uint.h"

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef Path_hpp
#define Path_hpp

class Path {

public:

  Model* model;
  arma::vec kernel_scale;
  arma::vec lambda;
  arma::uvec new_kernel;
  double lambda_factor;
  uint n_models;
  std::string mode = "grid_search";

  Path(){};

  Path(
    Model* model,
    arma::vec kernel_scale_,
    double lambda_factor_,
    uint n_lambda
  ){
    // run for a solution path at each given kernel_scale
    this->mode = "grid_search";
    this->model = model;
    this->n_models = kernel_scale_.n_elem * n_lambda;
    this->kernel_scale = arma::vec(this->n_models, arma::fill::zeros);
    this->lambda = arma::vec(this->n_models, arma::fill::zeros);
    this->new_kernel = arma::uvec(this->n_models, arma::fill::zeros);
    for (uint i = 0; i < kernel_scale_.n_elem; i++) {
      this->kernel_scale.subvec(i*n_lambda, (i+1)*n_lambda-1) = arma::vec(n_lambda, arma::fill::ones) * kernel_scale_(i);
      this->new_kernel(i*n_lambda) = 1;
    }
    this->lambda_factor = pow(lambda_factor_, 1./(n_lambda-1.));
    // Rcpp::Rcout << "[LSVCMM] Lambda factor: " << this->lambda_factor << " lambda_factor_=" <<
    //   lambda_factor_ << " n_lambda=" << n_lambda <<"\n";
  }

  Path(
    Model* model,
    arma::vec kernel_scale,
    arma::vec lambda
  ){
    // run for each pair of kernel scale and lambda
    this->mode = "path";
    this->model = model;
    this->kernel_scale = kernel_scale;
    this->lambda = lambda;

    // find out when to restart
    this->n_models = kernel_scale.n_elem;
    this->new_kernel = arma::uvec(this->n_models, arma::fill::zeros);
    this->new_kernel(0) = 1;
    for (uint i = 1; i < this->n_models; i++) {
      if (kernel_scale(i) != kernel_scale(i-1)) {
        this->new_kernel(i) = 1;
      }
    }
  }

  Rcpp::List run(Data &data){
    Rcpp::List models(this->n_models);
    for(uint m=0; m<this->n_models; m++){
      // new kernel scale
      if(this->new_kernel(m)){
        Rcpp::Rcout << "[LSVCMM] NEW KERNEL SCALE (" << this->kernel_scale(m) << ", "<< m+1 << "/" << this->n_models << ")\n";
        this->model->a.zeros();
        this->model->B.zeros();
        this->model->kernel->scale = this->kernel_scale(m);
        data.W = this->model->kernel->eval(this->model->estimated_time, data.t);

        // prepare Lipschitz constants
        Rcpp::Rcout << "         Computing Lipschitz constants\n";
        this->model->prepare_stepsize(data);

        // if(this->model->penalty->compute_penalty_weights){
        // we are not using this anymore, since initialization at the independent MLE is better in all cases.
        // anyway, there are very few cases where we don't want to do this
        // one example where this is important is adaptiveSGL with lambda=0 fixed; most other cases
        // need this anyways
          // compute adaptive weights
          Rcpp::Rcout << "         Initialization & Preparing penalty weights\n";
          // need to set to independent kernel with no penalty
          double old_vr = this->model->workingCovariance->variance_ratio;
          this->model->workingCovariance->variance_ratio = 0.;
          bool estimate_parameters = this->model->workingCovariance->estimate_parameters;
          this->model->workingCovariance->estimate_parameters = false;

          this->model->penalty->lambda = 0.;
          this->model->penalty->update_weights(); // to propagate
          this->model->fit(data);
          this->model->logger->reset();
          this->model->penalty->update_B0(this->model->B);
          // reset covariance
          this->model->workingCovariance->variance_ratio = old_vr;
          this->model->workingCovariance->estimate_parameters = estimate_parameters;
        // }

        // find largest lambda
        if(this->mode == "grid_search"){
          Rcpp::Rcout << "         Finding lambda max\n";
          this->lambda(m) = this->model->lambda_max(data);
        }
      }
      // same kernel scale
      // this->lambda.t().print("lambda for path");
      if(this->mode == "grid_search" and !this->new_kernel(m)) this->lambda(m) = this->lambda(m-1) * this->lambda_factor;
      // run
      Rcpp::Rcout << "         Fitting lambda=" << this->lambda(m) << "\n";
      this->model->penalty->lambda = this->lambda(m);
      this->model->penalty->update_weights();
      this->model->fit(data);

      models[m] = this->model->save();
      this->model->logger->reset();
    }
    return models;
  }


};

#endif
