#include "RcppArmadillo.h"
//[[Rcpp::depends(RcppArmadillo)]]

#ifndef Data_hpp
#define Data_hpp

class Data {

public:

  // Input data
  std::vector<arma::mat> X, U;
  std::vector<arma::colvec> y, t, o;
  uint px, pu, n, N;
  // Computed later
  std::vector<arma::mat> P, W, I;
  // storage during optimization
  // respectively: linear predict, mean, sd(mean), residuals, scaled residuals, scaled adjusted residuals
  std::vector<arma::colvec>  lp, m, s, r, sr, sPsr ;
  // bookkeeping for CV
  arma::uvec foldid;

  Data(
    const arma::colvec & response,
    const arma::ucolvec & subject,
    const arma::colvec & response_time,
    const arma::mat & vcm_covariates,
    const arma::mat & fixed_covariates,
    const arma::colvec & offset
  );

  // Helper to subset or resample
  Data(
    const std::vector<arma::colvec> & y,
    const std::vector<arma::colvec> & t,
    const std::vector<arma::colvec> & o,
    const std::vector<arma::mat> & X,
    const std::vector<arma::mat> & U,
    const std::vector<arma::mat> & P,
    const std::vector<arma::mat> & I,
    const std::vector<arma::mat> & W
  );

  void prepare_folds(uint nfolds);

  Data get(arma::uvec ids);

  Data get_fold(uint fold);

  Data get_other_folds(uint fold);

  Data resample();

};

#endif
