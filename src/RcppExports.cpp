// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// LSVCMM
Rcpp::List LSVCMM(arma::colvec& response, arma::ucolvec& subject, arma::colvec& response_time, arma::mat& vcm_covariates, arma::mat& fixed_covariates, arma::colvec& offset, std::string family_name, std::string link, std::string kernel_name, arma::rowvec& estimated_time, arma::vec& kernel_scale, unsigned int n_kernel_scale, bool penalize_intercept, double alpha, double adaptive, arma::vec& lambda, double lambda_factor, unsigned int n_lambda, std::string working_covariance, bool estimate_variance_components, double variance_ratio, unsigned int max_rounds, unsigned int max_iter, double rel_tol, unsigned int verbose, std::string update_method, double backtracking_fraction);
RcppExport SEXP _LSVCMM_LSVCMM(SEXP responseSEXP, SEXP subjectSEXP, SEXP response_timeSEXP, SEXP vcm_covariatesSEXP, SEXP fixed_covariatesSEXP, SEXP offsetSEXP, SEXP family_nameSEXP, SEXP linkSEXP, SEXP kernel_nameSEXP, SEXP estimated_timeSEXP, SEXP kernel_scaleSEXP, SEXP n_kernel_scaleSEXP, SEXP penalize_interceptSEXP, SEXP alphaSEXP, SEXP adaptiveSEXP, SEXP lambdaSEXP, SEXP lambda_factorSEXP, SEXP n_lambdaSEXP, SEXP working_covarianceSEXP, SEXP estimate_variance_componentsSEXP, SEXP variance_ratioSEXP, SEXP max_roundsSEXP, SEXP max_iterSEXP, SEXP rel_tolSEXP, SEXP verboseSEXP, SEXP update_methodSEXP, SEXP backtracking_fractionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< arma::ucolvec& >::type subject(subjectSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type response_time(response_timeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type vcm_covariates(vcm_covariatesSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type fixed_covariates(fixed_covariatesSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< std::string >::type family_name(family_nameSEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type estimated_time(estimated_timeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type kernel_scale(kernel_scaleSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_kernel_scale(n_kernel_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type penalize_intercept(penalize_interceptSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_factor(lambda_factorSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_lambda(n_lambdaSEXP);
    Rcpp::traits::input_parameter< std::string >::type working_covariance(working_covarianceSEXP);
    Rcpp::traits::input_parameter< bool >::type estimate_variance_components(estimate_variance_componentsSEXP);
    Rcpp::traits::input_parameter< double >::type variance_ratio(variance_ratioSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_rounds(max_roundsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type rel_tol(rel_tolSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< std::string >::type update_method(update_methodSEXP);
    Rcpp::traits::input_parameter< double >::type backtracking_fraction(backtracking_fractionSEXP);
    rcpp_result_gen = Rcpp::wrap(LSVCMM(response, subject, response_time, vcm_covariates, fixed_covariates, offset, family_name, link, kernel_name, estimated_time, kernel_scale, n_kernel_scale, penalize_intercept, alpha, adaptive, lambda, lambda_factor, n_lambda, working_covariance, estimate_variance_components, variance_ratio, max_rounds, max_iter, rel_tol, verbose, update_method, backtracking_fraction));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LSVCMM_LSVCMM", (DL_FUNC) &_LSVCMM_LSVCMM, 27},
    {NULL, NULL, 0}
};

RcppExport void R_init_LSVCMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}