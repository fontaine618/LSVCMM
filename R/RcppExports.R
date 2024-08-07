# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

LSVCMM_Bootstrap <- function(response, subject, response_time, vcm_covariates, fixed_covariates, offset, weight, family_name, family_power, link, kernel_name, estimated_time, kernel_scale, rescale_boundary, penalty_name, penalize_intercept, alpha, scad_a, adaptive, lambda, working_covariance, estimate_variance_components, variance_ratio, correlation, max_rounds, max_iter, rel_tol, verbose, update_method, backtracking_fraction, two_step_estimation, stepsize_factor, n_samples, resample_within_subject) {
    .Call(`_LSVCMM_LSVCMM_Bootstrap`, response, subject, response_time, vcm_covariates, fixed_covariates, offset, weight, family_name, family_power, link, kernel_name, estimated_time, kernel_scale, rescale_boundary, penalty_name, penalize_intercept, alpha, scad_a, adaptive, lambda, working_covariance, estimate_variance_components, variance_ratio, correlation, max_rounds, max_iter, rel_tol, verbose, update_method, backtracking_fraction, two_step_estimation, stepsize_factor, n_samples, resample_within_subject)
}

LSVCMM <- function(response, subject, response_time, vcm_covariates, fixed_covariates, offset, weight, family_name, family_power, link, kernel_name, estimated_time, kernel_scale, rescale_boundary, n_kernel_scale, penalty_name, penalize_intercept, alpha, scad_a, adaptive, lambda, lambda_factor, n_lambda, working_covariance, estimate_variance_components, variance_ratio, correlation, max_rounds, max_iter, rel_tol, verbose, update_method, backtracking_fraction, two_step_estimation, stepsize_factor) {
    .Call(`_LSVCMM_LSVCMM`, response, subject, response_time, vcm_covariates, fixed_covariates, offset, weight, family_name, family_power, link, kernel_name, estimated_time, kernel_scale, rescale_boundary, n_kernel_scale, penalty_name, penalize_intercept, alpha, scad_a, adaptive, lambda, lambda_factor, n_lambda, working_covariance, estimate_variance_components, variance_ratio, correlation, max_rounds, max_iter, rel_tol, verbose, update_method, backtracking_fraction, two_step_estimation, stepsize_factor)
}

