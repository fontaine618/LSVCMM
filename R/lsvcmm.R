#' Title
#'
#' @param data
#' @param response
#' @param subject
#' @param time
#' @param vc_covariates
#' @param nvc_covariates
#' @param offset
#' @param add_intercept
#' @param estimated_time
#' @param family
#' @param kernel
#' @param penalty
#' @param working_covariance
#' @param control
#' @param return_models
#'
#' @return
#' @export
#'
#' @examples
lsvcmm = function(
  data,
  response,
  subject,
  time,
  vc_covariates=NULL,
  nvc_covariates=NULL,
  offset=NULL,
  add_intercept=T,
  estimated_time=NULL,
  family=DEFAULT_FAMILY_ARGS,
  kernel=DEFAULT_KERNEL_ARGS,
  penalty=DEFAULT_PENALTY_ARGS,
  working_covariance=DEFAULT_WORKING_COVARIANCE_ARGS,
  control=DEFAULT_CONTROL_ARGS,
  return_models=F
){
  family = family_args(family)
  kernel = kernel_args(kernel)
  penalty = penalty_args(penalty)
  working_covariance = working_covariance_args(working_covariance)
  control = control_args(control)
  data = data_args(data, response, subject, time, vc_covariates, nvc_covariates, offset, add_intercept)
  time = time_args(data$t, estimated_time)

  # CALL
  fit = LSVCMM::LSVCMM(
    response=data$response,
    subject=data$subject,
    response_time=time$scaled_observed,
    vcm_covariates=data$vc_covariates,
    fixed_covariates=data$nvc_covariates,
    offset=data$offset,

    family_name=family$response,
    link=family$link,

    kernel_name=kernel$name,
    estimated_time=time$scaled,
    kernel_scale=kernel$scale,
    n_kernel_scale=kernel$n_scale,

    penalize_intercept=penalty$penalize_intercept,
    alpha=penalty$alpha,
    adaptive=penalty$adaptive,
    lambda=penalty$lambda,
    lambda_factor=penalty$lambda_factor,
    n_lambda=penalty$n_lambda,

    working_covariance=working_covariance$name,
    estimate_variance_components=working_covariance$estimate,
    variance_ratio=working_covariance$ratio,

    max_rounds=control$max_rounds,
    max_iter=control$max_iter,
    rel_tol=control$rel_tol,
    verbose=control$verbose,
    update_method=control$update_method,
    backtracking_fraction=control$backtracking_fraction
  )
  fit = lapply(fit, function(model){
    class(model) = "LSVCMM"
    model
  })

  # SUMMARIZE RESULTS
  models_path = lapply(fit, function(model) data.frame(model$results)) %>% dplyr::bind_rows()

  # AGGREGATE COEFFICIENTS
  pu = ncol(data$nvc_covariates)
  a = NULL
  if(pu > 0) a = matrix(sapply(fit, function(model) model$a, simplify="matrix"), nrow=pu)
  b = sapply(fit, function(model) model$B, simplify="array")

  # RESULTS PREPARATION
  res = list(
    family=family,
    kernel=kernel,
    penalty=penalty,
    working_covariance=working_covariance,
    control=control,
    results=models_path,
    nvc_path=a,
    vc_path=b,
    scaled_time=time$scaled, # the one used during optimization
    unscaled_time=time$unscaled, # the corresponding values in the original scale
    range_time=time$range # sufficient to map to scaled time
  )
  if(return_models) res$models = fit

  class(res) = c("LSVCMM.Path")
  return(res)
}

# TODO: need to store variable names
