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
#' @param n_samples
#'
#' @return
#' @export
#'
#' @examples
lsvcmm.boot = function(
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
    n_samples=1000
){
  family = family_args(family)
  kernel = kernel_args(kernel)
  if(length(kernel$scale) != 1) stop("lsvcmm.boot can only take one lamba value")
  penalty = penalty_args(penalty)
  if(length(penalty$lambda) != 1) stop("lsvcmm.boot can only take one lamba value")
  working_covariance = working_covariance_args(working_covariance)
  control = control_args(control)
  data = data_args(data, response, subject, time, vc_covariates, nvc_covariates, offset, add_intercept)
  time = time_args(data$t, estimated_time)
  boot = LSVCMM_Bootstrap(
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
    rescale_boundary=kernel$rescale_boundary,

    penalty_name=penalty$name,
    penalize_intercept=penalty$penalize_intercept,
    alpha=penalty$alpha,
    scad_a=penalty$a,
    adaptive=penalty$adaptive,
    lambda=penalty$lambda,

    working_covariance=working_covariance$name,
    estimate_variance_components=working_covariance$estimate,
    correlation=working_covariance$correlation,
    variance_ratio=working_covariance$ratio,

    max_rounds=control$max_rounds,
    max_iter=control$max_iter,
    rel_tol=control$rel_tol,
    verbose=control$verbose,
    update_method=control$update_method,
    backtracking_fraction=control$backtracking_fraction,
    two_step_estimation=control$two_step,
    stepsize_factor=control$stepsize_factor,

    n_samples=n_samples
  )
  # return(boot)
  # SUMMARIZE RESULTS
  results =lapply(boot$boot, function(model) data.frame(model$results)) %>% dplyr::bind_rows()

  # AGGREGATE COEFFICIENTS
  pu = ncol(data$nvc_covariates)
  a = NULL
  if(pu > 0) a = matrix(sapply(boot$boot, function(model) model$a, simplify="matrix"), nrow=pu)
  b = sapply(boot$boot, function(model) model$B, simplify="array")  # px*nt*_samples

  # RESULTS PREPARATION
  res = list(
    family=family,
    kernel=kernel,
    penalty=penalty,
    working_covariance=working_covariance,
    control=control,
    results=results,
    nvc_boot=a,
    vc_boot=b,
    scaled_time=time$scaled, # the one used during optimization
    unscaled_time=time$unscaled, # the corresponding values in the original scale
    range_time=time$range, # sufficient to map to scaled time
    full_model=boot$model,
    vc=boot$model$B,
    nvc=boot$model$a
  )
  class(res) = c("LSVCMM.Boot")
  return(res)

}
