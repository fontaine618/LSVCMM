#' Title
#'
#' @param data A data frame containing the variables in the model.
#' @param response The name of the response variable in \code{data}.
#' @param subject The name of the subject variable in \code{data}.
#' @param time The name of the time variable in \code{data}.
#' @param vc_covariates The names of the varying coefficient covariates in \code{data}.
#' @param nvc_covariates The names of the non-varying coefficient covariates in \code{data}.
#' @param offset The name of the offset variable in \code{data}.
#' @param add_intercept Whether to add an intercept to the model.
#' @param estimated_time The time points at which to estimate the varying coefficients. If missing, all observed time points are used.
#' @param family A list of arguments for the response distribution. See \code{\link{family_args}}.
#' @param kernel A list of arguments for the kernel. See \code{\link{kernel_args}}.
#' @param penalty A list of arguments for the penalty. See \code{\link{penalty_args}}.
#' @param working_covariance A list of arguments for the working covariance. See \code{\link{working_covariance_args}}.
#' @param control A list of arguments for the control parameters. See \code{\link{control_args}}.
#' @param n_samples Number of bootstrap samples to compute
#'
#' @return A list containing the following elements:
#' \describe{
#' \item{family}{A list of arguments for the response distribution. See \code{\link{family_args}}.}
#' \item{kernel}{A list of arguments for the kernel. See \code{\link{kernel_args}}.}
#' \item{penalty}{A list of arguments for the penalty. See \code{\link{penalty_args}}.}
#' \item{working_covariance}{A list of arguments for the working covariance. See \code{\link{working_covariance_args}}.}
#' \item{control}{A list of arguments for the control parameters. See \code{\link{control_args}}.}
#' \item{results}{A data frame containing the results of the optimization. Each row is a model resulting from a particular tuning parameter combination.}
#' \item{nvc_boot}{A matrix containing the estimated non-varying coefficient path of dimension (\code{p_u}, \code{n_samples}).}
#' \item{vc_boot}{An array containing the estimated varying coefficient path of dimension (\code{p_x}, \code{n_timepoints}, \code{n_samples}).}
#' \item{scaled_time}{The time points at which the varying coefficients were estimated.}
#' \item{unscaled_time}{The corresponding values in the original scale.}
#' \item{range_time}{The range of the time points.}
#' \item{full_model}{The model fitted on all data.}
#' \item{nvc}{A vector containing the estimated non-varying coefficient.}
#' \item{vc}{An matrix containing the estimated varying coefficient.}
#' }
#'
#' @import RcppProgress
#' @export
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
