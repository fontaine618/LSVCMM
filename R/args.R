# ==============================================================================
DEFAULT_FAMILY_ARGS = list(response="gaussian" , link="identity")
IMPLEMENTED_FAMILIES = c("gaussian")
IMPLEMENTED_LINKS = c("identity")


#' Prepare family arguments
#'
#' @param args A list of arguments. See See *Details*.
#'
#' @details
#' The following arguments are supported:
#' \describe{
#'  \item{\code{response}}{The response distribution. Currently, only \code{"gaussian"} is supported.}
#'  \item{\code{link}}{The link function. Currently, only \code{"identity"} is supported.}
#'  \item{\code{...}}{Additional arguments. Currently ignored.}
#' }
#'
#' If an argument is not specified, the default value is used.
#'
#' @return A list of arguments
#' @export
family_args = function(args){
  out = DEFAULT_FAMILY_ARGS
  for(name in names(out)) if(!is.null(args[[name]])) out[[name]] = args[[name]]
  stopifnot(out[["response"]] %in% IMPLEMENTED_FAMILIES)
  stopifnot(out[["link"]] %in% IMPLEMENTED_LINKS)
  return(out)
}

# ==============================================================================
DEFAULT_KERNEL_ARGS = list(name="gaussian", scale=NULL, n_scale=11L, rescale_boundary=T)
IMPLEMENTED_KERNELS = c("gaussian", "epa", "linear", "triweight")


#' Prepare kernel arguments
#'
#' @param args A list of arguments. See *Details*.
#'
#' @details
#' The following arguments are supported:
#' \describe{
#' \item{\code{name}}{The kernel function. Currently, only \code{"gaussian"} is supported.}
#' \item{\code{scale}}{The scale parameter. If \code{NULL}, a logarithmic grid search is used.}
#' \item{\code{n_scale}}{The number of scale parameters to use if \code{scale} is \code{NULL} (default.)}
#' \item{\code{rescale_boundary}}{Whether to resscale the kernel at the boundary (default: TRUE.)}
#' \item{\code{...}}{Additional arguments. Currently ignored.}
#' }
#'
#' For the grid search, the range is defined to be \eqn{[L, U]} where
#' \eqn{L} is half the smallest gap in estimated time points and
#' \eqn{U} is twice the range of estimated time points.
#'
#' @return A list of arguments
#' @export
kernel_args = function(args){
  out = DEFAULT_KERNEL_ARGS
  for(name in names(out)) if(!is.null(args[[name]])) out[[name]] = args[[name]]
  if(is.null(out[["scale"]])) out[["scale"]] = numeric(0L)
  stopifnot(out[["name"]] %in% IMPLEMENTED_KERNELS)
  stopifnot(all(out[["scale"]] >= 0))
  stopifnot(out[["n_scale"]] > 0)
  return(out)
}

# ==============================================================================
DEFAULT_PENALTY_ARGS = list(name="adaptive_sparse_group_lasso", alpha=1., a=3.7,
                            lambda=NULL, n_lambda=100L,
                            lambda_factor=0.001, adaptive=0., penalize_intercept=F)
IMPLEMENTED_PENALTIES = c("adaptive_sparse_group_lasso", "sparse_group_scad", "sparse_group_mcp")


#' Prepare penalty arguments
#'
#' @param args A list of arguments. See See *Details*.
#'
#' @details
#' The following arguments are supported:
#' \describe{
#' \item{\code{name}}{The penalty function. Currently, only \code{"adaptive_sparse_group_lasso"} is supported.}
#' \item{\code{alpha}}{The mixing parameter. Must be between 0 and 1. 1 is pure Lasso (default); 0 is pure group Lasso.}
#' \item{\code{a}}{The SCAD parameter, requires >=2. Default: 3.7.}
#' \item{\code{lambda}}{The penalty parameter. If \code{NULL} (default), a logarithmic grid search is used.}
#' \item{\code{n_lambda}}{The number of penalty parameters to use if \code{lambda} is \code{NULL} (default.)}
#' \item{\code{lambda_factor}}{The factor by which the penalty parameter is reduced overall (default: \code{1e-4}.)}
#' \item{\code{adaptive}}{The adaptive parameter. Must be non-negative (default: \code{0.}.)}
#' \item{\code{penalize_intercept}}{Whether to penalize the intercept. Must be \code{TRUE} (default) or \code{FALSE}.}
#' \item{\code{...}}{Additional arguments. Currently ignored.}
#' }
#'
#' For the grid search, we first obtain the minimum \eqn{\lambda_\max} for which the solution is zero.
#' Then, a logarithmically-space sequence of \eqn{\lambda} values is generated starting from  \eqn{\lambda_\max} and
#' ending at \eqn{\lambda_\max \times \texttt{lambda_factor}}.
#'
#'
#' @return A list of arguments
#' @export
penalty_args = function(args){
  out = DEFAULT_PENALTY_ARGS
  for(name in names(out)) if(!is.null(args[[name]])) out[[name]] = args[[name]]
  if(is.null(out[["lambda"]])) out[["lambda"]] = numeric(0L)
  stopifnot(out[["name"]] %in% IMPLEMENTED_PENALTIES)
  stopifnot(out[["alpha"]] >= 0 & out[["alpha"]] <= 1)
  stopifnot(out[["a"]] >= 2)
  stopifnot(all(out[["lambda"]] >= 0))
  stopifnot(out[["n_lambda"]] > 0)
  stopifnot(out[["lambda_factor"]] > 0 & out[["lambda_factor"]] < 1)
  stopifnot(out[["adaptive"]] >= 0)
  return(out)
}

# ==============================================================================
DEFAULT_WORKING_COVARIANCE_ARGS = list(name="compound_symmetry", estimate=T, ratio=1., correlation=0.5)
IMPLEMENTED_WORKING_COVARIANCES = c("compound_symmetry", "autoregressive", "independent")

#' Prepare working covariance arguments
#'
#' @param args A list of arguments. See See *Details*.
#'
#' @details
#' The following arguments are supported:
#' \describe{
#' \item{\code{name}}{The working covariance function. Currently, only \code{"compound_symmetry"} is supported.}
#' \item{\code{estimate}}{Whether to estimate the parameters. Must be \code{TRUE} (default) or \code{FALSE}.}
#' \item{\code{ratio}}{The ratio of the variance of the random effects to the variance of the error term.
#' If negative (default), the ratio is estimated to `log(n)`.}
#' \item{\code{correlation}}{The correlation for AR(1).}
#' }
#'
#' @return A list of arguments
#' @export
working_covariance_args = function(args){
  out = DEFAULT_WORKING_COVARIANCE_ARGS
  for(name in names(out)) if(!is.null(args[[name]])) out[[name]] = args[[name]]
  stopifnot(out[["name"]] %in% IMPLEMENTED_WORKING_COVARIANCES)
  return(out)
}

# ==============================================================================
DEFAULT_CONTROL_ARGS = list(max_rounds=50, max_iter=1000, rel_tol=1e-6, verbose=1, update_method="PGD",
                            backtracking_fraction=0.9)

#' Prepare control arguments
#'
#' @param args A list of arguments. See See *Details*.
#'
#' @details
#' The following arguments are supported:
#' \describe{
#' \item{\code{max_rounds}}{The maximum number of rounds of coordinate descent (default: 50.)}
#' \item{\code{max_iter}}{The maximum number of iterations per round (default: 1000.)}
#' \item{\code{rel_tol}}{The relative tolerance for convergence (default: \code{1e-6}.)}
#' \item{\code{verbose}}{The verbosity level. Must be 0, 1(default), 2, or 3.}
#' \item{\code{update_method}}{The update method. Currently, only \code{"PGD"} is allowed.}
#' \item{\code{backtracking_fraction}}{The backtracking fraction. Must be between 0 and 1 (default: \code{0.9}.)}
#' \item{\code{...}}{Additional arguments. Currently ignored.}
#' }
#'
#' @return A list of arguments
#' @export
control_args = function(args){
  out = DEFAULT_CONTROL_ARGS
  for(name in names(out)) if(!is.null(args[[name]])) out[[name]] = args[[name]]
  stopifnot(out[["max_rounds"]] > 0)
  stopifnot(out[["max_iter"]] > 0)
  stopifnot(out[["rel_tol"]] > 0)
  stopifnot(out[["verbose"]] %in% c(0, 1, 2, 3))
  stopifnot(out[["update_method"]] %in% c("PGD"))
  stopifnot(out[["backtracking_fraction"]] > 0 & out[["backtracking_fraction"]] < 1)
  return(out)
}

# ==============================================================================
time_args = function(observed_time, estimated_time){
  if(is.null(estimated_time)) estimated_time = unique(observed_time)
  estimated_time = sort(estimated_time)
  trange = range(observed_time)
  scaled_observed = (observed_time - trange[1]) / (trange[2] - trange[1])
  scaled_estimated = (estimated_time - trange[1]) / (trange[2] - trange[1])
  return(list(
    scaled=scaled_estimated,
    unscaled=estimated_time,
    scaled_observed=scaled_observed,
    range=trange
  ))
}

# ==============================================================================
data_args = function(
    data,
    response,
    subject,
    time,
    vc_covariates=NULL,
    nvc_covariates=NULL,
    offset=NULL,
    add_intercept=T
){
  n = nrow(data)
  response = matrix(data[[response]], nrow=n, ncol=1)
  subject = data[[subject]]
  subject_ids = sort(unique(subject))
  s = sapply(subject, function(ss) which.max(ss==subject_ids)) - 1
  subject = matrix(s, nrow=n, ncol=1)
  time = matrix(data[[time]], nrow=n, ncol=1)
  if(!is.null(offset)) offset = matrix(data[[offset]], nrow=n, ncol=1)
  else offset = matrix(0, nrow=n, ncol=1)
  if(!is.null(vc_covariates)) vc_covariates = as.matrix(data[, vc_covariates, drop=F])
  else vc_covariates = matrix(0, nrow=n, ncol=0)
  if(!is.null(nvc_covariates)) nvc_covariates = as.matrix(data[, nvc_covariates, drop=F])
  else nvc_covariates = matrix(0, nrow=n, ncol=0)
  if(add_intercept) vc_covariates = cbind(1, vc_covariates)
  return(list(
    response=response,
    subject=subject,
    time=time,
    offset=offset,
    vc_covariates=vc_covariates,
    nvc_covariates=nvc_covariates
  ))
}
