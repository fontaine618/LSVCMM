DEFAULT_FAMILY_ARGS = list(response="gaussian" , link="identity")
IMPLEMENTED_FAMILIES = c("gaussian")
IMPLEMENTED_LINKS = c("identity")


family_args = function(args){
  out = DEFAULT_FAMILY_ARGS
  for(name in names(out)) if(!is.null(args[[name]])) out[[name]] = args[[name]]
  stopifnot(out[["response"]] %in% IMPLEMENTED_FAMILIES)
  stopifnot(out[["link"]] %in% IMPLEMENTED_LINKS)
  return(out)
}

DEFAULT_KERNEL_ARGS = list(name="gaussian", scale=NULL, n_scale=11L)
IMPLEMENTED_KERNELS = c("gaussian")

kernel_args = function(args){
  out = DEFAULT_KERNEL_ARGS
  for(name in names(out)) if(!is.null(args[[name]])) out[[name]] = args[[name]]
  if(is.null(out[["scale"]])) out[["scale"]] = numeric(0L)
  stopifnot(out[["name"]] %in% IMPLEMENTED_KERNELS)
  stopifnot(all(out[["scale"]] >= 0))
  stopifnot(out[["n_scale"]] > 0)
  return(out)
}

DEFAULT_PENALTY_ARGS = list(name="sgl", alpha=1., lambda=NULL, n_lambda=100L,
                            lambda_factor=0.001, adaptive=0., penalize_intercept=F)
IMPLEMENTED_PENALTIES = c("sgl")

penalty_args = function(args){
  out = DEFAULT_PENALTY_ARGS
  for(name in names(out)) if(!is.null(args[[name]])) out[[name]] = args[[name]]
  if(is.null(out[["lambda"]])) out[["lambda"]] = numeric(0L)
  stopifnot(out[["name"]] %in% IMPLEMENTED_PENALTIES)
  stopifnot(out[["alpha"]] >= 0 & out[["alpha"]] <= 1)
  stopifnot(all(out[["lambda"]] >= 0))
  stopifnot(out[["n_lambda"]] > 0)
  stopifnot(out[["lambda_factor"]] > 0 & out[["lambda_factor"]] < 1)
  stopifnot(out[["adaptive"]] >= 0)
  return(out)
}

DEFAULT_WORKING_COVARIANCE_ARGS = list(name="compound_symmetry", estimate=T, ratio=-1)
IMPLEMENTED_WORKING_COVARIANCES = c("compound_symmetry")

working_covariance_args = function(args){
  out = DEFAULT_WORKING_COVARIANCE_ARGS
  for(name in names(out)) if(!is.null(args[[name]])) out[[name]] = args[[name]]
  stopifnot(out[["name"]] %in% IMPLEMENTED_WORKING_COVARIANCES)
  return(out)
}

DEFAULT_CONTROL_ARGS = list(max_rounds=50, max_iter=1000, rel_tol=1e-6, verbose=1, update_method="PGD")

control_args = function(args){
  out = DEFAULT_CONTROL_ARGS
  for(name in names(out)) if(!is.null(args[[name]])) out[[name]] = args[[name]]
  stopifnot(out[["max_rounds"]] > 0)
  stopifnot(out[["max_iter"]] > 0)
  stopifnot(out[["rel_tol"]] > 0)
  stopifnot(out[["verbose"]] %in% c(0, 1, 2, 3))
  stopifnot(out[["update_method"]] %in% c("PGD", "APGD"))
  return(out)
}

time_args = function(observed_time, estimated_time){
  if(is.null(estimated_time)) estimated_time = unique(observed_time)
  estimated_time = sort(estimated_time)
  trange = range(observed_time)
  scaled_observed = (observed_time - trange[1]) / (trange[2] - trange[1])
  scaled_estimated = (estimated_time - trange[1]) / (trange[2] - trange[1])
  return(list(
    scaled=scaled_estimated,
    unscaled=estimated_time,
    scaled_observed=scaled_observed
  ))
}

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
