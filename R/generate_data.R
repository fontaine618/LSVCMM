#' Generate synthetic data
#'
#' @param n_subjects
#' @param n_timepoints
#' @param prop_observed
#' @param observation_variance
#' @param random_effect_variance_ratio
#' @param random_effect_ar1_correlation
#' @param effect_size
#' @param seed
#'
#' @return
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom dplyr filter select
#' @importFrom magrittr %>%
generate_synthetic_data = function(
    n_subjects=100,
    n_timepoints=11,
    prop_observed=0.7,
    observation_variance=1.,
    random_effect_variance_ratio=1.,
    random_effect_ar1_correlation=0.9,
    effect_size=1,
    seed=1
){
  library(dplyr)
  library(magrittr)
  set.seed(seed)
  corr = random_effect_ar1_correlation^seq(0, n_timepoints-1)
  corrmat = stats::toeplitz(corr)
  thetamat = mvtnorm::rmvnorm(n_subjects, sigma=corrmat) * sqrt(random_effect_variance_ratio)
  t0 = seq(0, n_timepoints-1) / (n_timepoints-1)
  timemat = matrix(t0, n_subjects, n_timepoints, byrow=T)
  f0 = function(t) t*0.
  f1raw = function(t) 1/(1+exp((0.6-t)*20))
  f1 = function(t) ifelse(abs(f1raw(t)) < 0.1, 0, f1raw(t))
  group = sample(0:1, n_subjects, T)
  groupmat = matrix(group, n_subjects, n_timepoints)
  term0mat = f0(timemat)
  term1mat = f1(timemat * groupmat)
  errormat = matrix(stats::rnorm(n_timepoints*n_subjects), n_subjects, n_timepoints) * sqrt(observation_variance)
  ymat = term0mat + term1mat * effect_size + thetamat + errormat
  smat = matrix(seq(n_subjects), n_subjects, n_timepoints)
  omat = matrix(stats::runif(n_timepoints*n_subjects) < prop_observed, n_subjects, n_timepoints)

  data_full = data.frame(
    response=as.vector(ymat),
    time=as.vector(timemat),
    subject_id=as.vector(smat),
    group=as.vector(groupmat),
    observed=as.vector(omat)
  )
  data_long = data_full %>% dplyr::filter(observed) %>% dplyr::select(-observed)

  # data_wide = data_long %>% tidyr::pivot_wider(
  #   id_cols=c("subject_id", "group"),
  #   names_from="time",
  #   values_from="response",
  #   names_prefix="t",
  #   names_sort=TRUE
  # ) %>% arrange(subject_id)
  #
  # Y = data_wide %>% select(starts_with("t")) %>% as.matrix
  # Yhat = refund::fpca.sc(Y=Y, argvals=t0, nbasis=4)$Yhat
  # Yhat[!is.na(Y)] = Y[!is.na(Y)]
  #
  # data_wide_imputed = bind_cols(
  #   data_wide %>% select(subject_id, group),
  #   data.frame(Yhat)
  # )

  true_values = data.frame(
    time=t0,
    b0=f0(t0),
    b1=f1(t0) * effect_size
  )

  instance = list(
    data=data_long,
    # data_wide=data_wide,
    # data_wide_imputed=data_wide_imputed,
    true_values=true_values,
    times=t0,
    colnames=list(
      subject_id="subject_id",
      group="group",
      long_index="time",
      long_response="response"
      # wide_times=paste0("t", t0)
    )
  )

  return(instance)
}
