#' Generate synthetic data
#'
#' @param n_subjects Number of subjects (default 100). Half of them are in group 0 and the other half in group 1 (on average).
#' @param n_timepoints Number of time points at which the data is sampled spread regularly on the \eqn{[0,1]} interval (default 11).
#' @param prop_observed Proportion of observed time points (default 0.7).
#' @param observation_variance Variance of the observation error (default 1).
#' @param random_effect_variance_ratio Ratio of the variance of the random effect to the variance of the observation error (default 1).
#' @param random_effect_ar1_correlation Correlation between two random effects at two consecutive timepoints (default 0.9).
#' @param effect_size Effect size of the group (default 1).
#' @param seed Seed for the random number generator (default 1).
#'
#' @examples
#' instance = generate_synthetic_data()
#'
#' @details
#' The data is generated according to the following model:
#' \deqn{y_{ij} = \beta_0(t_{ij}) + \beta_1(t_{ij}) g_i + \theta_{ij} + \epsilon_{ij}}
#' where \eqn{y_{ij}} is the response of subject \eqn{i} at time \eqn{j} at time \eqn{(t_{ij})},
#' \eqn{g_i} is the group of subject \eqn{i},
#' \eqn{\theta_{ij}} is the random effect of subject \eqn{i} at time \eqn{j},
#' \eqn{\epsilon_{ij}} is the observation error of subject \eqn{i} at time \eqn{j},
#' \eqn{\beta_0(\cdot)} is the intercept function, and \eqn{\beta_1(\cdot)} is the effect size.
#' The random effects are generated according to a multivariate normal distribution
#' with mean 0 and covariance matrix \eqn{\Sigma} such that \eqn{\Sigma_{ij} = r\sigma^2\rho^{|i-j|}}
#' (\eqn{\rho} defined by \code{random_effect_ar1_correlation} and \eqn{r} by \code{random_effect_variance_ratio}).
#' The observation errors are generated according to a multivariate normal distribution
#' with mean 0 and covariance matrix \eqn{\sigma^2 I}
#' (\eqn{\sigma^2} defined by \code{observation_variance}).
#' The random effects and observation errors are independent.
#' The response is generated according to a logistic function of time and group:
#' \deqn{f(t) = \frac{1}{1 + \exp((0.6-t)20)}}
#' The observations are then sampled uniformly at a proportion \code{prop_observed} of the timepoints.
#' An imputed version of the data is also provided where the missing values are imputed using
#' the function [refund::fpca.sc].
#' The data is returned in three formats:
#' \describe{
#' \item{data}{A data frame containing the longitudinal data in long format.}
#' \item{data_wide}{A data frame containing the longitudinal data in wide format.}
#' \item{data_wide_imputed}{A data frame containing the longitudinal data in wide format with imputed missing values.}
#' }
#'
#' @return A list containing the following elements:
#' \describe{
#' \item{data}{A data frame containing the longitudinal data in long format.}
#' \item{data_wide}{A data frame containing the longitudinal data in wide format.}
#' \item{data_wide_imputed}{A data frame containing the longitudinal data in wide format with imputed missing values.}
#' \item{true_values}{A data frame containing the true values of the fixed effects.}
#' \item{times}{A vector of timepoints at which the data is sampled.}
#' \item{colnames}{A list containing the column names of the data frames.}
#' }
#'
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom dplyr filter select bind_cols arrange
#' @importFrom tidyselect starts_with
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom refund fpca.sc
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
  # to get R CMD CHECK to stop whining
  observed = NULL
  subject_id = NULL

  set.seed(seed)
  corr = random_effect_ar1_correlation^seq(0, n_timepoints-1)
  corrmat = stats::toeplitz(corr)
  thetamat = mvtnorm::rmvnorm(n_subjects, sigma=corrmat) * sqrt(random_effect_variance_ratio) * sqrt(observation_variance)
  t0 = seq(0, n_timepoints-1) / (n_timepoints-1)
  timemat = matrix(t0, n_subjects, n_timepoints, byrow=T)
  f0 = function(t) t*0.
  f1raw = function(t) 1/(1+exp((0.6-t)*20))
  f1 = function(t) ifelse(abs(f1raw(t)) < 0.1, 0, f1raw(t)) # make small values exact 0s
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

  data_wide = data_long %>% tidyr::pivot_wider(
    id_cols=c("subject_id", "group"),
    names_from="time",
    values_from="response",
    names_prefix="t",
    names_sort=TRUE
  ) %>% dplyr::arrange(subject_id)

  Y = data_wide %>% dplyr::select(dplyr::starts_with("t")) %>% as.matrix
  Yhat = refund::fpca.sc(Y=Y, argvals=t0, nbasis=4)$Yhat
  Yhat[!is.na(Y)] = Y[!is.na(Y)]

  data_wide_imputed = dplyr::bind_cols(
    data_wide %>% select(subject_id, group),
    data.frame(Yhat)
  )

  true_values = data.frame(
    time=t0,
    b0=f0(t0),
    b1=f1(t0) * effect_size
  )

  instance = list(
    data=data_long,
    data_wide=data_wide,
    data_wide_imputed=data_wide_imputed,
    true_values=true_values,
    times=t0,
    colnames=list(
      subject="subject_id",
      group="group",
      long_index="time",
      long_response="response",
      wide_times=paste0("t", t0)
    )
  )

  return(instance)
}
