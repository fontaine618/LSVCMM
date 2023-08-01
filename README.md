
<!-- README.md is generated from README.Rmd. Please edit that file -->

Last updated: *Aug-01-2023*

# LSVCMM

<!-- badges: start -->
[![R-CMD-check](https://github.com/fontaine618/LSVCMM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fontaine618/LSVCMM/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/LSVCMM)](https://cran.r-project.org/package=LSVCMM)
<!-- badges: end -->

## Installation

You can install the development version of LSVCMM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fontaine618/LSVCMM")
```

## Reference

For details on the LSVCMM method, refer to the [Draft
manuscript](bit.ly/lsvcmm).

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(LSVCMM)
data = generate_synthetic_data()
fit = lsvcmm(
  data=data$data, 
  response="response",
  subject="subject_id",
  time="time",
  vc_covariates="group",
  kernel=list(scale=c(0.3, 0.5, 1.))
)
```

``` r
library(ggplot2)
ggplot() + 
  geom_line(
    data=fit$results,
    mapping=aes(x=penalty.lambda, y=ebich, color=as.factor(kernel.scale), group=kernel.scale)
  ) + scale_x_log10() + 
  theme_minimal() + labs(color="Kernel scale", x="Reg. parameter", y="EBIC")
```

<img src="man/figures/README-results-1.png" width="100%" />

``` r

i = which.min(fit$results$ebich)
fit$vc_path[,,i]
#>             [,1]        [,2]       [,3]        [,4]        [,5]        [,6]
#> [1,] -0.06224119 -0.07737476 -0.0875152 -0.09274494 -0.09389309 -0.09237971
#> [2,]  0.00000000  0.00000000  0.0000000  0.00000000  0.00000000  0.10615007
#>             [,7]        [,8]        [,9]       [,10]      [,11]
#> [1,] -0.08995033 -0.08835181 -0.08901332 -0.09279479 -0.0998499
#> [2,]  0.48351439  0.80429355  1.01417988  1.06646627  0.9379095
fit$results[i, ]
#>           llk      rss family.dispersion family.name
#> 227 -1273.265 1088.203          1.452874    gaussian
#>                    penalty.name penalty.alpha penalty.lambda penalty.adaptive
#> 227 adaptive_sparse_group_lasso             1       3.006175                0
#>     penalty.penalize_intercept link_function.name working_covariance.estimate
#> 227                      FALSE           identity                        TRUE
#>     working_covariance.ratio working_covariance.name kernel.name kernel.scale
#> 227                0.4226973       compound_symmetry    gaussian            1
#>     control.max_iter control.max_rounds control.rel_tol control.verbose
#> 227             1000                 50           1e-06               1
#>     control.update_method control.backtracking_fraction df df_kernel  df_logn
#> 227                  BPGD                           0.9 17  1.545455 99.81993
#>     df_logn_kernel df_max      aic     aich      bic     bich     ebic    ebich
#> 227       9.074539     22 2580.531 2549.622 2646.351 2555.605 3020.351 2589.605
```
