
<!-- README.md is generated from README.Rmd. Please edit that file -->

Last updated: *Aug-12-2023*

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
#>             [,1]        [,2]       [,3]       [,4]       [,5]       [,6]
#> [1,] -0.05392229 -0.07821785 -0.0979816 -0.1124308 -0.1214756 -0.1256716
#> [2,]  0.00000000  0.00000000  0.0000000  0.0000000  0.0277699  0.1815235
#>            [,7]       [,8]       [,9]      [,10]      [,11]
#> [1,] -0.1260648 -0.1239628 -0.1206807 -0.1173126 -0.1145710
#> [2,]  0.5833630  0.9045200  1.0915827  1.1004405  0.9117964
data$true_values
#>    time b0        b1
#> 1   0.0  0 0.0000000
#> 2   0.1  0 0.0000000
#> 3   0.2  0 0.0000000
#> 4   0.3  0 0.0000000
#> 5   0.4  0 0.0000000
#> 6   0.5  0 0.1192029
#> 7   0.6  0 0.5000000
#> 8   0.7  0 0.8807971
#> 9   0.8  0 0.9820138
#> 10  0.9  0 0.9975274
#> 11  1.0  0 0.9996646
fit$results[i, ]
#>           llk      rss family.dispersion family.name
#> 234 -1272.443 1085.661          1.449481    gaussian
#>                    penalty.name penalty.alpha penalty.lambda penalty.adaptive
#> 234 adaptive_sparse_group_lasso             1       1.844561                0
#>     penalty.penalize_intercept link_function.name working_covariance.estimate
#> 234                      FALSE           identity                        TRUE
#>     working_covariance.ratio working_covariance.name kernel.name kernel.scale
#> 234                0.4232983       compound_symmetry    gaussian            1
#>     control.max_iter control.max_rounds control.rel_tol control.verbose
#> 234             1000                 50           1e-06               1
#>     control.update_method control.backtracking_fraction df df_kernel  df_logn
#> 234                  BPGD                           0.9 18  1.636364 105.7625
#>     df_logn_kernel df_max      aic     aich      bic     bich     ebic    ebich
#> 234       9.614768     22 2580.886 2548.159 2650.649 2554.501 2706.288 2559.559
```
