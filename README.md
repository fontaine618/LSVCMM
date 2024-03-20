

# LSVCMM

## Installation

From the command line

```
R CMD INSTALL <path to LSVCMM package folder>
```


## Example

``` r
library(LSVCMM)
instance = generate_synthetic_data()
fit = lsvcmm(
  data=instance$data, 
  response="response",
  subject="subject_id",
  time="time",
  vc_covariates="group",
  kernel=list(scale=c(0.1, 0.2, 0.3)),
  penalty=list(adaptive=0.5, penalize_intercept=T)
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

``` r
i = which.min(fit$results$ebich)
t(fit$vc_path[,,i])
instance$true_values
fit$results[i, ]
```
