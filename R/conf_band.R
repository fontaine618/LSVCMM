#' @importFrom dplyr bind_rows
confidence_band = function(
    obj,
    level=0.95,
    method=c("quantile"),
    var=2
){
  if(!("LSVCMM.Boot" %in% class(obj))) stop("obj must be of class LSVCMM.Boot, i.e., the output of a lsvcmm.boot call")
  stopifnot(0<level & level<1)
  level = 1-level
  method = match.arg(method)
  samples = obj$vc_boot[var,,,drop=F]  # nvars x nt x B

  out = switch(
    method,
    "quantile"=quantile_confidence_band(samples, level)
  )

  prop_in_band = proportion_samples_in_band(out$L, out$U, samples)
  pval = pointwise_pvalues(samples)
  nvars = length(var)
  nt = dim(samples)[2]
  pval.adj = matrix(stats::p.adjust(pval, method="BH"), nvars, nt)
  df = lapply(seq_along(var), function(j){
    data.frame(
      L=out$L[j,],
      U=out$U[j,],
      pointwise_confidence=1-2*out$p,
      estimate=obj$vc[j, ],
      mean=apply(samples[j,,], 1, mean),
      median=apply(samples[j,,], 1, stats::median),
      prop0=apply(samples[j,,], 1, function(x) mean(x==0)),
      estimated_time=obj$unscaled_time,
      level=prop_in_band,
      excludes_zero=(out$L[j,]>0) | (out$U[j,]<0),
      pval=pval[j,],
      pval.adj=pval.adj[j,],
      var=var[j]
    )
  }) %>% dplyr::bind_rows()

  return(df)
}

#' @importFrom pracma logseq
quantile_confidence_band = function(samples, level){
  B = dim(samples)[3]
  nc = dim(samples)[2] * dim(samples)[1]
  minp = level / (2*nc)
  maxp = level / 2
  # ps = c(seq(minp, maxp, 1/(B+1)), maxp)
  ps = pracma::logseq(minp, maxp, 200)
  prop = sapply(ps, function(p){
    out = pointwise_quantile_confidence_band(p, samples)
    proportion_samples_in_band(out$L, out$U, samples)
  })
  valid = prop >= 1-level
  which = length(valid) + 1 - which.max(rev(valid))
  p = ps[which]
  out = pointwise_quantile_confidence_band(p, samples)
  return(out)
}


pointwise_quantile_confidence_band = function(p, samples){
  qs = apply(samples, 1:2, stats::quantile, c(p, 1-p)) # 2 x len(var) x nt
  L = qs[1,,]
  U = qs[2,,]
  if(dim(samples)[1]==1){
    nt = dim(samples)[2]
    L = matrix(L, 1, nt)
    U = matrix(U, 1, nt)
  }
  return(list(L=L, U=U, p=p))
}

proportion_samples_in_band = function(L, U, samples){
  nvars = dim(samples)[1]
  between = sapply(seq(nvars), function(j){
    aL = apply(samples[j,,], 2, function(x) all(x>=L[j, ]))
    bU = apply(samples[j,,], 2, function(x) all(x<=U[j, ]))
    between = aL*bU
  })
  between = apply(between, 1, min)
  return(mean(between))
}

pointwise_pvalues = function(samples){
  B = dim(samples)[3]+1
  pos = apply(samples, 1:2, function(x) sum(x>0))+1
  neg = apply(samples, 1:2, function(x) sum(x<0))+1
  ppos = 1-pos/B
  pneg = 1-neg/B
  pval_npm = pmin(2 * pmin(ppos, pneg), 1)
  return(pval_npm)
}

transform_obj = function(obj, mat){
  if(!("LSVCMM.Boot" %in% class(obj))) stop("obj must be of class LSVCMM.Boot, i.e., the outpute of a lsvcmm.boot call")
  new_obj=list(
    nvc_boot=obj$nvc_boot,
    nvc=obj$nvc,
    unscaled_time=obj$unscaled_time
  )
  vc_boot = array(0, dim=c(dim(mat)[1], dim(obj$vc_boot)[2], dim(obj$vc_boot)[3]))
  for(b in seq(dim(obj$vc_boot)[3])) vc_boot[,,b] = mat %*% obj$vc_boot[,,b]
  vc = mat %*% obj$vc
  new_obj$vc = vc
  new_obj$vc_boot = vc_boot
  class(new_obj) = c("LSVCMM.Boot")
  return(new_obj)
}
