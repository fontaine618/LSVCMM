#' @importFrom dplyr bind_rows
confidence_band = function(
    obj,
    level=0.95,
    method=c("quantile"),
    var=2,
    pval.method=c("normal", "percentile")
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
  pval_normal = pointwise_pvalues(samples, "normal")
  pval_percentile = pointwise_pvalues(samples, "percentile")
  nvars = length(var)
  nt = dim(samples)[2]
  jpval = joint_pvalue(samples)
  df = lapply(seq_along(var), function(j){
    data.frame(
      L=out$L[j,],
      U=out$U[j,],
      pointwise_confidence=1-2*out$p,
      estimate=obj$vc[var[j], ],
      mean=apply(samples[j,,], 1, mean),
      median=apply(samples[j,,], 1, stats::median),
      prop0=apply(samples[j,,], 1, function(x) mean(x==0)),
      estimated_time=obj$unscaled_time,
      level=prop_in_band,
      excludes_zero=(out$L[j,]>0) | (out$U[j,]<0),
      pval_normal=pval_normal[j,],
      pval_percentile=pval_percentile[j,],
      pval_joint=joint_pvalue(samples[j,,,drop=F]),
      pval_joint_all=jpval,
      var=var[j]
    )
  }) %>% dplyr::bind_rows()

  return(df)
}

#' @importFrom pracma logseq
quantile_confidence_band = function(samples, level){
  B = dim(samples)[3]
  nc = dim(samples)[2] * dim(samples)[1]
  minp = level / (2*(nc+1))
  maxp = level / 2
  ps = rev(c(seq(minp, maxp, 1/(B+1)), maxp))
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

pointwise_pvalues = function(samples, method=c("normal", "percentile")){
  method = match.arg(method)
  out = switch(
    method,
    "normal"=pointwise_pvalues_normal(samples),
    "percentile"=pointwise_pvalues_percentile(samples)
  )
  return(out)
}


joint_pvalue = function(samples){
  B = dim(samples)[3]
  cs = c(1/(B+1), 0.5, 1.)
  diff = cs[3] - cs[1]
  while(diff > 1/(100*B)){
    # print(cs)
    incl0 = sapply(cs, function(cc){
      out = pointwise_quantile_confidence_band((1-cc)/2, samples)
      !any(out$L > 0 | out$U < 0)
    })
    # print(incl0)
    if(!any(incl0)) break
    firstTrue = which(incl0)[1]
    mid = (cs[firstTrue]+ cs[firstTrue-1])/2
    cs = c(cs[firstTrue-1], mid, cs[firstTrue])
    diff = diff/2
  }
  pval = 1-cs[1]
  pval
}

pointwise_pvalues_normal = function(samples){
  m = apply(samples, 1:2, mean)
  s = apply(samples, 1:2, sd)
  s = pmax(s, 1e-10) # when s==0, we have m=0, so it does not matter
  pvalues = 2*pnorm(-abs(m/s))
  pvalues
}

pointwise_pvalues_percentile = function(samples){
  B = dim(samples)[3]+1
  pos = apply(samples, 1:2, function(x) sum(x>0))
  neg = apply(samples, 1:2, function(x) sum(x<0))
  ppos = 1-pos/B
  pneg = 1-neg/B
  pvalues = pmin(2 * pmin(ppos, pneg), 1)
  pvalues
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
