confidence_band = function(
    obj,
    level=0.95,
    method=c("quantile"),
    var=2
){
  if(!("LSVCMM.Boot" %in% class(obj))) stop("obj must be of class LSVCMM.Boot, i.e., the outpute of a lsvcmm.boot call")
  stopifnot(0<level & level<1)
  level = 1-level
  method = match.arg(method)
  samples = obj$vc_boot[var,,,drop=F]  # 1 x nt x B

  out = switch(
    method,
    "quantile"=quantile_confidence_band(samples, level)
  )

  pval = pointwise_pvalues(samples[1,,])
  out = data.frame(
    L=out$L,
    U=out$U,
    pointwise_confidence=1-2*out$p,
    estimate=obj$vc[var, ],
    mean=apply(samples, 2, mean),
    median=apply(samples, 2, stats::median),
    prop0=apply(samples, 2, function(x) mean(x==0)),
    estimated_time=obj$unscaled_time,
    level=proportion_samples_in_band(out$L, out$U, samples),
    excludes_zero=(out$L>0) | (out$U<0),
    pval=pval,
    pval.adj=stats::p.adjust(pval, method="BH")
  )

  return(out)
}


#' @importFrom pracma logseq
quantile_confidence_band = function(samples, level){
  nc = dim(samples)[2] * dim(samples)[1]
  minp = level / (2*nc)
  maxp = level / 2
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
  L = qs[1,,] # NB this may be a vector or a matrix, cannot use drop=F ...
  U = qs[2,,]
  return(list(L=L, U=U, p=p))
}

proportion_samples_in_band = function(L, U, samples){
  aL = apply(samples, 3, function(x) all(x>=L))
  bU = apply(samples, 3, function(x) all(x<=U))
  between = aL*bU
  return(mean(between))
}

pointwise_pvalues = function(samples){
  B = dim(samples)[2]+1
  pos = apply(samples, 1, function(x) sum(x>0))+1
  neg = apply(samples, 1, function(x) sum(x<0))+1
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
