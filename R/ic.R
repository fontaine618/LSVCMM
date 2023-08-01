ebic_single = function(object, extended=0., kernel=T){
  llk = object$results$llk
  df = ifelse(kernel, object$results$df_kernel, object$results$df)
  df_logn = ifelse(kernel, object$results$df_logn_kernel, object$results$df_logn)
  maxdf = prod(dim(object$B))
  return(
    -2*llk + df_logn + extended * df * log(maxdf)
  )
}
