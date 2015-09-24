WAIC <- function(input){
  colVars <- function(M){
    vars <- M[1, ]
    for (n in seq_along(vars)) vars[n] <- var(M[, n])
    return(vars)
  }
  
  if ( class(input)=="stanfit") log_lik <- extract (input, "log_lik")$log_lik
  if ( class(input)!="stanfit") log_lik <- input
  
  dim(log_lik) <- if (length(dim(log_lik)) == 1) c(length(log_lik), 1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  
  ## --- WAIC
  lpd    <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  
  ## --- leave one out
  loo_weights_raw <- 1 / exp(log_lik - max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
    matrix(colMeans(loo_weights_raw), nrow=S, ncol=n, byrow=TRUE)
  loo_weights_regularized <- pmin(loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik) * loo_weights_regularized) / 
                    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  
  ## --- return
  pointwise <- cbind(waic, lpd, p_waic, elpd_waic, p_loo, elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n * colVars(pointwise))
  
  return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
              p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
              pointwise=pointwise, total=total, se=se))
}