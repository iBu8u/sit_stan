beta_corplot <- function(stanfit) {
  # beta_corplot() plots the correlation matrix of individuals' BETAs
  
  library(corrplot)
  
  if ( class(stanfit)=='stanfit' ) {
    stanfit <- stanfit
  } else {
    stanfit <- stanfit$fit
  }
  #### get BETAs from Stanfit output -------------------------------------------------------------
  ns <- 129
  parm <- get_posterior_mean(stanfit, 'beta')
  bMat <- matrix(parm[,5], nrow = ns)
  colnames(bMat) <- c('beta1','beta2','beta3','beta4','beta5','beta6')
  corrMat <- cor(bMat)
  
  #### plot ---------------------------------------------------------------------------------------
  corrplot(corrMat, method = 'square', diag = FALSE, addCoef.col = T)
}