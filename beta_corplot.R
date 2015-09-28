beta_corplot <- function(stanfit, splt = NULL) {
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
  
  if (is.null(splt) ) {
    corrMat <- cor(bMat)
    corrplot(corrMat, method = 'square', diag = FALSE, addCoef.col = T)
  } else if (splt == 6) {
    bMat_pos <- bMat[bMat[,6]>=0,]
    bMat_neg <- bMat[bMat[,6] <0,]
    corrMat_pos <- cor(bMat_pos)
    corrMat_neg <- cor(bMat_neg)
    par(mfrow = c(1,2))
    corrplot(corrMat_pos, method = 'square', diag = FALSE, addCoef.col = T)
    corrplot(corrMat_neg, method = 'square', diag = FALSE, addCoef.col = T)
  } else if (splt == 3) {
    bMat_pos <- bMat[bMat[,3]>=0,]
    bMat_neg <- bMat[bMat[,3] <0,]
    corrMat_pos <- cor(bMat_pos)
    corrMat_neg <- cor(bMat_neg)
    par(mfrow = c(1,2))
    corrplot(corrMat_pos, method = 'square', diag = FALSE, addCoef.col = T)
    corrplot(corrMat_neg, method = 'square', diag = FALSE, addCoef.col = T)
  }
}