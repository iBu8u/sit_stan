loglik_rpdf <- function(mcmc,m,fb,trialwise=TRUE) {
  

  mcmc <- as.matrix(mcmc)
  names <- colnames(mcmc)
  sz <- dim(mcmc)
  
  nTrials   <- 64
  nSamples  <- sz[1]
  nSubjects <- 23

  if (trialwise == TRUE) log_lik <- array(0,dim=c(nSamples,nTrials,nSubjects))
  if (trialwise ==FALSE) log_lik <- array(0,dim=c(nSamples,nSubjects))
  
  for (sub in 1:nSubjects) {
    cat("Subject",sub,"/",nSubjects,"\n")
    r <- mcmc[,which(names==paste0("r[",sub,"]"))]
    p <- mcmc[,which(names==paste0("p[",sub,"]"))]
    d <- mcmc[,which(names==paste0("d[",sub,"]"))]   
    f <- mcmc[,which(names==paste0("f[",sub,"]"))]
    
    for (i in 1:nSamples) {
      #cat("r =",r[i],"p =",p[i],"d =",d[i],"f =",f[i])
      a <- matrix(0,3,nTrials+1)
      s <- matrix(0,3,nTrials)
      prob <- matrix(0,1,nTrials)
  
      a[,1] <- c(1/3,1/3,1/3)
      
      lik <- matrix(0,1,nTrials)
      
      for (t in 1:nTrials) {
        
        prob[t] <- ( m[,t,sub] %*% a[,t]^d[i] ) / sum( a[,t]^d[i] + 1.0e-10)
        if ( prob[t] <= 0 ) prob[t] <- 1.0e-10
        if ( prob[t] > 1 ) prob[t] <- 1-1.0e-10
        
        if (fb[sub,t] == 1) {
          s[,t] <- (m[,t,sub] * a[,t]^f[i]) / sum( m[,t,sub] * a[,t]^f[i] + 1.0e-10)
          a[,t+1] <- (1-r[i])*a[,t] + r[i]*s[,t]
        }
        if (fb[sub,t] == 0) {
          s[,t] <- ((1-m[,t,sub]) * a[,t]^f[i]) / sum( (1-m[,t,sub]) * a[,t]^f[i] + 1.0e-10)
          a[,t+1] <- (1-p[i])*a[,t] + p[i]*s[,t]
        }
          
        lik[t] <- log(prob[t])
        
      }

      if (trialwise == TRUE) log_lik[i,,sub] <- lik
      if (trialwise ==FALSE) log_lik[i,sub] <- sum(lik)
      
    }
        
  }
  
  return(log_lik)
}