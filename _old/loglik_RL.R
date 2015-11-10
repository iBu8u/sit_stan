loglik_RL <- function(mcmc, choice, reward, trialwise=FALSE) {
  ## compute log likelihood using mcmc returned by JAGS or Stan
  #  input:
  #   mcmc - combined mcmc object 
  #   trialwise - for RL model, false, for tasks whose single trial is independent, TRUE
  
  mcmc   <- as.matrix(mcmc)
  pnames <- colnames(mcmc)  # parameter names
  sz  <- dim(mcmc)
  eps <- 1.0e-10
  
  nSubjects <- dim(choice)[1]
  nTrials   <- dim(choice)[2]
  nSamples  <- sz[1]
  
  if (trialwise == TRUE) log_lik <- array(0,dim=c(nSamples,nTrials,nSubjects))
  if (trialwise ==FALSE) log_lik <- array(0,dim=c(nSamples,nSubjects))
  
  
  cat("Calculating log_lik for simple RL model... \n")
  startTime = Sys.time()
  print(startTime)
  
  ## --- start calculating
  for (subj in 1:nSubjects) {
    cat("Subject",subj,"/",nSubjects,"\n")
    
    # --- extract samples for per parameter per subject
    lr   <- mcmc[,which(pnames==paste0("lr[",subj,"]"))]
    tau  <- mcmc[,which(pnames==paste0("tau[",subj,"]"))]
    
    for (i in 1:nSamples) {
      #cat("lr =",lr[i],"tau =",tau[i])
      v    <- matrix(0,nTrials+1,2)
      pe   <- matrix(0,nTrials,1)
      prob <- matrix(0,nTrials,2)
      
      lik <- matrix(0,nTrials,1)
      
      for (t in 1:nTrials) {
        
        c <- choice[subj,t]
        r <- reward[subj,t]
        
        prob[t,1] <- 1 / (1 + exp(tau[i] * (v[t,2]-v[t,1]) ))
        prob[t,2] <- 1 / (1 + exp(tau[i] * (v[t,1]-v[t,2]) ))
        
        pe[t,1]   = r - v[t,c];
        
        v[t+1,]  <- v[t,]
        v[t+1,c] <- v[t,c] + lr[i] * pe[t,1]
        
        lik[t,1] <- log(prob[t,c])
      } # trial loop
      
      if (trialwise == TRUE) log_lik[i,,subj] <- lik
      if (trialwise ==FALSE) log_lik[i,subj]  <- sum(lik)
      
    }   # sample loop
  }     # subject loop
  
  cat("Finishing log_lik calculation ... \n")
  endTime = Sys.time()
  print(endTime)  
  cat("It took",as.character.Date(endTime - startTime), "\n")
  
  return(log_lik)
}