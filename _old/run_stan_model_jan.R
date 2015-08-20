run_model <- function() {
  
  # fit hierarchical RPDF model using Stan
  
  library(rstan)
  # library(parallel)
  
  # setwd('/Users/glaescher/Projects/ppgdata/wcst/stan')

  modelStr <- "RL"
  
  sink("RL.log")
  
  # setup up Stan configuration
  # nCores   <- 4
  nSamples <- 2000
  nChains  <- 4
  # nSamplesPerChain <- ceiling(nSamples/nChains)
  nBurnin <- 2000
  nThin    <- 5
  
  # model file
  modelFile <- paste0(modelStr,".stan")

  # data
  source('data_small.dump') # or load('...')
  dataList <- list(choice = choice,
                   reward = reward,
                   nTrials = nTrials,
                   nSubjects = nSubjects)
  
  # create inits
  inits <- create_inits(modelStr,nChains,nSubjects)
  
  # Estimate Model
  cat("Estimating RL model... \n")
  flush.console()

  # this is for running stan on multiple processor cores (taken from Andy Gelman's website / blog)
  #para_result <- mclapply(1:nChains,mc.cores=nCores,FUN=function(chain) 
  #  {stan(modelFile,data=dataList,chains=1,thin=nThin,warmup=nBurnin,iter=nPerChain,verbose=FALSE)})
  #fit <- sflist2stanfit(para_results)

  fit <- stan(modelFile,data=dataList) # defaults: thin=1, iter=2000, warmup=(iter/2), chains=4, init="random"
  # fit <- stan(modelFile,data=dataList,chains=nChains,thin=nThin,iter=nSamples,warmup=nBurnin,init=inits,verbose=FALSE)
  
  
  return(fit)
} 

create_inits <- function(model,n.chains,n.subjects) {
  
  library(random)
  inits <- list()
  prior.a <- 5
  prior.b <- 5
  
  for (c in 1:n.chains)
    
    if ( model == "RL") {
      
      lr.mu <- rbeta(1,prior.a,prior.b)
      lr.kappa0 <- rbeta(1,prior.a,prior.b)
      
      temp.mu <- rbeta(1,prior.a,prior.b)
      temp.kappa0 <- rbeta(1,prior.a,prior.b)
      
      lr <- rbeta(n.subjects,prior.a,prior.b)
      temp0 <- rbeta(n.subjects,prior.a,prior.b)
      
      inits[[c]] <- list(lr_mu=lr.mu,
                         temp_mu=temp.mu,
                         lr_kappa0=lr.kappa0,
                         temp_kappa0=temp.kappa0,
                         lr=lr,temp0=temp0)
    }
  return(inits);  
}