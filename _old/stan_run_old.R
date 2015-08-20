run_model <- function(modelStr, runPara=FALSE, fitOBJ=NA) {

  library(rstan); library(parallel); library(loo)
  rstan_options(auto_write = TRUE)
  options(mc.cores = 4)
  
  #### prepare data #### ===========================================================================
  load("_data/sit_reversal_betnoshow_129.rdata")
  dataList <- list()
  sz <- dim(mydata)
  nt <- sz[1]; ns <- sz[3]
  dataList$nSubjects <- ns; dataList$nTrials <- nt
  
  choice1 <- array(0,dim = c(ns,nt)); choice2 <- array(0,dim = c(ns,nt)); reward <- array(0,dim = c(ns,nt))
  choice1 <- t(mydata[,3,])   # 1 OR  2
  choice2 <- t(mydata[,10,])  # 1 OR  2
  reward  <- t(mydata[,14,])  # 1 OR -1
  dataList$choice1 <- choice1
  dataList$choice2 <- choice2
  dataList$reward  <- reward
  
  if (modelStr == "RevLearn_RLcoh" || modelStr == "RevLearn_RLcoh_cfa" || 
             modelStr == "RevLearn_RLcoh_2lr" || modelStr == "RevLearn_RLcoh_2lr_cfa" ||
             modelStr == "RevLearn_RLcoh_modvalue" || modelStr == "RevLearn_RLcoh_modprob_tempin" || 
             modelStr == "RevLearn_RLcoh_modprob_tempout" || modelStr == "RevLearn_RLcoh_2lr2t_modvalue" ||
             modelStr == "RevLearn_RLcoh_2lr2t_modprob_tempin" || modelStr == "RevLearn_RLcoh_2lr2t_modprob_tempout" ||
             modelStr == "RevLearn_RLcoh_2lr2t_modprob_tempin_bern" ||
             modelStr == "RevLearn_RLcoh_2lr_bet") {
    
    chswtch <- array(0,dim = c(ns,nt))
    bet1    <- array(0,dim = c(ns,nt)); bet2    <- array(0,dim = c(ns,nt))
    with    <- array(0,dim = c(ns,nt)); against <- array(0,dim = c(ns,nt))
    my1     <- 0; other1  <- c(0,0,0,0)
    
    chswtch <- t(mydata[,5,])
    bet1    <- t(mydata[,13,]); bet2    <- t(mydata[,19,])

    for (s in 1:ns) {
      for (tr in 1:nt){
        my1 <- mydata[tr,3,s]; other1 <- mydata[tr,6:9,s]
        with[s,tr]    <- length(which(other1==my1)) /4  # count of with, either 1, 2, 3, or 4, divided by 4
        against[s,tr] <- length(which(other1!=my1)) /4  # count of against, either 1, 2, 3, or 4, divided by 4 
      }
    }
    
    dataList$chswtch <- chswtch
    dataList$bet1    <- bet1;    dataList$bet2    <- bet2
    dataList$with    <- with;    dataList$against <- against

  } else if (modelStr == "RevLearn_RLbeta_alt1" || modelStr == "RevLearn_RLbeta_alt2") { 
    
    chswtch <- array(0,dim = c(ns,nt))
    bet1    <- array(0,dim = c(ns,nt)); bet2    <- array(0,dim = c(ns,nt))
    with    <- array(0,dim = c(ns,nt)); against <- array(0,dim = c(ns,nt))
    my1     <- 0; other1  <- c(0,0,0,0)
    otherChoice <- array(0,dim = c(ns,nt,4));
    otherind    <- array(0,dim = c(ns,nt,4));
    
    chswtch <- t(mydata[,5,])
    bet1    <- t(mydata[,13,]); bet2    <- t(mydata[,19,])

    for (s in 1:ns) {
      otherChoice[s,,] <- mydata[,6:9,s]
    }

    for (s in 1:ns) {
      for (tr in 1:nt){
        my1 <- mydata[tr,3,s]; other1 <- mydata[tr,6:9,s]
        with[s,tr]    <- length(which(other1==my1)) /4  # count of with, either 1, 2, 3, or 4, divided by 4
        against[s,tr] <- length(which(other1!=my1)) /4  # count of against, either 1, 2, 3, or 4, divided by 4 
      }
    }
    
    dataList$chswtch <- chswtch
    dataList$bet1    <- bet1;    dataList$bet2    <- bet2
    dataList$with    <- with;    dataList$against <- against

    
    
  } else if (modelStr == "RevLearn_RLcumrew" || modelStr == "RevLearn_RLcumrew_cfa" || 
             modelStr == "RevLearn_RLcumrew_2lr" || modelStr == "RevLearn_RLcumrew_2lr_cfa" ) {
    
    otherChoice <- array(0,dim = c(ns,nt,4));  otherReward <- array(0,dim = c(ns,nt,4))

    for (s in 1:ns) {
      otherChoice[s,,] <- mydata[,6:9,s]
      otherReward[s,,] <- mydata[,24:27,s]
    }
    
    dataList$otherChoice <- otherChoice; dataList$otherReward <- otherReward
  }
    
  #### preparation for running stan #### ============================================================
  # model string in a separate .stan file
  modelFile <- paste0("_scripts/",modelStr,".stan")

  # setup up Stan configuration
  nSamples <- 16000 #12000
  nChains  <- 4 
  nBurnin  <- floor(nSamples/2)
  nThin    <- 8 #3
  
  # initial parameters, OR use inits from Stan
  inits <- create_inits(modelStr,nChains,ns,runPara)
  # parameter of interest (this could save both memory and space)
  poi <- create_pois(modelStr)

  #### run stan ####  ==============================================================================
  cat("Estimating", modelStr, "model... \n")
  startTime = Sys.time(); print(startTime)
  
  if ( runPara == TRUE ) {
    # this is for running stan on multiple processor cores (taken from Andy Gelman's website/blog)
    # [update!] from Stan 2.7.0 on, Stan automatically supports sampling in parallel
    cat("Calling", nChains, "simulations using parallel method... \n")
    
    foo  <- stan(modelFile, fit = fitOBJ, data = dataList, iter=1, chains=1)
    para_results <- mclapply(1:nChains, mc.cores=nChains, FUN=function(chain) {
      stan(fit     = foo,
           data    = dataList,
           pars    = poi,
           chains  = 1,
           iter    = nSamples,
           warmup  = nBurnin,
           thin    = nThin,
           # init    = inits,
           init    = "random",
           verbose = TRUE,
           refresh = -1,
           chain_id= chain)})
    stanfit <- sflist2stanfit(para_results)  # this takes only within 1 sec
    
  } else { 
    cat("Calling", nChains, "simulations using serial method... \n")
    
    stanfit <- stan(modelFile,
                    fit     = fitOBJ,
                    data    = dataList,
                    pars    = poi,
                    chains  = nChains,
                    iter    = nSamples,
                    warmup  = nBurnin,
                    thin    = nThin,
                    #init    = inits,
                    init    = "random",
                    verbose = FALSE)
  }
  
  cat("Finishing", modelStr, "model simulation ... \n")
  endTime = Sys.time(); print(endTime)  
  cat("It took",as.character.Date(endTime - startTime), "\n")
  
  return(stanfit)
}  # function run_model()


#### nested functions #### ===========================================================================
create_inits <- function(model, n.chains, n.subjects, paraType) {
  
  inits <-  list()
  if (paraType == TRUE)   n.chains = 1
  
  for (c in 1:n.chains) {
    
    lr_mu_pr  = qnorm(runif(1,0.45,0.55))
    tau_mu_pr = qnorm(runif(1,1.5,2.5)/10)
    lr_sd     = runif(1,0,1.5)
    tau_sd    = runif(1,0,1.5)
    lr_row    = qnorm(runif(n.subjects,0.4,0.6))
    tau_row   = qnorm(runif(n.subjects,1.5,2.5)/10)
    
    if ( model == "RevLearn_RL") {
      
      inits[[c]] <- list(lr_mu_pr  = lr_mu_pr,
                         tau_mu_pr = tau_mu_pr,
                         lr_sd     = lr_sd,
                         tau_sd    = tau_sd,
                         lr_row    = lr_row,
                         tau_row   = tau_row)
    } else if (model == "RevLearn_RLcoh"){
      coha_mu  = runif(1,0.02,0.1)
      cohw_mu  = runif(1,0.02,0.1)
      coha_sd  = runif(1,0.02,0.5)
      cohw_sd  = runif(1,0.02,0.5)
      coha     = runif(n.subjects,0.02,0.1)
      cohw     = runif(n.subjects,0.02,0.1)

      inits[[c]] <- list(lr_mu_pr=lr_mu_pr, tau_mu_pr=tau_mu_pr, coha_mu =coha_mu, cohw_mu =cohw_mu,
                         lr_sd   =lr_sd,    tau_sd   =tau_sd,    coha_sd =coha_sd, cohw_sd =cohw_sd,
                         lr_row  =lr_row,   tau_row  =tau_row,   coha    =coha,    cohw    =cohw   )
    
    }
  } # for

    return(inits)
} # function


create_pois <- function(model){
  pois <- list()
  
  if (model == "RevLearn_RL" || model == "RevLearn_RLnc"){
    pois <- c("lr_mu", "tau_mu", 
              "lr_sd", "tau_sd",
              "lr", "tau", 
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLnc_2lr") {
    pois <- c("lr1_mu", "lr2_mu", "tau_mu",  
              "lr1_sd", "lr2_sd", "tau_sd", 
              "lr1", "lr2", "tau",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLnc_cfa") {
    pois <- c("lr_mu", "tau_mu",  "cfa_mu", 
              "lr_sd", "tau_sd", "cfa_sd",
              "lr", "tau", "cfa",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLnc_2lr_cfa") {
    pois <- c("lr1_mu", "lr2_mu", "tau_mu", "cfa_mu",
              "lr1_sd", "lr2_sd", "tau_sd", "cfa_sd",
              "lr1", "lr2", "tau", "cfa",
              "c_rep",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcoh"){
    pois <- c("lr_mu", "tau_mu", "coha_mu", "cohw_mu", 
              "lr_sd", "tau_sd", "coha_sd", "cohw_sd",
              "lr", "tau", "coha", "cohw",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcoh_modvalue" || model == "RevLearn_RLcoh_modprob_tempin" || model == "RevLearn_RLcoh_modprob_tempout"){
    pois <- c("lr_mu", "tau_mu", "coha_mu", "cohw_mu", 
              "lr_sd", "tau_sd", "coha_sd", "cohw_sd",
              "lr", "tau", "coha", "cohw",
              "log_lik1", 
              "log_lik2", "lp__")
  } else if (model == "RevLearn_RLcoh_2lr2t_modvalue" || model == "RevLearn_RLcoh_2lr2t_modprob_tempin" || 
             model == "RevLearn_RLcoh_2lr2t_modprob_tempout" || model == "RevLearn_RLcoh_2lr2t_modprob_tempin_bern") {   
    pois <- c("lr1_mu", "tau1_mu", "lr2_mu", "tau2_mu", "coha_mu", "cohw_mu", 
              "lr1_sd", "tau1_sd", "lr2_sd", "tau2_sd", "coha_sd", "cohw_sd",
              "lr1", "tau1", "lr2", "tau2", "coha", "cohw",
              "log_lik1", 
              "log_lik2", "lp__")
  } else if (model == "RevLearn_RLcoh_2lr"){
    pois <- c("lr1_mu", "lr2_mu", "tau_mu", "coha_mu", "cohw_mu", 
              "lr1_sd", "lr2_sd", "tau_sd", "coha_sd", "cohw_sd",
              "lr1", "lr2", "tau", "coha", "cohw",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcoh_cfa") {
    pois <- c("lr_mu", "tau_mu", "coha_mu", "cohw_mu", "cfa_mu",
              "lr_sd", "tau_sd", "coha_sd", "cohw_sd", "cfa_sd",
              "lr", "tau", "coha", "cohw", "cfa",
              "log_lik", "lp__")
  } else if (model == "RevLearn_RLcoh_2lr_cfa") {
    pois <- c("lr1_mu", "lr2_mu", "tau_mu", "coha_mu", "cohw_mu", "cfa_mu",
              "lr1_sd", "lr2_sd", "tau_sd", "coha_sd", "cohw_sd", "cfa_sd",
              "lr1", "lr2", "tau", "coha", "cohw", "cfa", 
              "log_lik", "lp__")
  }
  
  return(pois)
} # function


#### end of function ####