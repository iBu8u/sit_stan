run_model <- function(modelStr, fitOBJ=NA) {

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

  } else if (modelStr == "RevLearn_RLbeta_alt1_c" || modelStr == "RevLearn_RLbeta_alt1_bc" || 
             modelStr == "RevLearn_RLbeta_alt2_c" || modelStr == "RevLearn_RLbeta_alt2_bc") { 
    
    chswtch <- array(0,dim = c(ns,nt))
    bet1    <- array(0,dim = c(ns,nt)); bet2    <- array(0,dim = c(ns,nt))
    with    <- array(0,dim = c(ns,nt)); against <- array(0,dim = c(ns,nt))
    my1     <- 0; other1  <- c(0,0,0,0)
    otherChoice <- array(0,dim = c(ns,nt,4)); otherReward <- array(0,dim = c(ns,nt,4))
    pref        <- array(0,dim = c(ns,nt,4)); wOthers     <- array(0,dim = c(ns,nt,4)) # others' weight [.75 .5 .25 .25]
    wghtValue   <- array(0,dim = c(ns,nt,2)) # others' value based on weight
    cfsC2       <- array(0,dim = c(ns,nt,4)) # cumulative-window frequency, same as my C2 
    cfoC2       <- array(0,dim = c(ns,nt,4)) # cumulative-window frequency, opposite to my C2
    
    chswtch <- t(mydata[,5,])
    bet1    <- t(mydata[,13,]); bet2    <- t(mydata[,19,])
    for (s in 1:ns) {
      otherChoice[s,,] <- mydata[,6:9,s]
      otherReward[s,,] <- mydata[,24:27,s]
      pref[s,,]        <- mydata[,47:50,s]
      wOthers[s,,]     <- mydata[,51:54,s]
      wghtValue[s,,]   <- mydata[,59:60,s]
      cfsC2[s,,]       <- mydata[,61:64,s]
      cfoC2[s,,]       <- mydata[,65:68,s]
      
      for (t in 1:nt){
        my1 <- mydata[t,3,s]; other1 <- mydata[t,6:9,s]
        with[s,t]    <- length(which(other1==my1))  # count of with, either 1, 2, 3, or 4
        against[s,t] <- length(which(other1!=my1))  # count of against, either 1, 2, 3, or 4
      }
    }
    
    dataList$chswtch <- chswtch;
    dataList$otherChoice <- otherChoice
    dataList$otherReward <- otherReward
    dataList$wghtValue   <- wghtValue
    dataList$bet1    <- bet1;   dataList$bet2    <- bet2
    dataList$with    <- with;   dataList$against <- against
    dataList$pref    <- pref;   dataList$wOthers <- wOthers
    dataList$cfsC2   <- cfsC2;  dataList$cfoC2   <- cfoC2;
    
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
  nSamples <- 10#4000
  nChains  <- 1#4 
  nBurnin  <- 0#floor(nSamples/2)
  nThin    <- 1
  
  # parameter of interest (this could save both memory and space)
  poi <- create_pois(modelStr)

  #### run stan ####  ==============================================================================
  cat("Estimating", modelStr, "model... \n")
  startTime = Sys.time(); print(startTime)
  
  cat("Calling", nChains, "simulations in Stan... \n")
  stanfit <- stan(modelFile,
                fit     = fitOBJ,
                data    = dataList,
                pars    = poi,
                chains  = nChains,
                iter    = nSamples,
                warmup  = nBurnin,
                thin    = nThin,
                init    = "random",
                verbose = FALSE)
  
  cat("Finishing", modelStr, "model simulation ... \n")
  endTime = Sys.time(); print(endTime)  
  cat("It took",as.character.Date(endTime - startTime), "\n")
  
  return(stanfit)
}  # function run_model()


#### nested functions #### ===========================================================================
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
  } else if (model == "RevLearn_RLbeta_alt1_bc") {
    pois <- c("lr_mu", "thrs_mu", "beta_mu",
              "lr_sd", "thrs_sd", "beta_sd",
              "lr", "thrs", "beta",
              "log_likc1", "log_likc2", "log_likb1", "log_likb2", "lp__")
  } else if (model == "RevLearn_RLbeta_alt1_c") {
    pois <- c("lr_mu", "beta_mu",
              "lr_sd", "beta_sd",
              "lr", "beta",
              "log_likc1", "log_likc2", "lp__")
  } else if (model == "RevLearn_RLbeta_alt2") {
    pois <- c("lr_mu", "thrs_mu", "evid_wght_mu", "beta_mu",
              "lr_sd", "thrs_sd", "evid_wght_sd", "beta_sd",
              "lr", "thrs", "evid_wght", "beta",
              "log_likc1", "log_likc2", "log_likb1", "log_likb2", "lp__")
  }
  
  return(pois)
} # function

#### end of function ####