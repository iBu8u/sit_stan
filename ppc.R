ppc <- function(stanfit, choice = NULL, sid = NULL, gid = NULL) {
  
  ## ppc: posterior predictive check=====================================
  
  #### data preparation and initialize -----------------------------------
  library(rstan); library(ggplot2); library(reshape2)
  load("_data/sit_reversal_betnoshow_129.rdata")
  sz <- dim(data)
  nt <- sz[1]; ns <- sz[3]
  choice1  <- array(0,dim = c(ns,nt)); choice2  <- array(0,dim = c(ns,nt)); reversal <- array(0,dim = c(ns,nt))
  choice1  <-  t(data[,3,]);           choice2  <-  t(data[,10,]);          reversal <- t(data[,2,])
  acc      <- array(0,dim = c(ns,nt))
  
  if ( is.null(choice) )  {
    c_rep  <- extract(stanfit, pars="c_rep", permuted=TRUE)$c_rep
    choice <- choice2 
  } else if (choice == 1) {
    c_rep  <- extract(stanfit, pars="c_rep1", permuted=TRUE)$c_rep1
    choice <- choice1
  } else if (choice == 2) {
    c_rep  <- extract(stanfit, pars="c_rep2", permuted=TRUE)$c_rep2
    choice <- choice2
  }
  
  niter <- dim(c_rep)[1]
  
  #### calculate ppc accuracy --------------------------------------------
  for (s in 1:ns)   {
    for (t in 1:nt) {
      acc[s,t] <- sum(c_rep[,s,t]==choice[s,t]) / niter
    }
  }
  
  acc_sub <- rowMeans(acc)
  acc_grd <- mean(acc)
  
  #### visualize the results =============================================================
  
  ## plot per subject, plus the reversal point -------------------------------------------
  if (!is.null(sid)) {
    p = qplot(1:nt, acc[sid,], geom ='line')
    p = p + geom_vline(xintercept=which(reversal[sid,]==1), colour="red", linetype="longdash")
    print(p)
  }
  
  ## plot for total (grand mean) ----- not too much sence...
  # qplot(1:100, colMeans(acc), geom ='line')
  
  ## plot each subject in a grid (facet), per subject group ------------------------------
  ## note that there are only 129 subjects!!! therefore this plot is not perfect!!!
  if (!is.null(gid)) {
    df1 <- as.data.frame(t(acc))
    colnames(df1) <- c(paste0("subj", 1:ns))
    df1$trial  <-  1:nt
    df2 <- as.data.frame(t(reversal))
    colnames(df2) <- c(paste0("subj", 1:ns))
    df2$trial  <-  1:nt
    
    # long-format data frame
    ldf1 <- melt(df1, id="trial")
    ldf2 <- melt(df2, id="trial")
    ldf  <- ldf1
    colnames(ldf) <- c("trial", "subj","accuracy")
    ldf$reversal  <- ldf2$value
    
    grp <- ldf[(1+(gid-1)*500):(500+(gid-1)*500),]
    
    p2 = ggplot(grp, aes(x = trial, y = accuracy) )
    p2 = p2 + geom_line() 
    p2 = p2 + facet_grid(subj ~ ., scales = "free_y")
    p2 = p2 + geom_vline(xintercept=which(grp$reversal[1:100]==1), colour="red", linetype="longdash")
    print(p2)
  }

  return(list(acc=acc,acc_sub=acc_sub,acc_grd=acc_grd))
}
#### end of fucntion