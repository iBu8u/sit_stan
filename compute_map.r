compute_map <- function(mcmc, meth = "density") {
  
  library(modeest)
  startTime = Sys.time()
  
  sz <- dim(mcmc) 
  map <- array(0,dim=c(sz[2],1))
  
  for (f in 1:sz[2]) {
    tmp <- mlv(mcmc[,f], method=meth)
    map[f] <- mean(tmp$M)
  }
  rownames(map) <- colnames(mcmc)
  
  endTime = Sys.time()
  cat("It took",as.character.Date(endTime - startTime), "\n")
  return(map)
}
