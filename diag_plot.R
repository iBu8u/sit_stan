source("_scripts/DBDA2E-utilities.R")
library(coda)

mcmcCoda = mcmc.list( lapply( 1:ncol(out1) , function(x) {mcmc(as.array(out1)[,x,])}))

fileNameRoot = "some_file_name"
graphFileType = "eps"

diagMCMC(mcmcCoda, parName=c("lr_mu")
         # , saveName = fileNameRoot
         # , saveType = graphFileType
         )


plotPost(mcmcCoda[,"lr_mu"], main = "lr_mu", xlab = bquote(lr_mu),
         cenTend = "mode", credMass = 0.95, showCurve = FALSE,
         col = "skyblue")
## optional plot input arguments for plotPost:
# border = "skyblue"
# xlim   = c(0,1)
# points( lr , 0 , pch="+" , col="red" , cex=3 )
# points() can be used for ploting the true parameter for simulated data



