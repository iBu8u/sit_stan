N <- 50
x <- rnorm(N)
y <- rnorm(N,3,1.2)
L <- list()
L$N <- N
L$x <- x
L$y <- y

library(rstan)
stanfit <- stan("_scripts/toy_model.stan",data= L, chains = 1,iter = 5, init = "random",verbose = FALSE)



