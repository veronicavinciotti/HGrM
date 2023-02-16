G <- matrix(rbinom(100, 1, 0.5), ncol = 10)
result <- Gmcmc(G, n.iter = 100, n.burnin = 500)
