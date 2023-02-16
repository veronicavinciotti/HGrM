#devtools::install_github("franciscorichter/HGrM")
library(HGrM)

G <- matrix(rbinom(100, 1, 0.5), ncol = 10)
result <- Gmcmc(G, n.iter = 100, n.burnin = 50)

G.loc<-Gmcmc(G,alpha=sample.alpha[,k],cloc=sample.cloc[,,k],n.iter=1,n.burnin = 0)

G.loc<-Gmcmc(G,n.iter=1,n.burnin = 0)
