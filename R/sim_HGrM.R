sim_HGrM <- function(D = 2, p = 81, B = 10, seed = 123) {
  
  set.seed(seed)
  
  # Simulating latent space model
  cond.loc1 <- rnorm(B, rep(c(0,0.5), each=5), 0.1)
  cond.loc2 <- rnorm(B, rep(c(0.5,0), each=5), 0.1)
  
  # Simulating one covariate
  n.edge <- p*(p-1)/2
  X <- runif(n.edge, -0.5, 0.5)
  X <- as.matrix(X)
  
  # True parameters
  cloc.true <- cbind(cond.loc1, cond.loc2)
  alpha.true <- rnorm(B, -2)
  beta.true <- 2.5
  
  # Simulating true graph
  G.true <- matrix(0, ncol = B, nrow = n.edge)
  for(i in 1:1000) {
    dist.cond <- matrix(ncol=B, nrow=n.edge)
    for (b in 1:B){
      # Updating condition-specific intercept
      dist.cond[,b] <- apply(G.true, 1, function(g, cloc, b) { crossprod(apply(cloc*g, 2, sum) - cloc[b,]*g[b], cloc[b,]) }, cloc = cloc.true, b = b)
    }
    m <- matrix(1:p, ncol=p, nrow=p)
    e1 <- t(m)[lower.tri(m)]
    e2 <- m[lower.tri(m)]
    Pi.true <- matrix(ncol = B, nrow = n.edge)
    for (b in 1:B) {
      for (k in 2:p) {
        for (j in 1:(k-1)) {
          ind <- e1 == j & e2 == k
          Pi.true[ind, b] <- pnorm(alpha.true[b] + dist.cond[ind, b] + X[ind,] %*% beta.true)
        }
      }
    }
    G.true.new <- matrix(rbinom(n.edge * B, 1, Pi.true), ncol = B, nrow = n.edge)
    G.true <- G.true.new
  }
  
 # hist(Pi.true)
  
  # Simulating data (Gaussian)
  n <- 1000
  data <- vector(mode = "list", length = B)
  for (j in 1:B) {
    A <- matrix(0, nrow = p, ncol = p)
    A[lower.tri(A)] <- G.true[,j]
    A <- A + t(A)
    data[[j]] <- BDgraph::bdgraph.sim( p = p, n = n, graph = A)$data
  }
  
  list(data = data, cloc_true = cloc.true, alpha_true = alpha.true, beta_true = beta.true,G.true=G.true)
}
