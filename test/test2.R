# Simulate data
set.seed(123)
n <- 100
p <- 5
Sigma <- diag(p)
Sigma[1,2] <- Sigma[2,1] <- 0.5
Sigma[2,3] <- Sigma[3,2] <- -0.3
X <- MASS::mvrnorm(n, rep(0, p), Sigma)
G <- cor(X) != 0

# Fit model using Gmcmc function
result <- Gmcmc(G, n.iter = 100, n.burnin = 10)

# Extract posterior samples
cloc_samples <- result$cloc
alpha_samples <- result$alpha

# Compute mean and standard deviation of posterior samples
cloc_mean <- apply(cloc_samples, c(1,2), mean)
cloc_sd <- apply(cloc_samples, c(1,2), sd)
alpha_mean <- apply(alpha_samples, 1, mean)
alpha_sd <- apply(alpha_samples, 1, sd)

# Print results
cat("Posterior means and standard deviations of latent condition locations:\n")
print(data.frame(mean = cloc_mean, sd = cloc_sd))
cat("Posterior means and standard deviations of condition-specific intercepts:\n")
print(data.frame(mean = alpha_mean, sd = alpha_sd))
