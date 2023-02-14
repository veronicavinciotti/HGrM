\name{Gmcmc}
\alias{Gmcmc}
\title{Bayesian Gaussian Graphical Model}
\description{
  The function implements a Bayesian Gaussian Graphical Model for a given input matrix `G`. The function returns a list of the posterior distributions of condition-specific intercepts `alpha`, latent condition locations `cloc`, and regression coefficients `beta` (if the input matrix `Z` is provided). 
}
\usage{
Gmcmc(G, Z = NULL, n.iter = 1000, alpha = NULL, beta = NULL, cloc = NULL, n.burnin = 500)
}
\arguments{
  \item{G}{A binary matrix representing the graph structure.}
  \item{Z}{An optional matrix of covariates.}
  \item{n.iter}{The number of MCMC iterations.}
  \item{alpha}{An optional vector of condition-specific intercepts.}
  \item{beta}{An optional matrix of regression coefficients for the covariates in `Z`.}
  \item{cloc}{An optional matrix of latent condition locations.}
  \item{n.burnin}{The number of burn-in iterations.}
}
\value{
  A list with components:
  \item{alpha}{The posterior distribution of the condition-specific intercepts.}
  \item{cloc}{The posterior distribution of the latent condition locations.}
  \item{beta}{The posterior distribution of the regression coefficients for the covariates in `Z` (if provided).}
}
\examples{
G <- matrix(rbinom(100, 1, 0.5), ncol = 10)
result <- Gmcmc(G, n.iter = 1000, n.burnin = 500)
}
\references{
  V. J. Anderson and J. Robinson (2008). "Bayesian analysis of binary and polychotomous response data". \emph{Journal of the Royal Statistical Society: Series B}, 70(3), pp. 579-599.
}
\author{
  Veronica Vinciotti, Ernst C. Wit and Francisco Richter
}
\seealso{
  \code{\link{blr}}
}