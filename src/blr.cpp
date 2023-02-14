#include <RcppArmadillo.h>
#include <math.h>
#include <random>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// Define a class to store the parameters of the BLR model
class BLR_Model {
public:
  int n;
  int D;
  int N1;
  int N0;
  mat X;
  vec y;
  vec offset;
  vec theta_0;
  mat Q_0;
};

// Define a class to store the posterior distribution of the parameters
class BLR_Posterior {
public:
  mat theta_chain;
};

// Function to draw from a truncated normal distribution
double rtruncnorm(double mean, double sd, double a, double b,
                  mt19937& rng) {
  normal_distribution<double> norm_dist(mean, sd);
  double x;
  do {
    x = norm_dist(rng);
  } while (x < a || x > b);
  return x;
}

// Function to perform Gibbs sampling
BLR_Posterior blr_gibbs_sampler(const BLR_Model& model, int N_sim,
                                unsigned long seed) {
  // Initialize the random number generator
  mt19937 rng(seed);
  
  // Initialize the latent variable z
  vec z(model.n);
  z.fill(NA_REAL);
  
  // Initialize the posterior distribution of the parameters
  BLR_Posterior posterior;
  posterior.theta_chain = mat(N_sim, model.D);
  
  // Compute the posterior variance of theta
  mat prec_0 = inv(model.Q_0);
  mat V = inv(prec_0 + model.X.t() * model.X);
  
  for (int t = 0; t < N_sim; t++) {
    // Update the mean of z
    vec mu_z = model.X * model.theta_0 + model.offset;
    
    // Draw the latent variable z from its full conditional
    for (int i = 0; i < model.n; i++) {
      if (model.y[i] == 0) {
        z[i] = rtruncnorm(mu_z[i], 1.0, -INFINITY, 0.0, rng);
      } else if (model.y[i] == 1) {
        z[i] = rtruncnorm(mu_z[i], 1.0, 0.0, INFINITY, rng);
      }
    }
    
    // Compute the posterior mean of theta
    vec M = V * (prec_0 * model.theta_0 + model.X.t() * (z - model.offset));
    
    // Draw theta from its full conditional
    normal_distribution<double> norm_dist;
    vec theta(model.D);
    for (int j = 0; j < model.D; j++) {
      norm_dist.param(normal_distribution<double>::param_type(M[j], sqrt(V(j, j))));
      theta[j] = norm_dist(rng);
    }
    posterior.theta_chain.row(t) = theta.t();
  }
  
  return posterior;
}
    

      
                   