#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix blr(NumericVector y, NumericMatrix X, double offset = 0, 
                  NumericVector theta, double theta_0 = 0, int N_sim = 1) {
  int n = y.length();
  int N1 = sum(y);
  int N0 = n - N1;
  int D = X.ncol();
  
  NumericMatrix Q_0(D, D);
  for (int i = 0; i < D; i++) {
    Q_0(i, i) = 10;
  }
  
  NumericVector z(n);
  NumericMatrix theta_chain(N_sim, D);
  
  NumericMatrix prec_0 = solve(Q_0);
  NumericMatrix V = solve(prec_0 + tcrossprod(X));
  
  NumericVector mu_z(n);
  NumericVector M(D);
  NumericVector temp(D);
  
  for (int t = 0; t < N_sim; t++) {
    mu_z = X * theta + offset;
    if (N0 > 0) {
      NumericVector z0 = rtruncnorm(N0, mu_z[y == 0], 1, -INFINITY, 0);
      int j = 0;
      for (int i = 0; i < n; i++) {
        if (y[i] == 0) {
          z[i] = z0[j];
          j++;
        }
      }
    }
    if (N1 > 0) {
      NumericVector z1 = rtruncnorm(N1, mu_z[y == 1], 1, 0, INFINITY);
      int j = 0;
      for (int i = 0; i < n; i++) {
        if (y[i] == 1) {
          z[i] = z1[j];
          j++;
        }
      }
    }
    
    M = V * (prec_0 * theta_0 + tcrossprod(X, z - offset));
    temp = rnorm(1, M, V);
    for (int i = 0; i < D; i++) {
      theta[i] = temp[i];
    }
    
    for (int i = 0; i < D; i++) {
      theta_chain(t, i) = theta[i];
    }
  }
  
  return theta_chain;
}