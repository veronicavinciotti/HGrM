#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List sample_data(List data, List discrete_data, List K, List tpoints){
  int B = data.length();
  List dat(B);
  for (int i = 0; i < B; i++) {
    mat S = as<mat>(K[i]);
    int p = Rf_ncols(data[i]);
    for (int j = 0; j < p; j++) {
      mat S22i = inv(S.submat(0,0,p-2,p-2));
      vec S12 = S.col(p-1);
      double S11 = S(p-1,p-1);
      mat data_mat = as<mat>(data[i]);
      mat data_submat = data_mat.submat(0,0,p-2,p-2);
      vec mu_j = S12.t()*S22i*data_submat.t();
      double var_j = S11 - S12.t()*S22i*S12;
      vec tpoint1 = as<vec>(tpoints[i][0]);
      vec tpoint2 = as<vec>(tpoints[i][1]);
      data_mat.col(j) = rtruncnorm(mu_j.n_elem,tpoint1[j],tpoint2[j],mu_j,sqrt(var_j));
    }
    dat[i] = data_mat;
  }
  return dat;
}


void gmcmc(const MatrixXd &G, MatrixXd &Z, int n_iter, VectorXd &alpha, VectorXd &beta, MatrixXd &cloc, int n_burnin) {
  int B = G.cols
  int n_edge = G.rows();
  int p = (sqrt(1 + 8 * n_edge) + 1) / 2;
  MatrixXd m = MatrixXd::LinSpaced(p, 1, p);
  VectorXd e1 = m.triangularView<Eigen::StrictlyLower>().transpose().colwise().flatten();
  VectorXd e2 = m.triangularView<Eigen::StrictlyLower>().colwise().flatten();
  if (cloc.size() == 0) {
    cloc = MatrixXd::Random(B, 2);
  }
  if (alpha.size() == 0) {
    alpha = VectorXd::Random(B);
  }
  int dim_cond = cloc.cols();
  MatrixXd cloc_save(B, dim_cond, n_iter - n_burnin);
  MatrixXd alpha_save(B, n_iter - n_burnin);
  MatrixXd beta_save;
  if (Z.size() > 0) {
    if (beta.size() == 0) {
      beta = VectorXd::Zero(Z.cols());
    }
    beta_save.resize(Z.cols(), n_iter - n_burnin);
  }
  for (int k = 0; k < n_iter; ++k) {
    VectorXd y = G.colwise().flatten();
    for (int b = 0; b < B; ++b) {
      MatrixXd X = cloc.replicate(1, n_edge) * G.col(b).replicate(B, 1);
      X.block(b * n_edge, 0, n_edge, dim_cond) = (G * cloc).rowwise().sum() - cloc.row(b) * G.col(b).transpose();
      VectorXd hlp2;
      for (int bb = 0; bb < B; ++b) {
        if (bb == b) {
          hlp2.conservativeResize(hlp2.size() + n_edge);
          hlp2.tail(n_edge).setZero();
        } else {
          VectorXd hlp3 = (G.colwise().sum().asDiagonal() * cloc - cloc.row(b) * G.col(b).transpose() - cloc.row(bb) * G.col(bb).transpose()).rowwise().squaredNorm();
          hlp2.conservativeResize(hlp2.size() + hlp3.size());
          hlp2.tail(hlp3.size()) = hlp3;
        }
      }
      VectorXd offset = hlp2 + alpha.replicate(n_edge);
      if (Z.size() > 0) {
        offset = hlp2 + alpha.replicate(n_edge) + (G * beta.replicate(B, 1)).colwise().flatten();
      }
      cloc.row(b) = blr(y, X, offset, cloc.row(            b).transpose(), VectorXd::Zero(dim_cond)).transpose();
    }
    MatrixXd dist_cond(n_edge, B);
    for (int b = 0; b < B; ++b) {
      dist_cond.col(b) = (G * cloc).rowwise().squaredNorm() - cloc.row(b) * G.col(b).transpose();
      VectorXd offset = dist_cond.col(b);
      if (Z.size() > 0) {
        offset = dist_cond.col(b) + Z * beta;
      }
      VectorXd y = G.col(b);
      MatrixXd X = MatrixXd::Ones(n_edge, 1);
      alpha(b) = blr(y, X, offset, alpha(b), 0.0)(0);
    }
    if (k >= n_burnin) {
      cloc_save.block(0, 0, B, dim_cond, k - n_burnin) = cloc;
      alpha_save.col(k - n_burnin) = alpha;
      if (Z.size() > 0) {
        beta_save.col(k - n_burnin) = beta;
      }
    }
  }
}




