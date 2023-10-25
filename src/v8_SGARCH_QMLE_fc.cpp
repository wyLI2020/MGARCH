#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

////////////////////////////////////////////////// Estimation \widehat\bm\vartheta //////////////////////////////////////////////////
////// Under low-rank constraints rank(G_{0k}) = 1 and rank(G_{1k}) = rank(G_{2k}) = 2, 
////// \bm\vartheta = (\bm\delta', \bm\beta')'
//////              = (\underline\bm\omega'(m*1), \bm\kappa', \beta_1, \beta_2, \underline\bmr'(m(m-1)/2*1))'
//////              = (\underline\bm\omega'(m*1), \bm\lambda'(r*1), \bm\gamma'(s*1), \bm\varphi'(s*1), \bm{g}_0'(2rm*1), \bm{g}_1'(4sm*1), \bm{g}_2'(4sm*1), \beta_1, \beta_2, \underline\bm{r}'(m(m-1)/2*1))'
////// \bm{g}_0 = (\bm{g}_{01}', ..., \bm{g}_{0r}')', \bm{g}_1 = (\bm{g}_{11}', ..., \bm{g}_{1s}')', \bm{g}_2 = (\bm{g}_{21}', ..., \bm{g}_{2s}')' with 
////// \bm{g}_{0k} = (\bm{g}_{0k1}', \bm{g}_{0k2}')' with G_0k = \bm{g}_{0k1} \bm{g}_{0k2}' such that rank(G_0k) = 1, 
////// \bm{g}_{1k} = (\bm{g}_{1k1}', \bm{g}_{1k2}', \bm{g}_{1k3}', \bm{g}_{1k4}')' with G_1k = \bm{g}_{1k1} \bm{g}_{1k2}' + \bm{g}_{1k3} \bm{g}_{1k4}' such that rank(G_1k) <= 2, 
////// \bm{g}_{2k} = (\bm{g}_{2k1}', \bm{g}_{2k2}', \bm{g}_{2k3}', \bm{g}_{2k4}')' with G_2k = \bm{g}_{2k1} \bm{g}_{2k2}' + \bm{g}_{2k3} \bm{g}_{2k4}' such that rank(G_2k) <= 2, 
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////// constraints of parameters //////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat fc_lowrank_linear_constr_ui(int m, int r, int s, int dimvartheta) {
  // linear inequality constraints, ui %*% theta - ci >=0
  int dimruline = m*(m-1)/2;
  int dimchi = dimvartheta - dimruline;
  arma::vec dim_inequality_lambda_choices(2); dim_inequality_lambda_choices(0) = r-1; dim_inequality_lambda_choices(1) = 0;
  int dim_inequality_lambda = max(dim_inequality_lambda_choices);
  int dim_constr_lambda = 2*r + 4*s + dim_inequality_lambda;
  int dim_constr_beta12 = 5;
  int dim_constr_ruline = 2*dimruline;
  arma::mat ui(dim_constr_lambda+dim_constr_beta12+dim_constr_ruline, dimvartheta); ui.fill(0.0);
  if(s == 0) {
    for(int k = 0; k < r; k++) {
      ui(2*k, m+k) = -1.0;
      ui(2*k+1, m+k) = 1.0;
    }
    for(int k = 1; k < r; k++) {
      ui(2*r+4*s+k-1, m+k-1) = 1.0; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
      ui(2*r+4*s+k-1, m+k) = -1.0; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
    }
  } else if(r == 0) {
    for(int k = 0; k < s; k++) {
      ui(2*r+2*k, m+r+k) = -1.0;
      ui(2*r+2*k+1, m+r+k) = 1.0;
      ui(2*r+2*s+2*k, m+r+s+k) = -1.0;
      ui(2*r+2*s+2*k+1, m+r+s+k) = 1.0;
    }
  } else {
    for(int k = 0; k < r; k++) {
      ui(2*k, m+k) = -1.0;
      ui(2*k+1, m+k) = 1.0;
    }
    for(int k = 0; k < s; k++) {
      ui(2*r+2*k, m+r+k) = -1.0;
      ui(2*r+2*k+1, m+r+k) = 1.0;
      ui(2*r+2*s+2*k, m+r+s+k) = -1.0;
      ui(2*r+2*s+2*k+1, m+r+s+k) = 1.0;
    }
    for(int k = 1; k < r; k++) {
      ui(2*r+4*s+k-1, m+k-1) = 1.0; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
      ui(2*r+4*s+k-1, m+k) = -1.0; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
    }
  }
  ui(dim_constr_lambda, dimchi-2) = -1.0;
  ui(dim_constr_lambda+1, dimchi-2) = 1.0;
  ui(dim_constr_lambda+2, dimchi-1) = -1.0;
  ui(dim_constr_lambda+3, dimchi-1) = 1.0;
  ui(dim_constr_lambda+4, dimchi-2) = -1.0; // \beta_1 + \beta_2 <= 0.999
  ui(dim_constr_lambda+4, dimchi-1) = -1.0; // \beta_1 + \beta_2 <= 0.999
  for(int k = 0; k < dimruline; k++) {
    ui(dim_constr_lambda+5+2*k, dimchi+k) = -1.0;
    ui(dim_constr_lambda+5+2*k+1, dimchi+k) = 1.0;
  }
  return ui;
}

// [[Rcpp::export]]
arma::vec fc_lowrank_linear_constr_ci(int m, int r, int s) {
  // linear inequality constraints, ui %*% theta - ci >=0
  int dimruline = m*(m-1)/2;
  arma::vec dim_inequality_lambda_choices(2); dim_inequality_lambda_choices(0) = r-1; dim_inequality_lambda_choices(1) = 0;
  int dim_inequality_lambda = max(dim_inequality_lambda_choices);
  int dim_constr_lambda = 2*r + 4*s + dim_inequality_lambda;
  int dim_constr_beta12 = 5;
  int dim_constr_ruline = 2*dimruline;
  arma::mat ci_mat(dim_constr_lambda+dim_constr_beta12+dim_constr_ruline, 1); ci_mat.fill(0.0);
  arma::vec ci = ci_mat.col(0);
  if(s == 0) {
    for(int k = 0; k < r; k++) {
      ci(2*k) = -0.999;
      ci(2*k+1) = -0.999;
    }
    for(int k = 1; k < r; k++) {
      ci(2*r+4*s+k-1) = 0.001; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
    }
  } else if(r == 0) {
    for(int k = 0; k < s; k++) {
      ci(2*r+2*k) = -0.999;
      ci(2*r+2*k+1) = 0.001;
      ci(2*r+2*s+2*k) = -3.141;
      ci(2*r+2*s+2*k+1) = 0.001;
    }
  } else {
    for(int k = 0; k < r; k++) {
      ci(2*k) = -0.999;
      ci(2*k+1) = -0.999;
    }
    for(int k = 0; k < s; k++) {
      ci(2*r+2*k) = -0.999;
      ci(2*r+2*k+1) = 0.001;
      ci(2*r+2*s+2*k) = -3.141;
      ci(2*r+2*s+2*k+1) = 0.001;
    }
    for(int k = 1; k < r; k++) {
      ci(2*r+4*s+k-1) = 0.001; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
    }
  }
  ci(dim_constr_lambda) = -0.999;
  ci(dim_constr_lambda+1) = 0.001;
  ci(dim_constr_lambda+2) = -0.999;
  ci(dim_constr_lambda+3) = 0.001;
  ci(dim_constr_lambda+4) = -0.999; // \beta_1 + \beta_2 <= 0.999
  for(int k = 0; k < dimruline; k++) {
    ci(dim_constr_lambda+5+2*k) = -0.999;
    ci(dim_constr_lambda+5+2*k+1) = -0.999;
  }
  return ci;
}

////////////////////////////////////////////////// summations using FFT algorithm //////////////////////////////////////////////////

arma::vec fc_FFT(double cst, arma::vec lnyuline_nminus1_l, arma::vec a) {
  // FFT-based method
  Function f("fr_FFT");
  NumericVector sum_iT = f(cst, lnyuline_nminus1_l, a);
  arma::vec sum_iT_arma(sum_iT.begin(), sum_iT.size(), false);
  return sum_iT_arma;
}

arma::Col<int> fc_seq(int a, int b) {
  Rcpp::IntegerVector seq_a_b_1 = seq(a, b);
  arma::Col<int> seq_a_b_arma(seq_a_b_1.begin(), seq_a_b_1.size(), false);
  return seq_a_b_arma;
}

arma::vec fc_vpow(arma::vec base, arma::vec power) {
  NumericVector base_std(base.begin(), base.end());
  NumericVector power_std(power.begin(), power.end());
  NumericVector res_std(base.size());
  std::transform(base.begin(), base.end(), power.begin(), res_std.begin(),
                 [&](double lhs, double rhs) -> double {
                   return std::pow(lhs, rhs);
                 });
  arma::colvec res(res_std.begin(), res_std.size(), false);
  return res;
}

arma::cube fc_sum_lambdaklnyuline_0_m_nminus1_r(int n, int m, int r, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}, for k=1,...,r and t=2,...,n
  arma::cube sum_lambdaklnyuline_0_m_nminus1_r(m, n-1, r); sum_lambdaklnyuline_0_m_nminus1_r.fill(0.0);
  arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < r; k++) {
    double lambda_k = lambda_vec(k); arma::vec lambda_k_rep = rep(lambda_k, n-1);
    arma::mat sum_lambdaklnyuline_0_nminus1_m_k(n-1, m); sum_lambdaklnyuline_0_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = fc_vpow(lambda_k_rep, i_seq_1_nminus1-1);
      sum_lambdaklnyuline_0_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_lambdaklnyuline_0_m_nminus1_r.slice(k) = sum_lambdaklnyuline_0_nminus1_m_k.t();
  }
  return sum_lambdaklnyuline_0_m_nminus1_r;
}

arma::cube fc_sum_lambdaklnyuline_1_m_nminus1_r(int n, int m, int r, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}, for k=1,...,r and t=2,...,n
  arma::cube sum_lambdaklnyuline_1_m_nminus1_r(m, n-1, r); sum_lambdaklnyuline_1_m_nminus1_r.fill(0.0);
  arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < r; k++) {
    double lambda_k = lambda_vec(k); arma::vec lambda_k_rep = rep(lambda_k, n-1);
    arma::mat sum_lambdaklnyuline_1_nminus1_m_k(n-1, m); sum_lambdaklnyuline_1_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % fc_vpow(lambda_k_rep, i_seq_1_nminus1-2);
      sum_lambdaklnyuline_1_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_lambdaklnyuline_1_m_nminus1_r.slice(k) = sum_lambdaklnyuline_1_nminus1_m_k.t();
  }
  return sum_lambdaklnyuline_1_m_nminus1_r;
}

arma::cube fc_sum_lambdaklnyuline_2_m_nminus1_r(int n, int m, int r, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=3}^{t-1} (i-1) (i-2) \lambda_k^{i-3} \ln\underline{\bm{y}}_{t-i}, for k=1,...,r and t=2,...,n
  arma::cube sum_lambdaklnyuline_2_m_nminus1_r(m, n-1, r); sum_lambdaklnyuline_2_m_nminus1_r.fill(0.0);
  arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < r; k++) {
    double lambda_k = lambda_vec(k); arma::vec lambda_k_rep = rep(lambda_k, n-1);
    arma::mat sum_lambdaklnyuline_2_nminus1_m_k(n-1, m); sum_lambdaklnyuline_2_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % (i_seq_1_nminus1-2) % fc_vpow(lambda_k_rep, i_seq_1_nminus1-3);
      sum_lambdaklnyuline_2_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_lambdaklnyuline_2_m_nminus1_r.slice(k) = sum_lambdaklnyuline_2_nminus1_m_k.t();
  }
  return sum_lambdaklnyuline_2_m_nminus1_r;
}

arma::cube fc_sum_gammaklnyuline_01_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_01_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_01_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_01_nminus1_m_k(n-1, m); sum_gammaklnyuline_01_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = fc_vpow(gamma_k_rep, i_seq_1_nminus1-1) % cos((i_seq_1_nminus1-1) * varphi_k);
      sum_gammaklnyuline_01_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_01_m_nminus1_s.slice(k) = sum_gammaklnyuline_01_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_01_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_02_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_02_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_02_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_02_nminus1_m_k(n-1, m); sum_gammaklnyuline_02_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = fc_vpow(gamma_k_rep, i_seq_1_nminus1-1) % sin((i_seq_1_nminus1-1) * varphi_k);
      sum_gammaklnyuline_02_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_02_m_nminus1_s.slice(k) = sum_gammaklnyuline_02_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_02_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_11_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_11_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_11_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_11_nminus1_m_k(n-1, m); sum_gammaklnyuline_11_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % fc_vpow(gamma_k_rep, i_seq_1_nminus1-2) % cos((i_seq_1_nminus1-1) * varphi_k);
      sum_gammaklnyuline_11_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_11_m_nminus1_s.slice(k) = sum_gammaklnyuline_11_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_11_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_12_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_12_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_12_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_12_nminus1_m_k(n-1, m); sum_gammaklnyuline_12_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % fc_vpow(gamma_k_rep, i_seq_1_nminus1-2) % sin((i_seq_1_nminus1-1) * varphi_k);
      sum_gammaklnyuline_12_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_12_m_nminus1_s.slice(k) = sum_gammaklnyuline_12_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_12_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_21_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_21_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_21_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_21_nminus1_m_k(n-1, m); sum_gammaklnyuline_21_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % fc_vpow(gamma_k_rep, i_seq_1_nminus1-1) % (-sin((i_seq_1_nminus1-1) * varphi_k));
      sum_gammaklnyuline_21_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_21_m_nminus1_s.slice(k) = sum_gammaklnyuline_21_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_21_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_22_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_22_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_22_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_22_nminus1_m_k(n-1, m); sum_gammaklnyuline_22_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % fc_vpow(gamma_k_rep, i_seq_1_nminus1-1) % cos((i_seq_1_nminus1-1) * varphi_k);
      sum_gammaklnyuline_22_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_22_m_nminus1_s.slice(k) = sum_gammaklnyuline_22_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_22_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_31_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_31_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_31_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_31_nminus1_m_k(n-1, m); sum_gammaklnyuline_31_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % (i_seq_1_nminus1-1) % fc_vpow(gamma_k_rep, i_seq_1_nminus1-2) % (-sin((i_seq_1_nminus1-1) * varphi_k));
      sum_gammaklnyuline_31_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_31_m_nminus1_s.slice(k) = sum_gammaklnyuline_31_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_31_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_32_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_32_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_32_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_32_nminus1_m_k(n-1, m); sum_gammaklnyuline_32_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % (i_seq_1_nminus1-1) % fc_vpow(gamma_k_rep, i_seq_1_nminus1-2) % cos((i_seq_1_nminus1-1) * varphi_k);
      sum_gammaklnyuline_32_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_32_m_nminus1_s.slice(k) = sum_gammaklnyuline_32_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_32_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_41_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\cos((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_41_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_41_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_41_nminus1_m_k(n-1, m); sum_gammaklnyuline_41_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % (i_seq_1_nminus1-1) % fc_vpow(gamma_k_rep, i_seq_1_nminus1-1) % (-cos((i_seq_1_nminus1-1) * varphi_k));
      sum_gammaklnyuline_41_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_41_m_nminus1_s.slice(k) = sum_gammaklnyuline_41_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_41_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_42_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_42_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_42_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_42_nminus1_m_k(n-1, m); sum_gammaklnyuline_42_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % (i_seq_1_nminus1-1) % fc_vpow(gamma_k_rep, i_seq_1_nminus1-1) % (-sin((i_seq_1_nminus1-1) * varphi_k));
      sum_gammaklnyuline_42_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_42_m_nminus1_s.slice(k) = sum_gammaklnyuline_42_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_42_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_51_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_51_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_51_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_51_nminus1_m_k(n-1, m); sum_gammaklnyuline_51_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % (i_seq_1_nminus1-2) % fc_vpow(gamma_k_rep, i_seq_1_nminus1-3) % cos((i_seq_1_nminus1-1) * varphi_k);
      sum_gammaklnyuline_51_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_51_m_nminus1_s.slice(k) = sum_gammaklnyuline_51_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_51_m_nminus1_s;
}

arma::cube fc_sum_gammaklnyuline_52_m_nminus1_s(int n, int m, int r, int s, arma::vec kappa_vec, arma::mat y_m_n) {
  // the (t-1)-th column of the k-th slice is \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}, for k=1,...,s and t=2,...,n
  arma::cube sum_gammaklnyuline_52_m_nminus1_s(m, n-1, s); sum_gammaklnyuline_52_m_nminus1_s.fill(0.0);
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::mat y_m_nminus1 = y_m_n.cols(0, n-2);
  arma::mat ty_nminus1_m = y_m_nminus1.t();
  arma::Col<int> i_seq_1_nminus1_int = fc_seq(1, n-1); arma::vec onesvec_nminus1 = rep(1.0, n-1); arma::vec i_seq_1_nminus1 = i_seq_1_nminus1_int % onesvec_nminus1;
  for(int k = 0; k < s; k++) {
    double gamma_k = gamma_vec(k); arma::vec gamma_k_rep = rep(gamma_k, n-1);
    double varphi_k = varphi_vec(k); arma::vec varphi_k_rep = rep(varphi_k, n-1);
    arma::mat sum_gammaklnyuline_52_nminus1_m_k(n-1, m); sum_gammaklnyuline_52_nminus1_m_k.fill(0.0);
    for(int l = 0; l < m; l++) {
      arma::vec y_nminus1_l = ty_nminus1_m.col(l);
      arma::vec lnyuline_nminus1_l = log(pow(y_nminus1_l, 2));
      arma::vec a = (i_seq_1_nminus1-1) % (i_seq_1_nminus1-2) % fc_vpow(gamma_k_rep, i_seq_1_nminus1-3) % sin((i_seq_1_nminus1-1) * varphi_k);
      sum_gammaklnyuline_52_nminus1_m_k.col(l) = fc_FFT(0.0, lnyuline_nminus1_l, a);
    }
    sum_gammaklnyuline_52_m_nminus1_s.slice(k) = sum_gammaklnyuline_52_nminus1_m_k.t();
  }
  return sum_gammaklnyuline_52_m_nminus1_s;
}

////////////////////////////////////////////////// loss function //////////////////////////////////////////////////

arma::mat fc_asmat(arma::vec vec1, int nrow, int ncol) {
  // Fill matrix with elements of vector
  arma::mat vec1_mat(nrow*ncol, 1); vec1_mat.col(0) = vec1;
  vec1_mat.reshape(nrow, ncol);
  return vec1_mat;
}

arma::vec fc_asvec(arma::mat mat1) {
  // Matrix straighten
  int nrow = mat1.n_rows;
  int ncol = mat1.n_cols;
  mat1.reshape(nrow*ncol, 1);
  return mat1.col(0);
}

arma::mat fc_Ruline(int m, arma::vec ruline) {
  // \underline{R}
  arma::mat Ruline(m, m); Ruline.fill(1.0);
  int l = 0;
  for(int j = 0; j < (m-1); j++) { // \underline{R}_{2,1}, ..., \underline{R}_{m,1}; \underline{R}_{3,2}, ..., \underline{R}_{m,2}; ...; \underline{R}_{m,m-1}
    for(int i = (j+1); i < m; i++) {
      Ruline(i,j) = Ruline(j,i) = ruline(l);
      l = l + 1;
    }
  }
  return Ruline;
}

// arma::mat fc_lowrank_Phi_i(int i, int m, int r, int s, arma::vec kappa_vec) {
//   // r>0 and s>0, \Phi_i(\bm\kappa) = \sum_{k=1}^r \lambda_k^{i-1} G_{0,k} + \sum_{k=1}^s \gamma_k^{i-1} [\cos((i-1)\varphi_k) G_{1,k} + \sin((i-1)\varphi_k) G_{2,k}]
//   arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
//   arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
//   arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
//   arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, 2*m, r);
//   arma::vec g_1_vec = kappa_vec.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, 4*m, s);
//   arma::vec g_2_vec = kappa_vec.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_mat = fc_asmat(g_2_vec, 4*m, s);
//   arma::mat sum_lambdaG0(m, m); sum_lambdaG0.fill(0.0);
//   arma::mat sum_gammaG12(m, m); sum_gammaG12.fill(0.0);
//   for(int k = 0; k < r; k++) {
//     double lambda_k = lambda_vec(k);
//     arma::vec g_0k = g_0_mat.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2); 
//     arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
//     arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
//     arma::mat G_0k = g_0k1 * g_0k2.t();
//     sum_lambdaG0 = sum_lambdaG0 + pow(lambda_k, i-1) * G_0k;
//   }
//   for(int k = 0; k < s; k++) {
//     double gamma_k = gamma_vec(k);
//     double varphi_k = varphi_vec(k);
//     arma::vec g_1k = g_1_mat.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4); 
//     arma::vec g_2k = g_2_mat.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4); 
//     arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
//     arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
//     arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
//     arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
//     arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
//     arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
//     arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
//     arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
//     arma::mat G_1k = g_1k1 * g_1k2.t() + g_1k3 * g_1k4.t();
//     arma::mat G_2k = g_2k1 * g_2k2.t() + g_2k3 * g_2k4.t();
//     sum_gammaG12 = sum_gammaG12 + pow(gamma_k, i-1) * (cos((i-1) * varphi_k) * G_1k + sin((i-1) * varphi_k) * G_2k);
//   }
//   arma::mat Phi_i = sum_lambdaG0 + sum_gammaG12; 
//   return Phi_i;
// }
// 
// arma::mat fc_lowrank_Phi_i_requal0(int i, int m, int r, int s, arma::vec kappa_vec) {
//   // r=0 and s>0, \Phi_i(\bm\kappa) = \sum_{k=1}^s \gamma_k^{i-1} [\cos((i-1)\varphi_k) G_{1,k} + \sin((i-1)\varphi_k) G_{2,k}]
//   // arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
//   arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
//   arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
//   // arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, 2*m, r);
//   arma::vec g_1_vec = kappa_vec.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, 4*m, s);
//   arma::vec g_2_vec = kappa_vec.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_mat = fc_asmat(g_2_vec, 4*m, s);
//   arma::mat sum_lambdaG0(m, m); sum_lambdaG0.fill(0.0);
//   arma::mat sum_gammaG12(m, m); sum_gammaG12.fill(0.0);
//   // for(int k = 0; k < r; k++) {
//   //   double lambda_k = lambda_vec(k);
//   //   arma::vec g_0k = g_0_mat.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2); 
//   //   arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
//   //   arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
//   //   arma::mat G_0k = g_0k1 * g_0k2.t();
//   //   sum_lambdaG0 = sum_lambdaG0 + pow(lambda_k, i-1) * G_0k;
//   // }
//   for(int k = 0; k < s; k++) {
//     double gamma_k = gamma_vec(k);
//     double varphi_k = varphi_vec(k);
//     arma::vec g_1k = g_1_mat.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4); 
//     arma::vec g_2k = g_2_mat.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4); 
//     arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
//     arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
//     arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
//     arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
//     arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
//     arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
//     arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
//     arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
//     arma::mat G_1k = g_1k1 * g_1k2.t() + g_1k3 * g_1k4.t();
//     arma::mat G_2k = g_2k1 * g_2k2.t() + g_2k3 * g_2k4.t();
//     sum_gammaG12 = sum_gammaG12 + pow(gamma_k, i-1) * (cos((i-1) * varphi_k) * G_1k + sin((i-1) * varphi_k) * G_2k);
//   }
//   arma::mat Phi_i = sum_lambdaG0 + sum_gammaG12; 
//   return Phi_i;
// }
// 
// arma::mat fc_lowrank_Phi_i_sequal0(int i, int m, int r, int s, arma::vec kappa_vec) {
//   // r>0 and s=0, \Phi_i(\bm\kappa) = \sum_{k=1}^r \lambda_k^{i-1} G_{0,k}
//   arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
//   // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
//   // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
//   arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, 2*m, r);
//   // arma::vec g_1_vec = kappa_vec.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, 4*m, s);
//   // arma::vec g_2_vec = kappa_vec.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_mat = fc_asmat(g_2_vec, 4*m, s);
//   arma::mat sum_lambdaG0(m, m); sum_lambdaG0.fill(0.0);
//   arma::mat sum_gammaG12(m, m); sum_gammaG12.fill(0.0);
//   for(int k = 0; k < r; k++) {
//     double lambda_k = lambda_vec(k);
//     arma::vec g_0k = g_0_mat.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2); 
//     arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
//     arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
//     arma::mat G_0k = g_0k1 * g_0k2.t();
//     sum_lambdaG0 = sum_lambdaG0 + pow(lambda_k, i-1) * G_0k;
//   }
//   // for(int k = 0; k < s; k++) {
//   //   double gamma_k = gamma_vec(k);
//   //   double varphi_k = varphi_vec(k);
//   //   arma::vec g_1k = g_1_mat.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4); 
//   //   arma::vec g_2k = g_2_mat.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4); 
//   //   arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
//   //   arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
//   //   arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
//   //   arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
//   //   arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
//   //   arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
//   //   arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
//   //   arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
//   //   arma::mat G_1k = g_1k1 * g_1k2.t() + g_1k3 * g_1k4.t();
//   //   arma::mat G_2k = g_2k1 * g_2k2.t() + g_2k3 * g_2k4.t();
//   //   sum_gammaG12 = sum_gammaG12 + pow(gamma_k, i-1) * (cos((i-1) * varphi_k) * G_1k + sin((i-1) * varphi_k) * G_2k);
//   // }
//   arma::mat Phi_i = sum_lambdaG0 + sum_gammaG12; 
//   return Phi_i;
// }
// 
// arma::mat fc_lowrank_Phi_i_general(int i, int m, int r, int s, arma::vec kappa_vec) {
//   // \Phi_i(\bm\kappa) = \sum_{k=1}^r \lambda_k^{i-1} G_{0,k} + \sum_{k=1}^s \gamma_k^{i-1} [\cos((i-1)\varphi_k) G_{1,k} + \sin((i-1)\varphi_k) G_{2,k}]
//   arma::mat Phi_i(m, m); Phi_i.fill(0.0);
//   if(s == 0) {
//     Phi_i = fc_lowrank_Phi_i_sequal0(i, m, r, s, kappa_vec);
//   } else if(r == 0) {
//     Phi_i = fc_lowrank_Phi_i_requal0(i, m, r, s, kappa_vec);
//   } else {
//     Phi_i = fc_lowrank_Phi_i(i, m, r, s, kappa_vec);
//   }
//   return Phi_i;
// }

// arma::vec fc_lowrank_lnhulinetilde_t(int t, int m, int r, int s, arma::vec omegauline, arma::vec kappa_vec, arma::mat y_m_n) {
//   // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = (\ln\widetilde{h}_{11,t}(\bm\delta), ..., \ln\widetilde{h}_{mm,t}(\bm\delta))'
//   arma::mat Philnyuline_mat(m, t-1); Philnyuline_mat.fill(0.0); // (\Phi_1(\bm\kappa) \ln\underline{\bm{y}}_{t-1}, ..., \Phi_{t-1}(\bm\kappa) \ln\underline{\bm{y}}_{1})
//   for(int i = 1; i < t; i++) {
//     arma::mat Phi_i = fc_lowrank_Phi_i_general(i, m, r, s, kappa_vec); // \Phi_i(\bm\kappa)
//     arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
//     arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln y_{1,t-i}^2, ..., \ln y_{m,t-i}^2)'
//     Philnyuline_mat.col(i-1) = Phi_i * lnyuline_tminusi; // \Phi_i(\bm\kappa) \ln\underline{\bm{y}}_{t-i}
//   }
//   arma::vec lnhulinetilde_t = omegauline + sum(Philnyuline_mat, 1); // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = \underline{\bm{\omega}} + \sum_{i=1}^{t-1} \Phi_i(\bm\kappa) \ln\underline{\bm{y}}_{t-i}
//   return lnhulinetilde_t;
// }

arma::vec fc_lowrank_lnhulinetilde_t(int t, int m, int r, int s, arma::vec omegauline, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_0_m_nminus1_r, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s) {
  // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = (\ln\widetilde{h}_{11,t}(\bm\delta), ..., \ln\widetilde{h}_{mm,t}(\bm\delta))'
  arma::vec lnhulinetilde_t = omegauline; // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = \underline{\bm{\omega}} + \sum_{i=1}^{t-1} \Phi_i(\bm\kappa) \ln\underline{\bm{y}}_{t-i}
  // arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
  // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, 2*m, r);
  arma::vec g_1_vec = kappa_vec.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, 4*m, s);
  arma::vec g_2_vec = kappa_vec.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_mat = fc_asmat(g_2_vec, 4*m, s);
  arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  for(int k = 0; k < r; k++) {
    // double lambda_k = lambda_vec(k);
    arma::vec g_0k = g_0_mat.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2); 
    arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
    arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
    arma::mat G_0k = g_0k1 * g_0k2.t();
    arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    lnhulinetilde_t = lnhulinetilde_t + G_0k * sum_lambdaklnyuline_0;
  }
  for(int k = 0; k < s; k++) {
    // double gamma_k = gamma_vec(k);
    // double varphi_k = varphi_vec(k);
    arma::vec g_1k = g_1_mat.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4); 
    arma::vec g_2k = g_2_mat.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4); 
    arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
    arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
    arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
    arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
    arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
    arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
    arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
    arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
    arma::mat G_1k = g_1k1 * g_1k2.t() + g_1k3 * g_1k4.t();
    arma::mat G_2k = g_2k1 * g_2k2.t() + g_2k3 * g_2k4.t();
    arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    lnhulinetilde_t = lnhulinetilde_t + G_1k * sum_gammaklnyuline_01 + G_2k * sum_gammaklnyuline_02;
  }
  return lnhulinetilde_t;
}

arma::vec fc_lowrank_lnhulinetilde_t_requal0(int t, int m, int r, int s, arma::vec omegauline, arma::vec kappa_vec, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s) {
  // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = (\ln\widetilde{h}_{11,t}(\bm\delta), ..., \ln\widetilde{h}_{mm,t}(\bm\delta))'
  arma::vec lnhulinetilde_t = omegauline; // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = \underline{\bm{\omega}} + \sum_{i=1}^{t-1} \Phi_i(\bm\kappa) \ln\underline{\bm{y}}_{t-i}
  // arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
  // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  // arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, 2*m, r);
  arma::vec g_1_vec = kappa_vec.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, 4*m, s);
  arma::vec g_2_vec = kappa_vec.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_mat = fc_asmat(g_2_vec, 4*m, s);
  // arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  // for(int k = 0; k < r; k++) {
  //   // double lambda_k = lambda_vec(k);
  //   arma::vec g_0k = g_0_mat.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2); 
  //   arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
  //   arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
  //   arma::mat G_0k = g_0k1 * g_0k2.t();
  //   arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
  //   lnhulinetilde_t = lnhulinetilde_t + G_0k * sum_lambdaklnyuline_0;
  // }
  for(int k = 0; k < s; k++) {
    // double gamma_k = gamma_vec(k);
    // double varphi_k = varphi_vec(k);
    arma::vec g_1k = g_1_mat.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4); 
    arma::vec g_2k = g_2_mat.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4); 
    arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
    arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
    arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
    arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
    arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
    arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
    arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
    arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
    arma::mat G_1k = g_1k1 * g_1k2.t() + g_1k3 * g_1k4.t();
    arma::mat G_2k = g_2k1 * g_2k2.t() + g_2k3 * g_2k4.t();
    arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    lnhulinetilde_t = lnhulinetilde_t + G_1k * sum_gammaklnyuline_01 + G_2k * sum_gammaklnyuline_02;
  }
  return lnhulinetilde_t;
}

arma::vec fc_lowrank_lnhulinetilde_t_sequal0(int t, int m, int r, int s, arma::vec omegauline, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_0_m_nminus1_r) {
  // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = (\ln\widetilde{h}_{11,t}(\bm\delta), ..., \ln\widetilde{h}_{mm,t}(\bm\delta))'
  arma::vec lnhulinetilde_t = omegauline; // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = \underline{\bm{\omega}} + \sum_{i=1}^{t-1} \Phi_i(\bm\kappa) \ln\underline{\bm{y}}_{t-i}
  // arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
  // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, 2*m, r);
  // arma::vec g_1_vec = kappa_vec.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, 4*m, s);
  // arma::vec g_2_vec = kappa_vec.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_mat = fc_asmat(g_2_vec, 4*m, s);
  arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  // arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  for(int k = 0; k < r; k++) {
    // double lambda_k = lambda_vec(k);
    arma::vec g_0k = g_0_mat.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2); 
    arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
    arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
    arma::mat G_0k = g_0k1 * g_0k2.t();
    arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    lnhulinetilde_t = lnhulinetilde_t + G_0k * sum_lambdaklnyuline_0;
  }
  // for(int k = 0; k < s; k++) {
  //   // double gamma_k = gamma_vec(k);
  //   // double varphi_k = varphi_vec(k);
  //   arma::vec g_1k = g_1_mat.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4); 
  //   arma::vec g_2k = g_2_mat.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4); 
  //   arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
  //   arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
  //   arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
  //   arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
  //   arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
  //   arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
  //   arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
  //   arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
  //   arma::mat G_1k = g_1k1 * g_1k2.t() + g_1k3 * g_1k4.t();
  //   arma::mat G_2k = g_2k1 * g_2k2.t() + g_2k3 * g_2k4.t();
  //   arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   lnhulinetilde_t = lnhulinetilde_t + G_1k * sum_gammaklnyuline_01 + G_2k * sum_gammaklnyuline_02;
  // }
  return lnhulinetilde_t;
}

arma::mat fc_lowrank_lnhulinetilde_t_general(int t, int m, int r, int s, arma::vec omegauline, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_0_m_nminus1_r, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s) {
  // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = (\ln\widetilde{h}_{11,t}(\bm\delta), ..., \ln\widetilde{h}_{mm,t}(\bm\delta))'
  arma::vec lnhulinetilde_t(m);
  if(s == 0) {
    lnhulinetilde_t = fc_lowrank_lnhulinetilde_t_sequal0(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r);
  } else if(r == 0) {
    lnhulinetilde_t = fc_lowrank_lnhulinetilde_t_requal0(t, m, r, s, omegauline, kappa_vec, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s);
  } else {
    lnhulinetilde_t = fc_lowrank_lnhulinetilde_t(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s);
  }
  return lnhulinetilde_t;
}

arma::mat fc_Dtilde_1(arma::vec omegauline) {
  // \widetilde{D}_{1}(\bm\delta)
  arma::mat Dtilde_1 = diagmat(exp(0.5*omegauline),0); // \ln\underline{\bm{h}}_1 = \underline{\bm{\omega}} and then \widetilde{D}_{1}(\bm\delta) = \Diag{\exp{1/2*\underline{\omega}_1}, ..., \exp{1/2*\underline{\omega}_m}}
  return Dtilde_1;
}

arma::mat fc_inverse_Dtilde_1(arma::vec omegauline) {
  // \widetilde{D}_{1}(\bm\delta)
  arma::mat inverse_Dtilde_1 = diagmat(exp(-0.5*omegauline),0); // \widetilde{D}_{1}^{-1}(\bm\delta)
  return inverse_Dtilde_1;
}

arma::mat fc_Dtilde_t(arma::vec lnhulinetilde_t) {
  // \widetilde{D}_{t}(\bm\delta)
  arma::vec hulinetilde_t = exp(lnhulinetilde_t); // \widetilde{h}_{ii,t}(\bm\delta) = \exp{\ln\widetilde{h}_{ii,t}(\bm\delta)}
  arma::mat Dtilde_t = diagmat(sqrt(hulinetilde_t),0); // \widetilde{D}_{t}(\bm\delta) = \Diag{\sqrt{\widetilde{h}_{11,t}(\bm\delta)}, ..., \sqrt{\widetilde{h}_{mm,t}(\bm\delta)}}
  return Dtilde_t;
}

arma::mat fc_inverse_Dtilde_t(arma::vec lnhulinetilde_t) {
  // \widetilde{D}_{t}(\bm\delta)
  arma::vec hulinetilde_t = exp(lnhulinetilde_t); // \widetilde{h}_{ii,t}(\bm\delta) = \exp{\ln\widetilde{h}_{ii,t}(\bm\delta)}
  arma::mat inverse_Dtilde_t = diagmat(1.0/sqrt(hulinetilde_t),0); // \widetilde{D}_{t}^{-1}(\bm\delta)
  return inverse_Dtilde_t;
}

arma::mat fc_Psi(int m, arma::mat x_m_Bk) {
  // the sample correlation matrix of random vectors {x_1,...,x_{\Bbbk}}
  arma::mat Psi(m, m); Psi.fill(1.0);
  for(int i = 1; i < m; i++) { // Psi_{2,1}; Psi_{3,1}, Psi_{3,2}; ...; Psi_{m,1}, ..., Psi_{m,m-1}
    for(int j = 0; j < i; j++) {
      arma::rowvec x_ith_Bk = x_m_Bk.row(i); 
      arma::rowvec x_jth_Bk = x_m_Bk.row(j); 
      double numerator = dot(x_ith_Bk, x_jth_Bk);
      double denominator = sqrt(dot(x_ith_Bk, x_ith_Bk) * dot(x_jth_Bk, x_jth_Bk));
      Psi(i,j) = Psi(j,i) = numerator / denominator;
    }
  }
  return Psi;
}

arma::mat fc_Rtilde_1(int m, double beta_1, double beta_2, arma::mat Ruline) {
  // \widetilde{R}_1(\bm\theta)
  arma::mat onesmat_m(m, m); onesmat_m.fill(1.0); // m*m matrix of ones
  arma::mat Rtilde_1 = (1.0 - beta_1 - beta_2) / (1.0 - beta_2) * Ruline + beta_1 / (1.0 - beta_2) * onesmat_m; // \widetilde{R}_{1}(\bm\theta) = (1 - \beta_1 - \beta_2) / (1 - \beta_2) \underline{R} + \beta_1 / (1 - \beta_2) 1_m
  return Rtilde_1;
}

arma::mat fc_Rtilde_t(int m, int Bk, double beta_1, double beta_2, arma::mat Ruline, arma::mat Psitilde_tminus1, arma::mat Rtilde_tminus1) {
  // \widetilde{R}_t(\bm\theta)
  arma::mat Rtilde_t = (1.0 - beta_1 - beta_2) * Ruline + beta_1 * Psitilde_tminus1 + beta_2 * Rtilde_tminus1; // \widetilde{R}_t(\bm\theta) = (1-\beta_1-\beta_2) \underline{R} + \beta_1 \widetilde{\Psi}_{t-1}(\bm\delta) + \beta_2 \widetilde{R}_{t-1}(\bm\theta)
  return Rtilde_t;
}

arma::mat fc_Htilde_t(arma::mat Dtilde_t, arma::mat Rtilde_t) {
  // \widetilde{H}_t(\bm\theta)
  arma::mat Htilde_t = Dtilde_t * Rtilde_t * Dtilde_t; // \widetilde{H}_t(\bm\theta) = \widetilde{D}_t(\bm\delta) \widetilde{R}_t(\bm\theta) \widetilde{D}_t(\bm\delta)
  return Htilde_t;
}

double fc_ltilde_t(arma::vec y_t, arma::mat Htilde_t) { 
  // \widetilde{\ell}_{t}(\bm\theta)
  arma::mat inverse_Htilde_t = Htilde_t.i(); // \widetilde{H}_t^{-1}(\bm\theta)
  // double logdet_Htilde_t, sign_logdet; log_det(logdet_Htilde_t, sign_logdet, Htilde_t); // \ln|\widetilde{H}_t(\bm\theta)| 
  double det_Htilde_t = det(Htilde_t); 
  double logdet_Htilde_t = log(det_Htilde_t); // \ln|\widetilde{H}_t(\bm\theta)| 
  double ltilde_t = dot(y_t, inverse_Htilde_t * y_t) + logdet_Htilde_t; // \widetilde{\ell}_{t}(\bm\theta) = \bm{y}_t' \widetilde{H}_t^{-1}(\bm\theta) \bm{y}_t + \ln|\widetilde{H}_t(\bm\theta)|
  return ltilde_t;
}

// [[Rcpp::export]]
double fc_lowrank_mathcalLtilde_n(int n, int m, int r, int s, int Bk, arma::vec vartheta_vec, arma::mat y_m_n) {
  // \widetilde{\mathcal{L}}_{n}(\bm\theta), under initial values \widetilde{y}_s = 1_m for s <= 0
  int dimkappa = r + 2*s + 2*r*m + 8*s*m;
  int dimbeta = 2 + m*(m-1)/2;
  arma::vec omegauline = vartheta_vec.subvec(0, m-1); // \underline{\bm\omega}
  arma::vec kappa_vec = vartheta_vec.subvec(m, m+dimkappa-1); // \bm\kappa
  arma::vec beta_vec = vartheta_vec.subvec(m+dimkappa, m+dimkappa+dimbeta-1); // \bm\beta
  double beta_1 = beta_vec(0); double beta_2 = beta_vec(1); // \beta_1 and \beta_2
  arma::vec ruline = beta_vec.tail(dimbeta-2); // \underline\bm{r}
  arma::mat Ruline = fc_Ruline(m, ruline); // \underline{R}
  // summations using FFT algorithm
  int max_r_1 = r; if(max_r_1 == 0) {max_r_1 = 1;}
  int max_s_1 = s; if(max_s_1 == 0) {max_s_1 = 1;}
  arma::cube sum_lambdaklnyuline_0_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_0_m_nminus1_r.fill(0.0);
  arma::cube sum_gammaklnyuline_01_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_01_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_02_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_02_m_nminus1_s.fill(0.0);
  if(s == 0) {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
  } else if(r == 0) {
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  } else {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  }
  // \widetilde{\ell}_{t}(\bm\theta) and \widetilde\mathcal{L}_n(\bm\theta)
  double sum_ltilde = 0.0; // \sum_{t=1}^{n} \widetilde{\ell}_{t}(\bm\theta)
  // t=1
  arma::vec y_1 = y_m_n.col(0); // \bm{y}_{1}
  arma::mat Dtilde_1 = fc_Dtilde_1(omegauline); // \widetilde{D}_{1}(\bm\delta)
  arma::mat inverse_Dtilde_1 = fc_inverse_Dtilde_1(omegauline); // \widetilde{D}_{1}^{-1}(\bm\delta)
  arma::mat Rtilde_1 = fc_Rtilde_1(m, beta_1, beta_2, Ruline); // \widetilde{R}_{1}(\bm\theta)
  arma::mat Htilde_1 = fc_Htilde_t(Dtilde_1, Rtilde_1); // \widetilde{H}_{1}(\bm\theta)
  arma::vec varepsilontilde_1 = inverse_Dtilde_1 * y_1; // \widetilde{\varepsilon}_{1}(\bm\delta) = \widetilde{D}_{1}^{-1}(\bm\delta) \bm{y}_{1}
  // \widetilde{\ell}_1(\bm\theta)
  double ltilde_1 = fc_ltilde_t(y_1, Htilde_1); // \widetilde{\ell}_{1}(\bm\theta)
  sum_ltilde = sum_ltilde + ltilde_1;
  // t=2,...,n
  // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
  arma::mat varepsilontilde_tminusBkminus1TOtminus2_mat(m, Bk); varepsilontilde_tminusBkminus1TOtminus2_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta))
  varepsilontilde_tminusBkminus1TOtminus2_mat.each_col() = exp(-0.5*omegauline); // for s <= 0, \widetilde{\varepsilon}_s(\bm\delta) = (\exp{-1/2*\underliner{\omega}_1}, ..., \exp{-1/2*\underliner{\omega}_m})'
  arma::vec varepsilontilde_tminus1 = varepsilontilde_1; // \widetilde{\varepsilon}_{t-1}(\bm\delta)
  // \widetilde{R}_{t-1}(\bm\theta)
  arma::mat Rtilde_tminus1 = Rtilde_1; // \widetilde{R}_{t-1}(\bm\theta)
  for(int t = 2; t < (n+1); t++) {
    // \bm{y}_2, ..., \bm{y}_n
    arma::vec y_t = y_m_n.col(t-1); // \bm{y}_{t}
    // \widetilde{D}_2(\bm\delta), ..., \widetilde{D}_n(\bm\delta)
    arma::vec lnhulinetilde_t = fc_lowrank_lnhulinetilde_t_general(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s); // \ln\widetilde\underline{\bm{h}}_t(\bm\delta)
    arma::mat Dtilde_t = fc_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}(\bm\delta)
    arma::mat inverse_Dtilde_t = fc_inverse_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}^{-1}(\bm\delta)
    // \widetilde{R}_2(\bm\theta), ..., \widetilde{R}_n(\bm\theta)
    arma::vec varepsilontilde_t = inverse_Dtilde_t * y_t; // \widetilde{\varepsilon}_{t}(\bm\delta) = \widetilde{D}_{t}^{-1}(\bm\delta) \bm{y}_{t}
    arma::mat varepsilontilde_tminusBkTOtminus1_mat(m, Bk); varepsilontilde_tminusBkTOtminus1_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk}(\bm\delta), ..., \widetilde{\varepsilon}_{t-1}(\bm\delta))
    varepsilontilde_tminusBkTOtminus1_mat.cols(0, Bk-2) = varepsilontilde_tminusBkminus1TOtminus2_mat.cols(1, Bk-1);
    varepsilontilde_tminusBkTOtminus1_mat.col(Bk-1) = varepsilontilde_tminus1;
    arma::mat Psitilde_tminus1 = fc_Psi(m, varepsilontilde_tminusBkTOtminus1_mat); // \widetilde{\Psi}_{t-1}(\bm\delta)
    arma::mat Rtilde_t = fc_Rtilde_t(m, Bk, beta_1, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1); // \widetilde{R}_{t}(\bm\theta)
    // \widetilde{H}_2(\bm\theta), ..., \widetilde{H}_n(\bm\theta)
    arma::mat Htilde_t = fc_Htilde_t(Dtilde_t, Rtilde_t); // \widetilde{H}_{t}(\bm\theta)
    // \widetilde{\ell}_2(\bm\theta), ..., \widetilde{\ell}_n(\bm\theta)
    double ltilde_t = fc_ltilde_t(y_t, Htilde_t); // \widetilde{\ell}_{t}(\bm\theta)
    sum_ltilde = sum_ltilde + ltilde_t;
    // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
    varepsilontilde_tminusBkminus1TOtminus2_mat = varepsilontilde_tminusBkTOtminus1_mat;
    varepsilontilde_tminus1 = varepsilontilde_t;
    // \widetilde{R}_{t-1}(\bm\theta)
    Rtilde_tminus1 = Rtilde_t;
  }
  double mathcalLtilde_n = sum_ltilde / (n * 1.0); // \widetilde\mathcal{L}_n(\bm\theta)
  return mathcalLtilde_n;
}

////////////////////////////////////////////////// gradient //////////////////////////////////////////////////

arma::mat fc_lowrank_dlnhulinetilde_ddelta_1_mat(int m, int dimdelta) {
  // the \ell-th column is \partial\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\delta_\ell
  arma::mat dlnhulinetilde_ddelta_1_mat(m, dimdelta); dlnhulinetilde_ddelta_1_mat.fill(0.0);
  arma::mat I_m(m,m); I_m.eye(m,m);
  dlnhulinetilde_ddelta_1_mat.cols(0, m-1) = I_m; // \partial\ln\widetilde\underline{\bm{h}}_1(\bm\delta)/\partial\underline\omega_\ell = \bm{e}_\ell
  return dlnhulinetilde_ddelta_1_mat;
}

arma::mat fc_lowrank_dlnhulinetilde_ddelta_t_mat(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::mat y_m_n, arma::cube sum_lambdaklnyuline_0_m_nminus1_r, arma::cube sum_lambdaklnyuline_1_m_nminus1_r, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s, arma::cube sum_gammaklnyuline_11_m_nminus1_s, arma::cube sum_gammaklnyuline_12_m_nminus1_s, arma::cube sum_gammaklnyuline_21_m_nminus1_s, arma::cube sum_gammaklnyuline_22_m_nminus1_s) {
  // r>0 and s>0, the \ell-th column is \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
  arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, 2*m, r);
  arma::vec g_1_vec = kappa_vec.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, 4*m, s);
  arma::vec g_2_vec = kappa_vec.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_mat = fc_asmat(g_2_vec, 4*m, s);
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat zerosmat_m_1(m, 1); zerosmat_m_1.fill(0.0);
  int Dim1 = m+r+2*s; int Dim2 = Dim1 + 2*r*m; int Dim3 = Dim2 + 4*s*m;
  arma::vec y_tminus1 = y_m_n.col(t-2); // \bm{y}_{t-1}
  arma::vec lnyuline_tminus1 = log(pow(y_tminus1,2)); // \ln\underline{\bm{y}}_{t-1} = (\ln{y_{1,t-1}^2}, ..., \ln{y_{m,t-1}^2})'
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta
  arma::mat dlnhulinetilde_ddelta_t_mat(m, dimdelta); dlnhulinetilde_ddelta_t_mat.fill(0.0);
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\underline\omega
  dlnhulinetilde_ddelta_t_mat.cols(0, m-1) = I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta)/\partial\underline\omega_\ell = \bm{e}_\ell
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\lambda, \partial\bm{g}_{0}'
  arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  arma::mat sum_lambdaklnyuline_1_tminus1_m_r = sum_lambdaklnyuline_1_m_nminus1_r.col(t-2);
  for(int k = 0; k < r; k++) {
    // double lambda_k = lambda_vec(k);
    arma::vec g_0k = g_0_mat.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2); 
    arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
    arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
    arma::mat G_0k = g_0k1 * g_0k2.t();
    // arma::vec sum_lambdaklnyuline_0 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_lambdaklnyuline_1 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    // for(int i = 2; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0 + pow(lambda_k, i-1) * lnyuline_tminusi;
    //   sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1 + (i-1) * pow(lambda_k, i-2) * lnyuline_tminusi;
    // }
    arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1_tminus1_m_r.col(k); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    dlnhulinetilde_ddelta_t_mat.col(m+k) = G_0k * sum_lambdaklnyuline_1; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k
    arma::mat sum_lambdaklnyuline_0_mat = fc_asmat(sum_lambdaklnyuline_0, m, 1);
    dlnhulinetilde_ddelta_t_mat.cols(Dim1 +k*2*m,   Dim1 +k*2*m+m-1)   = dot(g_0k2, sum_lambdaklnyuline_0) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{0,k,1}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim1 +k*2*m+m, Dim1 +(k+1)*2*m-1) = g_0k1 * sum_lambdaklnyuline_0_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{0,k,2}'
  }
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\gamma, \partial\bm\varphi, \partial\bm{g}_{1}', \partial\bm{g}_{2}'
  arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_11_tminus1_m_s = sum_gammaklnyuline_11_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_12_tminus1_m_s = sum_gammaklnyuline_12_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_21_tminus1_m_s = sum_gammaklnyuline_21_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_22_tminus1_m_s = sum_gammaklnyuline_22_m_nminus1_s.col(t-2);
  for(int k = 0; k < s; k++) {
    // double gamma_k = gamma_vec(k);
    // double varphi_k = varphi_vec(k);
    arma::vec g_1k = g_1_mat.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4); 
    arma::vec g_2k = g_2_mat.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4); 
    arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
    arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
    arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
    arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
    arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
    arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
    arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
    arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
    arma::mat G_1k = g_1k1 * g_1k2.t() + g_1k3 * g_1k4.t();
    arma::mat G_2k = g_2k1 * g_2k2.t() + g_2k3 * g_2k4.t();
    // arma::vec sum_gammaklnyuline_01 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_02 = zerosmat_m_1.col(0); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_11 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_12 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_21 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_22 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // for(int i = 2; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_gammaklnyuline_01 = sum_gammaklnyuline_01 + pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_02 = sum_gammaklnyuline_02 + pow(gamma_k, i-1) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_11 = sum_gammaklnyuline_11 + (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_12 = sum_gammaklnyuline_12 + (i-1) * pow(gamma_k, i-2) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_21 = sum_gammaklnyuline_21 + (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_22 = sum_gammaklnyuline_22 + (i-1) * pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    // }
    arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_11 = sum_gammaklnyuline_11_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_12 = sum_gammaklnyuline_12_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_21 = sum_gammaklnyuline_21_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_22 = sum_gammaklnyuline_22_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    dlnhulinetilde_ddelta_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_11 + G_2k * sum_gammaklnyuline_12; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k
    dlnhulinetilde_ddelta_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_21 + G_2k * sum_gammaklnyuline_22; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k
    arma::mat sum_gammaklnyuline_01_mat = fc_asmat(sum_gammaklnyuline_01, m, 1);
    arma::mat sum_gammaklnyuline_02_mat = fc_asmat(sum_gammaklnyuline_02, m, 1);
    dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m,     Dim2 +k*4*m+m-1)   = dot(g_1k2, sum_gammaklnyuline_01) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,1}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m+m,   Dim2 +k*4*m+2*m-1) = g_1k1 * sum_gammaklnyuline_01_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,2}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m+2*m, Dim2 +k*4*m+3*m-1) = dot(g_1k4, sum_gammaklnyuline_01) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,3}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m+3*m, Dim2 +(k+1)*4*m-1) = g_1k3 * sum_gammaklnyuline_01_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,4}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m,     Dim3 +k*4*m+m-1)   = dot(g_2k2, sum_gammaklnyuline_02) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,1}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m+m,   Dim3 +k*4*m+2*m-1) = g_2k1 * sum_gammaklnyuline_02_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,2}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m+2*m, Dim3 +k*4*m+3*m-1) = dot(g_2k4, sum_gammaklnyuline_02) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,3}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m+3*m, Dim3 +(k+1)*4*m-1) = g_2k3 * sum_gammaklnyuline_02_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,4}'
  }
  return dlnhulinetilde_ddelta_t_mat;
}

arma::mat fc_lowrank_dlnhulinetilde_ddelta_t_mat_requal0(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::mat y_m_n, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s, arma::cube sum_gammaklnyuline_11_m_nminus1_s, arma::cube sum_gammaklnyuline_12_m_nminus1_s, arma::cube sum_gammaklnyuline_21_m_nminus1_s, arma::cube sum_gammaklnyuline_22_m_nminus1_s) {
  // r=0 and s>0, the \ell-th column is \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
  // arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  // arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, 2*m, r);
  arma::vec g_1_vec = kappa_vec.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, 4*m, s);
  arma::vec g_2_vec = kappa_vec.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_mat = fc_asmat(g_2_vec, 4*m, s);
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat zerosmat_m_1(m, 1); zerosmat_m_1.fill(0.0);
  int Dim1 = m+r+2*s; int Dim2 = Dim1 + 2*r*m; int Dim3 = Dim2 + 4*s*m;
  arma::vec y_tminus1 = y_m_n.col(t-2); // \bm{y}_{t-1}
  arma::vec lnyuline_tminus1 = log(pow(y_tminus1,2)); // \ln\underline{\bm{y}}_{t-1} = (\ln{y_{1,t-1}^2}, ..., \ln{y_{m,t-1}^2})'
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta
  arma::mat dlnhulinetilde_ddelta_t_mat(m, dimdelta); dlnhulinetilde_ddelta_t_mat.fill(0.0);
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\underline\omega
  dlnhulinetilde_ddelta_t_mat.cols(0, m-1) = I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta)/\partial\underline\omega_\ell = \bm{e}_\ell
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\lambda, \partial\bm{g}_{0}'
  // arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  // arma::mat sum_lambdaklnyuline_1_tminus1_m_r = sum_lambdaklnyuline_1_m_nminus1_r.col(t-2);
  // for(int k = 0; k < r; k++) {
  //   // double lambda_k = lambda_vec(k);
  //   arma::vec g_0k = g_0_mat.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2);
  //   arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
  //   arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
  //   arma::mat G_0k = g_0k1 * g_0k2.t();
  //   // arma::vec sum_lambdaklnyuline_0 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_lambdaklnyuline_1 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
  //   // for(int i = 2; i < t; i++) {
  //   //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
  //   //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
  //   //   sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0 + pow(lambda_k, i-1) * lnyuline_tminusi;
  //   //   sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1 + (i-1) * pow(lambda_k, i-2) * lnyuline_tminusi;
  //   // }
  //   arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1_tminus1_m_r.col(k); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
  //   dlnhulinetilde_ddelta_t_mat.col(m+k) = G_0k * sum_lambdaklnyuline_1; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k
  //   arma::mat sum_lambdaklnyuline_0_mat = fc_asmat(sum_lambdaklnyuline_0, m, 1);
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim1 +k*2*m,   Dim1 +k*2*m+m-1)   = dot(g_0k2, sum_lambdaklnyuline_0) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{0,k,1}'
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim1 +k*2*m+m, Dim1 +(k+1)*2*m-1) = g_0k1 * sum_lambdaklnyuline_0_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{0,k,2}'
  // }
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\gamma, \partial\bm\varphi, \partial\bm{g}_{1}', \partial\bm{g}_{2}'
  arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_11_tminus1_m_s = sum_gammaklnyuline_11_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_12_tminus1_m_s = sum_gammaklnyuline_12_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_21_tminus1_m_s = sum_gammaklnyuline_21_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_22_tminus1_m_s = sum_gammaklnyuline_22_m_nminus1_s.col(t-2);
  for(int k = 0; k < s; k++) {
    // double gamma_k = gamma_vec(k);
    // double varphi_k = varphi_vec(k);
    arma::vec g_1k = g_1_mat.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4); 
    arma::vec g_2k = g_2_mat.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4); 
    arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
    arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
    arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
    arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
    arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
    arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
    arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
    arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
    arma::mat G_1k = g_1k1 * g_1k2.t() + g_1k3 * g_1k4.t();
    arma::mat G_2k = g_2k1 * g_2k2.t() + g_2k3 * g_2k4.t();
    // arma::vec sum_gammaklnyuline_01 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_02 = zerosmat_m_1.col(0); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_11 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_12 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_21 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_22 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // for(int i = 2; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_gammaklnyuline_01 = sum_gammaklnyuline_01 + pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_02 = sum_gammaklnyuline_02 + pow(gamma_k, i-1) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_11 = sum_gammaklnyuline_11 + (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_12 = sum_gammaklnyuline_12 + (i-1) * pow(gamma_k, i-2) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_21 = sum_gammaklnyuline_21 + (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_22 = sum_gammaklnyuline_22 + (i-1) * pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    // }
    arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_11 = sum_gammaklnyuline_11_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_12 = sum_gammaklnyuline_12_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_21 = sum_gammaklnyuline_21_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_22 = sum_gammaklnyuline_22_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    dlnhulinetilde_ddelta_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_11 + G_2k * sum_gammaklnyuline_12; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k
    dlnhulinetilde_ddelta_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_21 + G_2k * sum_gammaklnyuline_22; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k
    arma::mat sum_gammaklnyuline_01_mat = fc_asmat(sum_gammaklnyuline_01, m, 1);
    arma::mat sum_gammaklnyuline_02_mat = fc_asmat(sum_gammaklnyuline_02, m, 1);
    dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m,     Dim2 +k*4*m+m-1)   = dot(g_1k2, sum_gammaklnyuline_01) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,1}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m+m,   Dim2 +k*4*m+2*m-1) = g_1k1 * sum_gammaklnyuline_01_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,2}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m+2*m, Dim2 +k*4*m+3*m-1) = dot(g_1k4, sum_gammaklnyuline_01) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,3}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m+3*m, Dim2 +(k+1)*4*m-1) = g_1k3 * sum_gammaklnyuline_01_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,4}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m,     Dim3 +k*4*m+m-1)   = dot(g_2k2, sum_gammaklnyuline_02) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,1}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m+m,   Dim3 +k*4*m+2*m-1) = g_2k1 * sum_gammaklnyuline_02_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,2}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m+2*m, Dim3 +k*4*m+3*m-1) = dot(g_2k4, sum_gammaklnyuline_02) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,3}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m+3*m, Dim3 +(k+1)*4*m-1) = g_2k3 * sum_gammaklnyuline_02_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,4}'
  }
  return dlnhulinetilde_ddelta_t_mat;
}

arma::mat fc_lowrank_dlnhulinetilde_ddelta_t_mat_sequal0(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::mat y_m_n, arma::cube sum_lambdaklnyuline_0_m_nminus1_r, arma::cube sum_lambdaklnyuline_1_m_nminus1_r) {
  // r>0 and s=0, the \ell-th column is \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
  arma::vec lambda_vec = kappa_vec.subvec(0, r-1); 
  // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, 2*m, r);
  // arma::vec g_1_vec = kappa_vec.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, 4*m, s);
  // arma::vec g_2_vec = kappa_vec.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_mat = fc_asmat(g_2_vec, 4*m, s);
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat zerosmat_m_1(m, 1); zerosmat_m_1.fill(0.0);
  int Dim1 = m+r+2*s; // int Dim2 = Dim1 + 2*r*m; int Dim3 = Dim2 + 4*s*m;
  arma::vec y_tminus1 = y_m_n.col(t-2); // \bm{y}_{t-1}
  arma::vec lnyuline_tminus1 = log(pow(y_tminus1,2)); // \ln\underline{\bm{y}}_{t-1} = (\ln{y_{1,t-1}^2}, ..., \ln{y_{m,t-1}^2})'
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta
  arma::mat dlnhulinetilde_ddelta_t_mat(m, dimdelta); dlnhulinetilde_ddelta_t_mat.fill(0.0);
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\underline\omega
  dlnhulinetilde_ddelta_t_mat.cols(0, m-1) = I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta)/\partial\underline\omega_\ell = \bm{e}_\ell
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\lambda, \partial\bm{g}_{0}'
  arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  arma::mat sum_lambdaklnyuline_1_tminus1_m_r = sum_lambdaklnyuline_1_m_nminus1_r.col(t-2);
  for(int k = 0; k < r; k++) {
    // double lambda_k = lambda_vec(k);
    arma::vec g_0k = g_0_mat.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2); 
    arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
    arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
    arma::mat G_0k = g_0k1 * g_0k2.t();
    // arma::vec sum_lambdaklnyuline_0 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_lambdaklnyuline_1 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    // for(int i = 2; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0 + pow(lambda_k, i-1) * lnyuline_tminusi;
    //   sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1 + (i-1) * pow(lambda_k, i-2) * lnyuline_tminusi;
    // }
    arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1_tminus1_m_r.col(k); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    dlnhulinetilde_ddelta_t_mat.col(m+k) = G_0k * sum_lambdaklnyuline_1; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k
    arma::mat sum_lambdaklnyuline_0_mat = fc_asmat(sum_lambdaklnyuline_0, m, 1);
    dlnhulinetilde_ddelta_t_mat.cols(Dim1 +k*2*m,   Dim1 +k*2*m+m-1)   = dot(g_0k2, sum_lambdaklnyuline_0) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{0,k,1}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim1 +k*2*m+m, Dim1 +(k+1)*2*m-1) = g_0k1 * sum_lambdaklnyuline_0_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{0,k,2}'
  }
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\gamma, \partial\bm\varphi, \partial\bm{g}_{1}', \partial\bm{g}_{2}'
  // arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_11_tminus1_m_s = sum_gammaklnyuline_11_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_12_tminus1_m_s = sum_gammaklnyuline_12_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_21_tminus1_m_s = sum_gammaklnyuline_21_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_22_tminus1_m_s = sum_gammaklnyuline_22_m_nminus1_s.col(t-2);
  // for(int k = 0; k < s; k++) {
  //   // double gamma_k = gamma_vec(k);
  //   // double varphi_k = varphi_vec(k);
  //   arma::vec g_1k = g_1_mat.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4);
  //   arma::vec g_2k = g_2_mat.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4);
  //   arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
  //   arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
  //   arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
  //   arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
  //   arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
  //   arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
  //   arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
  //   arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
  //   arma::mat G_1k = g_1k1 * g_1k2.t() + g_1k3 * g_1k4.t();
  //   arma::mat G_2k = g_2k1 * g_2k2.t() + g_2k3 * g_2k4.t();
  //   // arma::vec sum_gammaklnyuline_01 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_02 = zerosmat_m_1.col(0); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_11 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_12 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_21 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_22 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // for(int i = 2; i < t; i++) {
  //   //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
  //   //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
  //   //   sum_gammaklnyuline_01 = sum_gammaklnyuline_01 + pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_02 = sum_gammaklnyuline_02 + pow(gamma_k, i-1) * sin((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_11 = sum_gammaklnyuline_11 + (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_12 = sum_gammaklnyuline_12 + (i-1) * pow(gamma_k, i-2) * sin((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_21 = sum_gammaklnyuline_21 + (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_22 = sum_gammaklnyuline_22 + (i-1) * pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
  //   // }
  //   arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_11 = sum_gammaklnyuline_11_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_12 = sum_gammaklnyuline_12_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_21 = sum_gammaklnyuline_21_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_22 = sum_gammaklnyuline_22_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   dlnhulinetilde_ddelta_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_11 + G_2k * sum_gammaklnyuline_12; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k
  //   dlnhulinetilde_ddelta_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_21 + G_2k * sum_gammaklnyuline_22; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k
  //   arma::mat sum_gammaklnyuline_01_mat = fc_asmat(sum_gammaklnyuline_01, m, 1);
  //   arma::mat sum_gammaklnyuline_02_mat = fc_asmat(sum_gammaklnyuline_02, m, 1);
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m,     Dim2 +k*4*m+m-1)   = dot(g_1k2, sum_gammaklnyuline_01) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,1}'
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m+m,   Dim2 +k*4*m+2*m-1) = g_1k1 * sum_gammaklnyuline_01_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,2}'
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m+2*m, Dim2 +k*4*m+3*m-1) = dot(g_1k4, sum_gammaklnyuline_01) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,3}'
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*4*m+3*m, Dim2 +(k+1)*4*m-1) = g_1k3 * sum_gammaklnyuline_01_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k,4}'
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m,     Dim3 +k*4*m+m-1)   = dot(g_2k2, sum_gammaklnyuline_02) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,1}'
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m+m,   Dim3 +k*4*m+2*m-1) = g_2k1 * sum_gammaklnyuline_02_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,2}'
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m+2*m, Dim3 +k*4*m+3*m-1) = dot(g_2k4, sum_gammaklnyuline_02) * I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,3}'
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*4*m+3*m, Dim3 +(k+1)*4*m-1) = g_2k3 * sum_gammaklnyuline_02_mat.t(); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k,4}'
  // }
  return dlnhulinetilde_ddelta_t_mat;
}

arma::mat fc_lowrank_dlnhulinetilde_ddelta_t_mat_general(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::mat y_m_n, arma::cube sum_lambdaklnyuline_0_m_nminus1_r, arma::cube sum_lambdaklnyuline_1_m_nminus1_r, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s, arma::cube sum_gammaklnyuline_11_m_nminus1_s, arma::cube sum_gammaklnyuline_12_m_nminus1_s, arma::cube sum_gammaklnyuline_21_m_nminus1_s, arma::cube sum_gammaklnyuline_22_m_nminus1_s) {
  // the \ell-th column is \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
  arma::mat dlnhulinetilde_ddelta_t_mat(m, dimdelta); dlnhulinetilde_ddelta_t_mat.fill(0.0);
  if(s == 0) {
    dlnhulinetilde_ddelta_t_mat = fc_lowrank_dlnhulinetilde_ddelta_t_mat_sequal0(t, m, r, s, dimdelta, kappa_vec, y_m_n, sum_lambdaklnyuline_0_m_nminus1_r, sum_lambdaklnyuline_1_m_nminus1_r);
  } else if(r == 0) {
    dlnhulinetilde_ddelta_t_mat = fc_lowrank_dlnhulinetilde_ddelta_t_mat_requal0(t, m, r, s, dimdelta, kappa_vec, y_m_n, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s);
  } else {
    dlnhulinetilde_ddelta_t_mat = fc_lowrank_dlnhulinetilde_ddelta_t_mat(t, m, r, s, dimdelta, kappa_vec, y_m_n, sum_lambdaklnyuline_0_m_nminus1_r, sum_lambdaklnyuline_1_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s);
  }
  return dlnhulinetilde_ddelta_t_mat;
}

arma::mat fc_dDtilde_ddelta_t_l(arma::mat Dtilde_t, arma::vec dlnhulinetilde_ddelta_t_l) {
  // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
  arma::vec dDtilde_ddelta_t_l_diag = 0.5 * Dtilde_t.diag(0) % dlnhulinetilde_ddelta_t_l; // \partial\widetilde{D}_t(\bm\delta)/\partial\delta_\ell = 1/2 \widetilde{D}_t(\bm\delta) \Diag{\partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta)/\partial\delta_\ell}
  arma::mat dDtilde_ddelta_t_l = diagmat(dDtilde_ddelta_t_l_diag, 0);
  return dDtilde_ddelta_t_l;
}

arma::mat fc_dPsitilde_ddelta_tminus1_l(int m, arma::mat varepsilontilde_tminusBkTOtminus1_mat, arma::mat Psitilde_tminus1, arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat) {
  // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
  arma::mat dPsitilde_ddelta_tminus1_l(m, m); dPsitilde_ddelta_tminus1_l.fill(0.0);
  for(int i = 1; i < m; i++) { // {2,1}; {3,1}, {3,2}; ...; {m,1}, ..., {m,m-1}
    for(int j = 0; j < i; j++) {
      double Psitilde_tminus1_ij = Psitilde_tminus1(i,j);
      arma::rowvec varepsilon_ith_Bk = varepsilontilde_tminusBkTOtminus1_mat.row(i); // (\widetilde\varepsilon_{i,t-\Bbbk}(\bm\delta), ..., \widetilde\varepsilon_{i,t-1}(\bm\delta))'
      arma::rowvec varepsilon_jth_Bk = varepsilontilde_tminusBkTOtminus1_mat.row(j); 
      arma::rowvec dlnhuline_ith_Bk = dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat.row(i); // (\patial\ln\widetilde{h}_{ii,t-\Bbbk}/\partial\delta_\ell, ..., \patial\ln\widetilde{h}_{ii,t-1}/\partial\delta_\ell)'
      arma::rowvec dlnhuline_jth_Bk = dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat.row(j);
      arma::rowvec dlnhuline_varepsilon_i_Bk = dlnhuline_ith_Bk % varepsilon_ith_Bk; // (\patial\ln\widetilde{h}_{ii,t-\Bbbk}/\partial\delta_\ell * \widetilde\varepsilon_{i,t-\Bbbk}(\bm\delta), ..., \patial\ln\widetilde{h}_{ii,t-1}/\partial\delta_\ell) * \widetilde\varepsilon_{i,t-1}(\bm\delta))'
      arma::rowvec dlnhuline_varepsilon_j_Bk = dlnhuline_jth_Bk % varepsilon_jth_Bk; // (\patial\ln\widetilde{h}_{ii,t-\Bbbk}/\partial\delta_\ell * \widetilde\varepsilon_{i,t-\Bbbk}(\bm\delta), ..., \patial\ln\widetilde{h}_{ii,t-1}/\partial\delta_\ell) * \widetilde\varepsilon_{i,t-1}(\bm\delta))'
      double denominator_i = dot(varepsilon_ith_Bk, varepsilon_ith_Bk); 
      double denominator_j = dot(varepsilon_jth_Bk, varepsilon_jth_Bk);
      double denominator_ij = sqrt(denominator_i * denominator_j);
      double numerator_i = dot(dlnhuline_varepsilon_i_Bk, varepsilon_ith_Bk);
      double numerator_j = dot(dlnhuline_varepsilon_j_Bk, varepsilon_jth_Bk);
      double numerator_ij = dot(dlnhuline_varepsilon_i_Bk, varepsilon_jth_Bk) + dot(dlnhuline_varepsilon_j_Bk, varepsilon_ith_Bk);
      dPsitilde_ddelta_tminus1_l(i,j) = dPsitilde_ddelta_tminus1_l(j,i) = -0.5 * numerator_ij / denominator_ij + 0.5 * Psitilde_tminus1_ij * (numerator_i / denominator_i + numerator_j / denominator_j);
    }
  }
  return dPsitilde_ddelta_tminus1_l;
}

arma::cube fc_dRtilde_ddelta_1_cube(int m, int dimdelta) {
  // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\delta_\ell = 0_m
  arma::cube dRtilde_ddelta_1_cube(m, m, dimdelta); dRtilde_ddelta_1_cube.fill(0.0); 
  return dRtilde_ddelta_1_cube; 
}

arma::cube fc_dRtilde_ddelta_t_cube(int m, int dimdelta, double beta_1, double beta_2, arma::mat varepsilontilde_tminusBkTOtminus1_mat, arma::mat Psitilde_tminus1, arma::cube dlnhulinetilde_ddelta_tminusBkTOtminus1_cube, arma::cube dRtilde_ddelta_tminus1_cube) {
  // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
  arma::cube dRtilde_ddelta_t_cube(m, m, dimdelta); dRtilde_ddelta_t_cube.fill(0.0);
  for(int l = 0; l < dimdelta; l++) {
    arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.col(l); // the k-th column is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-k)}(\bm\delta) / \partial\delta_\ell
    arma::mat dPsitilde_ddelta_tminus1_l = fc_dPsitilde_ddelta_tminus1_l(m, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
    arma::mat dRtilde_ddelta_tminus1_l = dRtilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
    dRtilde_ddelta_t_cube.slice(l) = beta_1 * dPsitilde_ddelta_tminus1_l + beta_2 * dRtilde_ddelta_tminus1_l;
  }
  return dRtilde_ddelta_t_cube;
}

arma::cube fc_dRuline_dbeta_cube(int m, int dimbeta) {
  // the \ell-th slice is \partial\underline{R} / \partial\beta_\ell
  arma::cube dRuline_dbeta_cube(m, m, dimbeta); dRuline_dbeta_cube.fill(0.0);
  arma::mat zerosmat_m_m(m, m); zerosmat_m_m.fill(0.0);
  dRuline_dbeta_cube.slice(0) = zerosmat_m_m;
  dRuline_dbeta_cube.slice(1) = zerosmat_m_m;
  int l = 2;
  for(int j = 0; j < (m-1); j++) { // \partial\underline{R}_{2,1}, ..., \partial\underline{R}_{m,1}; \partial\underline{R}_{3,2}, ..., \partial\underline{R}_{m,2}; ...; \partial\underline{R}_{m,m-1}
    for(int i = (j+1); i < m; i++) {
      dRuline_dbeta_cube(i,j,l) = dRuline_dbeta_cube(j,i,l) = 1.0;
      l = l + 1;
    }
  }
  return dRuline_dbeta_cube;
}

arma::cube fc_dRtilde_dbeta_1_cube(int m, int dimbeta, double beta_1, double beta_2, arma::mat Ruline, arma::cube dRuline_dbeta_cube) {
  // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\beta_\ell
  arma::cube dRtilde_dbeta_1_cube(m, m, dimbeta); dRtilde_dbeta_1_cube.fill(0.0); 
  arma::mat onesmat_m(m,m); onesmat_m.fill(1.0);
  dRtilde_dbeta_1_cube.slice(0) = - 1.0 / (1.0 - beta_2) * Ruline + 1.0 / (1.0 - beta_2) * onesmat_m; 
  dRtilde_dbeta_1_cube.slice(1) = - beta_1 / pow(1.0-beta_2, 2) * Ruline + beta_1 / pow(1.0-beta_2, 2) * onesmat_m; 
  dRtilde_dbeta_1_cube.slices(2, dimbeta-1) = (1.0 - beta_1 / (1.0 - beta_2)) * dRuline_dbeta_cube.slices(2, dimbeta-1); 
  return dRtilde_dbeta_1_cube;
}

arma::cube fc_dRtilde_dbeta_t_cube(int m, int dimbeta, double beta_2, arma::mat Ruline, arma::mat Psitilde_tminus1, arma::mat Rtilde_tminus1, arma::cube dRtilde_dbeta_tminus1_cube) {
  // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\beta_\ell
  arma::cube dRtilde_dbeta_t_cube(m, m, dimbeta); dRtilde_dbeta_t_cube.fill(0.0);
  dRtilde_dbeta_t_cube.slice(0) = - Ruline + Psitilde_tminus1 + beta_2 * dRtilde_dbeta_tminus1_cube.slice(0);
  dRtilde_dbeta_t_cube.slice(1) = - Ruline + Rtilde_tminus1 + beta_2 * dRtilde_dbeta_tminus1_cube.slice(1);
  dRtilde_dbeta_t_cube.slices(2, dimbeta-1) = dRtilde_dbeta_tminus1_cube.slices(2, dimbeta-1);
  return dRtilde_dbeta_t_cube;
}

arma::mat fc_dHtilde_ddelta_t_l(arma::mat Dtilde_t, arma::mat Rtilde_t, arma::mat dDtilde_ddelta_t_l, arma::mat dRtilde_ddelta_t_l) {
  // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
  arma::mat dHtilde_ddelta_t_l = dDtilde_ddelta_t_l * Rtilde_t * Dtilde_t + 
                                 Dtilde_t * dRtilde_ddelta_t_l * Dtilde_t + 
                                 Dtilde_t * Rtilde_t * dDtilde_ddelta_t_l;
  return dHtilde_ddelta_t_l; 
}

arma::mat fc_dHtilde_dbeta_t_l(arma::mat Dtilde_t, arma::mat dRtilde_dbeta_t_l) {
  // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
  arma::mat dHtilde_dbeta_t_l = Dtilde_t * dRtilde_dbeta_t_l * Dtilde_t; 
  return dHtilde_dbeta_t_l;
}

arma::vec fc_dltilde_dvartheta_t(int m, int dimdelta, int dimbeta, int dimvartheta, arma::vec y_t, arma::mat Htilde_t, arma::mat Dtilde_t, arma::mat Rtilde_t, arma::mat dlnhulinetilde_ddelta_t_mat, arma::cube dRtilde_ddelta_t_cube, arma::cube dRtilde_dbeta_t_cube) {
  // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\vartheta
  arma::vec dltilde_dvartheta_t(dimvartheta); // \partial\widetilde{\ell}_t(\bm\theta) / \partial\vartheta_\ell = \tr[(I_m - \widetilde{H}_t^{-1}(\bm\theta) \bm{y}_t \bm{y}_t') \widetilde{H}_t^{-1}(\bm\theta) \partial\widetilde{H}_t(\bm\theta)/\partial\vartheta_\ell]
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat y_t_mat = fc_asmat(y_t, m, 1);
  arma::mat inverse_Htilde_t = Htilde_t.i();
  arma::mat left_mat = (I_m - inverse_Htilde_t * y_t * y_t_mat.t()) * inverse_Htilde_t;
  for(int l = 0; l < dimdelta; l++) {
    arma::vec dlnhulinetilde_ddelta_t_l = dlnhulinetilde_ddelta_t_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
    arma::mat dDtilde_ddelta_t_l = fc_dDtilde_ddelta_t_l(Dtilde_t, dlnhulinetilde_ddelta_t_l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
    arma::mat dRtilde_ddelta_t_l = dRtilde_ddelta_t_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
    arma::mat dHtilde_ddelta_t_l = fc_dHtilde_ddelta_t_l(Dtilde_t, Rtilde_t, dDtilde_ddelta_t_l, dRtilde_ddelta_t_l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
    dltilde_dvartheta_t(l) = trace(left_mat * dHtilde_ddelta_t_l); // \partial\widetilde{\ell}_t(\bm\theta) / \partial\delta_\ell
  }
  for(int l = 0; l < dimbeta; l++) {
    arma::mat dRtilde_dbeta_t_l = dRtilde_dbeta_t_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_\ell
    arma::mat dHtilde_dbeta_t_l = fc_dHtilde_dbeta_t_l(Dtilde_t, dRtilde_dbeta_t_l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
    dltilde_dvartheta_t(dimdelta+l) = trace(left_mat * dHtilde_dbeta_t_l); // \partial\widetilde{\ell}_t(\bm\theta) / \partial\beta_\ell
  }
  return dltilde_dvartheta_t;
}

// [[Rcpp::export]]
arma::vec fc_dmathcalLtilde_dvartheta(int n, int m, int r, int s, int Bk, arma::vec vartheta_vec, arma::mat y_m_n) { 
  // \partial\widetilde{\mathcal{L}}_{n}(\bm\theta) / \partial\bm\vartheta, under initial values \widetilde{y}_s = 1_m for s <= 0
  int dimkappa = r + 2*s + 2*r*m + 8*s*m;
  int dimbeta = 2 + m*(m-1)/2;
  int dimdelta = m + dimkappa;
  int dimvartheta = dimdelta + dimbeta;
  arma::vec omegauline = vartheta_vec.subvec(0, m-1); // \underline{\bm\omega}
  arma::vec kappa_vec = vartheta_vec.subvec(m, m+dimkappa-1); // \bm\kappa
  arma::vec beta_vec = vartheta_vec.subvec(m+dimkappa, m+dimkappa+dimbeta-1); // \bm\beta
  double beta_1 = beta_vec(0); double beta_2 = beta_vec(1); // \beta_1 and \beta_2
  arma::vec ruline = beta_vec.tail(dimbeta-2); // \underline\bm{r}
  arma::mat Ruline = fc_Ruline(m, ruline); // \underline{R}
  // summations using FFT algorithm
  int max_r_1 = r; if(max_r_1 == 0) {max_r_1 = 1;}
  int max_s_1 = s; if(max_s_1 == 0) {max_s_1 = 1;}
  arma::cube sum_lambdaklnyuline_0_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_0_m_nminus1_r.fill(0.0);
  arma::cube sum_lambdaklnyuline_1_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_1_m_nminus1_r.fill(0.0);
  arma::cube sum_gammaklnyuline_01_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_01_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_02_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_02_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_11_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_11_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_12_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_12_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_21_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_21_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_22_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_22_m_nminus1_s.fill(0.0);
  if(s == 0) {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
  } else if(r == 0) {
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  } else {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  }
  // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\vartheta
  arma::mat sum_dltilde_dvartheta_mat(dimvartheta, 1); sum_dltilde_dvartheta_mat.fill(0.0);
  arma::vec sum_dltilde_dvartheta = sum_dltilde_dvartheta_mat.col(0); // \sum_{t=1}^{n} \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\vartheta
  // t=1
  arma::vec y_1 = y_m_n.col(0); // \bm{y}_{1}
  arma::mat Dtilde_1 = fc_Dtilde_1(omegauline); // \widetilde{D}_{1}(\bm\delta)
  arma::mat inverse_Dtilde_1 = fc_inverse_Dtilde_1(omegauline); // \widetilde{D}_{1}^{-1}(\bm\delta)
  arma::mat Rtilde_1 = fc_Rtilde_1(m, beta_1, beta_2, Ruline); // \widetilde{R}_{1}(\bm\theta)
  arma::mat Htilde_1 = fc_Htilde_t(Dtilde_1, Rtilde_1); // \widetilde{H}_{1}(\bm\theta)
  arma::vec varepsilontilde_1 = inverse_Dtilde_1 * y_1; // \widetilde{\varepsilon}_{1}(\bm\delta) = \widetilde{D}_{1}^{-1}(\bm\delta) \bm{y}_{1}
  // derivatives
  arma::mat dlnhulinetilde_ddelta_1_mat = fc_lowrank_dlnhulinetilde_ddelta_1_mat(m, dimdelta); // \partial\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\bm\delta'_\ell, which is equal to \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'_\ell for t <= 0
  arma::cube dRtilde_ddelta_1_cube = fc_dRtilde_ddelta_1_cube(m, dimdelta); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\delta_\ell
  arma::cube dRuline_dbeta_cube = fc_dRuline_dbeta_cube(m, dimbeta); // the \ell-th slice is \partial\underline{R} / \partial\beta_\ell
  arma::cube dRtilde_dbeta_1_cube = fc_dRtilde_dbeta_1_cube(m, dimbeta, beta_1, beta_2, Ruline, dRuline_dbeta_cube); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\beta_\ell
  // \partial\widetilde{\ell}_1(\bm\theta)/\partial\bm\vartheta
  arma::vec dltilde_dvartheta_1 = fc_dltilde_dvartheta_t(m, dimdelta, dimbeta, dimvartheta, y_1, Htilde_1, Dtilde_1, Rtilde_1, dlnhulinetilde_ddelta_1_mat, dRtilde_ddelta_1_cube, dRtilde_dbeta_1_cube); // \partial\widetilde{\ell}_{1}(\bm\theta) / \partial\bm\vartheta
  sum_dltilde_dvartheta = sum_dltilde_dvartheta + dltilde_dvartheta_1;
  // t=2,...,n
  // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
  arma::mat varepsilontilde_tminusBkminus1TOtminus2_mat(m, Bk); varepsilontilde_tminusBkminus1TOtminus2_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta))
  varepsilontilde_tminusBkminus1TOtminus2_mat.each_col() = exp(-0.5*omegauline); // for s <= 0, \widetilde{\varepsilon}_s(\bm\delta) = (\exp{-1/2*\underliner{\omega}_1}, ..., \exp{-1/2*\underliner{\omega}_m})'
  arma::vec varepsilontilde_tminus1 = varepsilontilde_1; // \widetilde{\varepsilon}_{t-1}(\bm\delta)
  // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
  arma::cube dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.fill(0.0); // the i-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-2+i}(\bm\delta) / \partial\bm\delta'_\ell
  dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.each_slice() = dlnhulinetilde_ddelta_1_mat;
  arma::mat dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_1_mat;
  // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell and \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell
  arma::mat Rtilde_tminus1 = Rtilde_1; // \widetilde{R}_{t-1}(\bm\theta)
  arma::cube dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
  arma::cube dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_\ell
  for(int t = 2; t < (n+1); t++) {
    arma::vec y_t = y_m_n.col(t-1); // \bm{y}_{t}
    arma::vec lnhulinetilde_t = fc_lowrank_lnhulinetilde_t_general(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s); // \ln\widetilde\underline{\bm{h}}_t(\bm\delta)
    arma::mat Dtilde_t = fc_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}(\bm\delta)
    arma::mat inverse_Dtilde_t = fc_inverse_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}^{-1}(\bm\delta)
    arma::vec varepsilontilde_t = inverse_Dtilde_t * y_t; // \widetilde{\varepsilon}_{t}(\bm\delta) = \widetilde{D}_{t}^{-1}(\bm\delta) \bm{y}_{t}
    arma::mat varepsilontilde_tminusBkTOtminus1_mat(m, Bk); varepsilontilde_tminusBkTOtminus1_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk}(\bm\delta), ..., \widetilde{\varepsilon}_{t-1}(\bm\delta))
    varepsilontilde_tminusBkTOtminus1_mat.cols(0, Bk-2) = varepsilontilde_tminusBkminus1TOtminus2_mat.cols(1, Bk-1);
    varepsilontilde_tminusBkTOtminus1_mat.col(Bk-1) = varepsilontilde_tminus1;
    arma::mat Psitilde_tminus1 = fc_Psi(m, varepsilontilde_tminusBkTOtminus1_mat); // \widetilde{\Psi}_{t-1}(\bm\delta)
    arma::mat Rtilde_t = fc_Rtilde_t(m, Bk, beta_1, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1); // \widetilde{R}_{t}(\bm\theta)
    arma::mat Htilde_t = fc_Htilde_t(Dtilde_t, Rtilde_t); // \widetilde{H}_{t}(\bm\theta)
    // derivatives
    arma::mat dlnhulinetilde_ddelta_t_mat = fc_lowrank_dlnhulinetilde_ddelta_t_mat_general(t, m, r, s, dimdelta, kappa_vec, y_m_n, sum_lambdaklnyuline_0_m_nminus1_r, sum_lambdaklnyuline_1_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'
    arma::cube dlnhulinetilde_ddelta_tminusBkTOtminus1_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.fill(0.0); // the k-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-k)}(\bm\delta) / \partial\bm\delta'
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slices(0, Bk-2) = dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.slices(1, Bk-1);
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slice(Bk-1) = dlnhulinetilde_ddelta_tminus1_mat;
    arma::cube dRtilde_ddelta_t_cube = fc_dRtilde_ddelta_t_cube(m, dimdelta, beta_1, beta_2, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dlnhulinetilde_ddelta_tminusBkTOtminus1_cube, dRtilde_ddelta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
    arma::cube dRtilde_dbeta_t_cube = fc_dRtilde_dbeta_t_cube(m, dimbeta, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1, dRtilde_dbeta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\beta_\ell
    // \partial\widetilde{\ell}_2(\bm\theta)/\partial\bm\vartheta, ..., \partial\widetilde{\ell}_n(\bm\theta)/\partial\bm\vartheta
    arma::vec dltilde_dvartheta_t = fc_dltilde_dvartheta_t(m, dimdelta, dimbeta, dimvartheta, y_t, Htilde_t, Dtilde_t, Rtilde_t, dlnhulinetilde_ddelta_t_mat, dRtilde_ddelta_t_cube, dRtilde_dbeta_t_cube); // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\vartheta
    sum_dltilde_dvartheta = sum_dltilde_dvartheta + dltilde_dvartheta_t;
    // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
    varepsilontilde_tminusBkminus1TOtminus2_mat = varepsilontilde_tminusBkTOtminus1_mat;
    varepsilontilde_tminus1 = varepsilontilde_t;
    // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
    dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube; 
    dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_t_mat; 
    // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell and \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell
    Rtilde_tminus1 = Rtilde_t;
    dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_t_cube;
    dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_t_cube;
  }
  arma::vec dmathcalLtilde_dvartheta = sum_dltilde_dvartheta / (n * 1.0); // \partial\widetilde\mathcal{L}_n(\bm\theta) / \partial\bm\vartheta
  return dmathcalLtilde_dvartheta;
}

////////////////////////////////////////////////// QMLE \widehat\bm\theta and ASD //////////////////////////////////////////////////
////// \bm\theta = (\bm\delta', \bm\beta')'
//////           = (\underline\bm\omega'(m*1), \bm\kappa', \beta_1, \beta_2, \underline\bmr'(m(m-1)/2*1))'
//////           = (\underline\bm\omega'(m*1), \bm\lambda'(r*1), \bm\gamma'(s*1), \bm\varphi'(s*1), \bm{g}_0'(rm^2*1), \bm{g}_1'(sm^2*1), \bm{g}_2'(sm^2*1), \beta_1, \beta_2, \underline\bm{r}'(m(m-1)/2*1))'
////// \bm{g}_0 = (\bm{g}_{01}', ..., \bm{g}_{0r}')', \bm{g}_1 = (\bm{g}_{11}', ..., \bm{g}_{1s}')', \bm{g}_2 = (\bm{g}_{21}', ..., \bm{g}_{2s}')' with 
////// \bm{g}_{0k} = vec(G_{0k}),
////// \bm{g}_{1k} = vec(G_{1k}),
////// \bm{g}_{2k} = vec(G_{2k}).
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////// QMLE \widehat\bm\theta //////////////////////////////////////////////////

arma::vec fc_varthetaTOtheta(int m, int r, int s, arma::vec varthetahat) {
  // r>0 and s>0, transfer \bm\vartheta to \bm\theta
  int dimkappa_in_vartheta = r + 2*s + 2*r*m + 8*s*m;
  int dimbeta = 2 + m*(m-1)/2;
  // int dimdelta_in_vartheta = m + dimkappa_in_vartheta;
  // int dimvartheta = dimdelta_in_vartheta + dimbeta;
  arma::vec omegauline = varthetahat.head(m); // \underline{\bm\omega}
  arma::vec kappa_in_vartheta = varthetahat.subvec(m, m+dimkappa_in_vartheta-1); // \bm\kappa in \bm\vartheta
  arma::vec lambda_gamma_varphi = kappa_in_vartheta.subvec(0, r+2*s-1); // (\bm\lambda', \bm\gamma', \bm\varphi')'
  arma::vec g_0_in_vartheta = kappa_in_vartheta.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_in_vartheta_mat = fc_asmat(g_0_in_vartheta, 2*m, r); // vector \bm{g}_0 = (\bm{g}_{01}', ..., \bm{g}_{0r}')' and matrix (\bm{g}_{01}, ..., \bm{g}_{0r}) in \bm\vartheta
  arma::vec g_1_in_vartheta = kappa_in_vartheta.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_in_vartheta_mat = fc_asmat(g_1_in_vartheta, 4*m, s); // vector \bm{g}_1 = (\bm{g}_{11}', ..., \bm{g}_{1s}')' and matrix (\bm{g}_{11}, ..., \bm{g}_{1s}) in \bm\vartheta
  arma::vec g_2_in_vartheta = kappa_in_vartheta.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_in_vartheta_mat = fc_asmat(g_2_in_vartheta, 4*m, s); // vector \bm{g}_2 = (\bm{g}_{21}', ..., \bm{g}_{2s}')' and matrix (\bm{g}_{21}, ..., \bm{g}_{2s}) in \bm\vartheta
  arma::vec beta_vec = varthetahat.tail(dimbeta); // \bm\beta
  int dimkappa_in_theta = r + 2*s + r*m*m + 2*s*m*m;
  int dimdelta_in_theta = m + dimkappa_in_theta;
  int dimtheta = dimdelta_in_theta + dimbeta;
  arma::mat g_0_in_theta_mat(m*m, r); g_0_in_theta_mat.fill(0.0); // matrix (\bm{g}_{01}, ..., \bm{g}_{0r}) in \bm\theta
  arma::mat g_1_in_theta_mat(m*m, s); g_1_in_theta_mat.fill(0.0); // matrix (\bm{g}_{11}, ..., \bm{g}_{1s}) in \bm\theta
  arma::mat g_2_in_theta_mat(m*m, s); g_2_in_theta_mat.fill(0.0); // matrix (\bm{g}_{21}, ..., \bm{g}_{2s}) in \bm\theta
  for(int k = 0; k < r; k++) {
    arma::vec g_0k_in_vartheta = g_0_in_vartheta_mat.col(k); // \bm{g}_{0k} = (\bm{g}_{0k1}', \bm{g}_{0k2}')' in \bm\vartheta
    arma::mat g_0k1_in_vartheta(m, 1); g_0k1_in_vartheta.col(0) = g_0k_in_vartheta.head(m); // \bm{g}_{0k1} in \bm\vartheta
    arma::mat g_0k2_in_vartheta(m, 1); g_0k2_in_vartheta.col(0) = g_0k_in_vartheta.tail(m); // \bm{g}_{0k2} in \bm\vartheta
    arma::mat G_0k = g_0k1_in_vartheta * g_0k2_in_vartheta.t(); // G_{0,k}
    arma::vec g_0k_in_theta = fc_asvec(G_0k); // \bm{g}_{0k} = vec(G_{0k}) in \bm\theta
    g_0_in_theta_mat.col(k) = g_0k_in_theta;
  }
  for(int k = 0; k < s; k++) {
    arma::vec g_1k_in_vartheta = g_1_in_vartheta_mat.col(k); // \bm{g}_{1k} = (\bm{g}_{1k1}', \bm{g}_{1k2}', \bm{g}_{1k3}', \bm{g}_{1k4}')' in \bm\vartheta
    arma::vec g_2k_in_vartheta = g_2_in_vartheta_mat.col(k); // \bm{g}_{2k} = (\bm{g}_{2k1}', \bm{g}_{2k2}', \bm{g}_{2k3}', \bm{g}_{2k4}')' in \bm\vartheta
    arma::mat g_1k1_in_vartheta(m, 1); g_1k1_in_vartheta.col(0) = g_1k_in_vartheta.head(m); // \bm{g}_{1k1} in \bm\vartheta
    arma::mat g_1k2_in_vartheta(m, 1); g_1k2_in_vartheta.col(0) = g_1k_in_vartheta.subvec(m, 2*m-1); // \bm{g}_{1k2} in \bm\vartheta
    arma::mat g_1k3_in_vartheta(m, 1); g_1k3_in_vartheta.col(0) = g_1k_in_vartheta.subvec(2*m, 3*m-1); // \bm{g}_{1k3} in \bm\vartheta
    arma::mat g_1k4_in_vartheta(m, 1); g_1k4_in_vartheta.col(0) = g_1k_in_vartheta.tail(m); // \bm{g}_{1k4} in \bm\vartheta
    arma::mat g_2k1_in_vartheta(m, 1); g_2k1_in_vartheta.col(0) = g_2k_in_vartheta.head(m); // \bm{g}_{2k1} in \bm\vartheta
    arma::mat g_2k2_in_vartheta(m, 1); g_2k2_in_vartheta.col(0) = g_2k_in_vartheta.subvec(m, 2*m-1); // \bm{g}_{2k2} in \bm\vartheta
    arma::mat g_2k3_in_vartheta(m, 1); g_2k3_in_vartheta.col(0) = g_2k_in_vartheta.subvec(2*m, 3*m-1); // \bm{g}_{2k3} in \bm\vartheta
    arma::mat g_2k4_in_vartheta(m, 1); g_2k4_in_vartheta.col(0) = g_2k_in_vartheta.tail(m); // \bm{g}_{2k4} in \bm\vartheta
    arma::mat G_1k = g_1k1_in_vartheta * g_1k2_in_vartheta.t() + g_1k3_in_vartheta * g_1k4_in_vartheta.t(); // G_{1,k}
    arma::mat G_2k = g_2k1_in_vartheta * g_2k2_in_vartheta.t() + g_2k3_in_vartheta * g_2k4_in_vartheta.t(); // G_{2,k}
    arma::vec g_1k_in_theta = fc_asvec(G_1k); // \bm{g}_{1k} = vec(G_{1k}) in \bm\theta
    arma::vec g_2k_in_theta = fc_asvec(G_2k); // \bm{g}_{2k} = vec(G_{2k}) in \bm\theta
    g_1_in_theta_mat.col(k) = g_1k_in_theta;
    g_2_in_theta_mat.col(k) = g_2k_in_theta;
  }
  arma::vec g_0_in_theta = fc_asvec(g_0_in_theta_mat); // (\bm{g}_{01}', ..., \bm{g}_{0r}')'
  arma::vec g_1_in_theta = fc_asvec(g_1_in_theta_mat); // (\bm{g}_{11}', ..., \bm{g}_{1s}')'
  arma::vec g_2_in_theta = fc_asvec(g_2_in_theta_mat); // (\bm{g}_{21}', ..., \bm{g}_{2s}')'
  arma::vec thetahat(dimtheta);
  thetahat.head(m) = omegauline;
  thetahat.subvec(m, m+r+2*s-1) = lambda_gamma_varphi;
  thetahat.subvec(m+r+2*s, m+r+2*s+r*m*m-1) = g_0_in_theta;
  thetahat.subvec(m+r+2*s+r*m*m, m+r+2*s+r*m*m+s*m*m-1) = g_1_in_theta;
  thetahat.subvec(m+r+2*s+r*m*m+s*m*m, m+r+2*s+r*m*m+2*s*m*m-1) = g_2_in_theta;
  thetahat.tail(dimbeta) = beta_vec;
  return thetahat;
}

arma::vec fc_varthetaTOtheta_requal0(int m, int r, int s, arma::vec varthetahat) {
  // r=0 and s>0, transfer \bm\vartheta to \bm\theta
  int dimkappa_in_vartheta = r + 2*s + 2*r*m + 8*s*m;
  int dimbeta = 2 + m*(m-1)/2;
  // int dimdelta_in_vartheta = m + dimkappa_in_vartheta;
  // int dimvartheta = dimdelta_in_vartheta + dimbeta;
  arma::vec omegauline = varthetahat.head(m); // \underline{\bm\omega}
  arma::vec kappa_in_vartheta = varthetahat.subvec(m, m+dimkappa_in_vartheta-1); // \bm\kappa in \bm\vartheta
  arma::vec lambda_gamma_varphi = kappa_in_vartheta.subvec(0, r+2*s-1); // (\bm\lambda', \bm\gamma', \bm\varphi')'
  // arma::vec g_0_in_vartheta = kappa_in_vartheta.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_in_vartheta_mat = fc_asmat(g_0_in_vartheta, 2*m, r); // vector \bm{g}_0 = (\bm{g}_{01}', ..., \bm{g}_{0r}')' and matrix (\bm{g}_{01}, ..., \bm{g}_{0r}) in \bm\vartheta
  arma::vec g_1_in_vartheta = kappa_in_vartheta.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_in_vartheta_mat = fc_asmat(g_1_in_vartheta, 4*m, s); // vector \bm{g}_1 = (\bm{g}_{11}', ..., \bm{g}_{1s}')' and matrix (\bm{g}_{11}, ..., \bm{g}_{1s}) in \bm\vartheta
  arma::vec g_2_in_vartheta = kappa_in_vartheta.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_in_vartheta_mat = fc_asmat(g_2_in_vartheta, 4*m, s); // vector \bm{g}_2 = (\bm{g}_{21}', ..., \bm{g}_{2s}')' and matrix (\bm{g}_{21}, ..., \bm{g}_{2s}) in \bm\vartheta
  arma::vec beta_vec = varthetahat.tail(dimbeta); // \bm\beta
  int dimkappa_in_theta = r + 2*s + r*m*m + 2*s*m*m;
  int dimdelta_in_theta = m + dimkappa_in_theta;
  int dimtheta = dimdelta_in_theta + dimbeta;
  // arma::mat g_0_in_theta_mat(m*m, r); g_0_in_theta_mat.fill(0.0); // matrix (\bm{g}_{01}, ..., \bm{g}_{0r}) in \bm\theta
  arma::mat g_1_in_theta_mat(m*m, s); g_1_in_theta_mat.fill(0.0); // matrix (\bm{g}_{11}, ..., \bm{g}_{1s}) in \bm\theta
  arma::mat g_2_in_theta_mat(m*m, s); g_2_in_theta_mat.fill(0.0); // matrix (\bm{g}_{21}, ..., \bm{g}_{2s}) in \bm\theta
  // for(int k = 0; k < r; k++) {
  //   arma::vec g_0k_in_vartheta = g_0_in_vartheta_mat.col(k); // \bm{g}_{0k} = (\bm{g}_{0k1}', \bm{g}_{0k2}')' in \bm\vartheta
  //   arma::mat g_0k1_in_vartheta(m, 1); g_0k1_in_vartheta.col(0) = g_0k_in_vartheta.head(m); // \bm{g}_{0k1} in \bm\vartheta
  //   arma::mat g_0k2_in_vartheta(m, 1); g_0k2_in_vartheta.col(0) = g_0k_in_vartheta.tail(m); // \bm{g}_{0k2} in \bm\vartheta
  //   arma::mat G_0k = g_0k1_in_vartheta * g_0k2_in_vartheta.t(); // G_{0,k}
  //   arma::vec g_0k_in_theta = fc_asvec(G_0k); // \bm{g}_{0k} = vec(G_{0k}) in \bm\theta
  //   g_0_in_theta_mat.col(k) = g_0k_in_theta;
  // }
  for(int k = 0; k < s; k++) {
    arma::vec g_1k_in_vartheta = g_1_in_vartheta_mat.col(k); // \bm{g}_{1k} = (\bm{g}_{1k1}', \bm{g}_{1k2}', \bm{g}_{1k3}', \bm{g}_{1k4}')' in \bm\vartheta
    arma::vec g_2k_in_vartheta = g_2_in_vartheta_mat.col(k); // \bm{g}_{2k} = (\bm{g}_{2k1}', \bm{g}_{2k2}', \bm{g}_{2k3}', \bm{g}_{2k4}')' in \bm\vartheta
    arma::mat g_1k1_in_vartheta(m, 1); g_1k1_in_vartheta.col(0) = g_1k_in_vartheta.head(m); // \bm{g}_{1k1} in \bm\vartheta
    arma::mat g_1k2_in_vartheta(m, 1); g_1k2_in_vartheta.col(0) = g_1k_in_vartheta.subvec(m, 2*m-1); // \bm{g}_{1k2} in \bm\vartheta
    arma::mat g_1k3_in_vartheta(m, 1); g_1k3_in_vartheta.col(0) = g_1k_in_vartheta.subvec(2*m, 3*m-1); // \bm{g}_{1k3} in \bm\vartheta
    arma::mat g_1k4_in_vartheta(m, 1); g_1k4_in_vartheta.col(0) = g_1k_in_vartheta.tail(m); // \bm{g}_{1k4} in \bm\vartheta
    arma::mat g_2k1_in_vartheta(m, 1); g_2k1_in_vartheta.col(0) = g_2k_in_vartheta.head(m); // \bm{g}_{2k1} in \bm\vartheta
    arma::mat g_2k2_in_vartheta(m, 1); g_2k2_in_vartheta.col(0) = g_2k_in_vartheta.subvec(m, 2*m-1); // \bm{g}_{2k2} in \bm\vartheta
    arma::mat g_2k3_in_vartheta(m, 1); g_2k3_in_vartheta.col(0) = g_2k_in_vartheta.subvec(2*m, 3*m-1); // \bm{g}_{2k3} in \bm\vartheta
    arma::mat g_2k4_in_vartheta(m, 1); g_2k4_in_vartheta.col(0) = g_2k_in_vartheta.tail(m); // \bm{g}_{2k4} in \bm\vartheta
    arma::mat G_1k = g_1k1_in_vartheta * g_1k2_in_vartheta.t() + g_1k3_in_vartheta * g_1k4_in_vartheta.t(); // G_{1,k}
    arma::mat G_2k = g_2k1_in_vartheta * g_2k2_in_vartheta.t() + g_2k3_in_vartheta * g_2k4_in_vartheta.t(); // G_{2,k}
    arma::vec g_1k_in_theta = fc_asvec(G_1k); // \bm{g}_{1k} = vec(G_{1k}) in \bm\theta
    arma::vec g_2k_in_theta = fc_asvec(G_2k); // \bm{g}_{2k} = vec(G_{2k}) in \bm\theta
    g_1_in_theta_mat.col(k) = g_1k_in_theta;
    g_2_in_theta_mat.col(k) = g_2k_in_theta;
  }
  // arma::vec g_0_in_theta = fc_asvec(g_0_in_theta_mat); // (\bm{g}_{01}', ..., \bm{g}_{0r}')'
  arma::vec g_1_in_theta = fc_asvec(g_1_in_theta_mat); // (\bm{g}_{11}', ..., \bm{g}_{1s}')'
  arma::vec g_2_in_theta = fc_asvec(g_2_in_theta_mat); // (\bm{g}_{21}', ..., \bm{g}_{2s}')'
  arma::vec thetahat(dimtheta);
  thetahat.head(m) = omegauline;
  thetahat.subvec(m, m+r+2*s-1) = lambda_gamma_varphi;
  // thetahat.subvec(m+r+2*s, m+r+2*s+r*m*m-1) = g_0_in_theta;
  thetahat.subvec(m+r+2*s+r*m*m, m+r+2*s+r*m*m+s*m*m-1) = g_1_in_theta;
  thetahat.subvec(m+r+2*s+r*m*m+s*m*m, m+r+2*s+r*m*m+2*s*m*m-1) = g_2_in_theta;
  thetahat.tail(dimbeta) = beta_vec;
  return thetahat;
}

arma::vec fc_varthetaTOtheta_sequal0(int m, int r, int s, arma::vec varthetahat) {
  // r>0 and s=0, transfer \bm\vartheta to \bm\theta
  int dimkappa_in_vartheta = r + 2*s + 2*r*m + 8*s*m;
  int dimbeta = 2 + m*(m-1)/2;
  // int dimdelta_in_vartheta = m + dimkappa_in_vartheta;
  // int dimvartheta = dimdelta_in_vartheta + dimbeta;
  arma::vec omegauline = varthetahat.head(m); // \underline{\bm\omega}
  arma::vec kappa_in_vartheta = varthetahat.subvec(m, m+dimkappa_in_vartheta-1); // \bm\kappa in \bm\vartheta
  arma::vec lambda_gamma_varphi = kappa_in_vartheta.subvec(0, r+2*s-1); // (\bm\lambda', \bm\gamma', \bm\varphi')'
  arma::vec g_0_in_vartheta = kappa_in_vartheta.subvec(r+2*s, r+2*s+2*r*m-1); arma::mat g_0_in_vartheta_mat = fc_asmat(g_0_in_vartheta, 2*m, r); // vector \bm{g}_0 = (\bm{g}_{01}', ..., \bm{g}_{0r}')' and matrix (\bm{g}_{01}, ..., \bm{g}_{0r}) in \bm\vartheta
  // arma::vec g_1_in_vartheta = kappa_in_vartheta.subvec(r+2*s+2*r*m, r+2*s+2*r*m+4*s*m-1); arma::mat g_1_in_vartheta_mat = fc_asmat(g_1_in_vartheta, 4*m, s); // vector \bm{g}_1 = (\bm{g}_{11}', ..., \bm{g}_{1s}')' and matrix (\bm{g}_{11}, ..., \bm{g}_{1s}) in \bm\vartheta
  // arma::vec g_2_in_vartheta = kappa_in_vartheta.subvec(r+2*s+2*r*m+4*s*m, r+2*s+2*r*m+8*s*m-1); arma::mat g_2_in_vartheta_mat = fc_asmat(g_2_in_vartheta, 4*m, s); // vector \bm{g}_2 = (\bm{g}_{21}', ..., \bm{g}_{2s}')' and matrix (\bm{g}_{21}, ..., \bm{g}_{2s}) in \bm\vartheta
  arma::vec beta_vec = varthetahat.tail(dimbeta); // \bm\beta
  int dimkappa_in_theta = r + 2*s + r*m*m + 2*s*m*m;
  int dimdelta_in_theta = m + dimkappa_in_theta;
  int dimtheta = dimdelta_in_theta + dimbeta;
  arma::mat g_0_in_theta_mat(m*m, r); g_0_in_theta_mat.fill(0.0); // matrix (\bm{g}_{01}, ..., \bm{g}_{0r}) in \bm\theta
  // arma::mat g_1_in_theta_mat(m*m, s); g_1_in_theta_mat.fill(0.0); // matrix (\bm{g}_{11}, ..., \bm{g}_{1s}) in \bm\theta
  // arma::mat g_2_in_theta_mat(m*m, s); g_2_in_theta_mat.fill(0.0); // matrix (\bm{g}_{21}, ..., \bm{g}_{2s}) in \bm\theta
  for(int k = 0; k < r; k++) {
    arma::vec g_0k_in_vartheta = g_0_in_vartheta_mat.col(k); // \bm{g}_{0k} = (\bm{g}_{0k1}', \bm{g}_{0k2}')' in \bm\vartheta
    arma::mat g_0k1_in_vartheta(m, 1); g_0k1_in_vartheta.col(0) = g_0k_in_vartheta.head(m); // \bm{g}_{0k1} in \bm\vartheta
    arma::mat g_0k2_in_vartheta(m, 1); g_0k2_in_vartheta.col(0) = g_0k_in_vartheta.tail(m); // \bm{g}_{0k2} in \bm\vartheta
    arma::mat G_0k = g_0k1_in_vartheta * g_0k2_in_vartheta.t(); // G_{0,k}
    arma::vec g_0k_in_theta = fc_asvec(G_0k); // \bm{g}_{0k} = vec(G_{0k}) in \bm\theta
    g_0_in_theta_mat.col(k) = g_0k_in_theta;
  }
  // for(int k = 0; k < s; k++) {
  //   arma::vec g_1k_in_vartheta = g_1_in_vartheta_mat.col(k); // \bm{g}_{1k} = (\bm{g}_{1k1}', \bm{g}_{1k2}', \bm{g}_{1k3}', \bm{g}_{1k4}')' in \bm\vartheta
  //   arma::vec g_2k_in_vartheta = g_2_in_vartheta_mat.col(k); // \bm{g}_{2k} = (\bm{g}_{2k1}', \bm{g}_{2k2}', \bm{g}_{2k3}', \bm{g}_{2k4}')' in \bm\vartheta
  //   arma::mat g_1k1_in_vartheta(m, 1); g_1k1_in_vartheta.col(0) = g_1k_in_vartheta.head(m); // \bm{g}_{1k1} in \bm\vartheta
  //   arma::mat g_1k2_in_vartheta(m, 1); g_1k2_in_vartheta.col(0) = g_1k_in_vartheta.subvec(m, 2*m-1); // \bm{g}_{1k2} in \bm\vartheta
  //   arma::mat g_1k3_in_vartheta(m, 1); g_1k3_in_vartheta.col(0) = g_1k_in_vartheta.subvec(2*m, 3*m-1); // \bm{g}_{1k3} in \bm\vartheta
  //   arma::mat g_1k4_in_vartheta(m, 1); g_1k4_in_vartheta.col(0) = g_1k_in_vartheta.tail(m); // \bm{g}_{1k4} in \bm\vartheta
  //   arma::mat g_2k1_in_vartheta(m, 1); g_2k1_in_vartheta.col(0) = g_2k_in_vartheta.head(m); // \bm{g}_{2k1} in \bm\vartheta
  //   arma::mat g_2k2_in_vartheta(m, 1); g_2k2_in_vartheta.col(0) = g_2k_in_vartheta.subvec(m, 2*m-1); // \bm{g}_{2k2} in \bm\vartheta
  //   arma::mat g_2k3_in_vartheta(m, 1); g_2k3_in_vartheta.col(0) = g_2k_in_vartheta.subvec(2*m, 3*m-1); // \bm{g}_{2k3} in \bm\vartheta
  //   arma::mat g_2k4_in_vartheta(m, 1); g_2k4_in_vartheta.col(0) = g_2k_in_vartheta.tail(m); // \bm{g}_{2k4} in \bm\vartheta
  //   arma::mat G_1k = g_1k1_in_vartheta * g_1k2_in_vartheta.t() + g_1k3_in_vartheta * g_1k4_in_vartheta.t(); // G_{1,k}
  //   arma::mat G_2k = g_2k1_in_vartheta * g_2k2_in_vartheta.t() + g_2k3_in_vartheta * g_2k4_in_vartheta.t(); // G_{2,k}
  //   arma::vec g_1k_in_theta = fc_asvec(G_1k); // \bm{g}_{1k} = vec(G_{1k}) in \bm\theta
  //   arma::vec g_2k_in_theta = fc_asvec(G_2k); // \bm{g}_{2k} = vec(G_{2k}) in \bm\theta
  //   g_1_in_theta_mat.col(k) = g_1k_in_theta;
  //   g_2_in_theta_mat.col(k) = g_2k_in_theta;
  // }
  arma::vec g_0_in_theta = fc_asvec(g_0_in_theta_mat); // (\bm{g}_{01}', ..., \bm{g}_{0r}')'
  // arma::vec g_1_in_theta = fc_asvec(g_1_in_theta_mat); // (\bm{g}_{11}', ..., \bm{g}_{1s}')'
  // arma::vec g_2_in_theta = fc_asvec(g_2_in_theta_mat); // (\bm{g}_{21}', ..., \bm{g}_{2s}')'
  arma::vec thetahat(dimtheta);
  thetahat.head(m) = omegauline;
  thetahat.subvec(m, m+r+2*s-1) = lambda_gamma_varphi;
  thetahat.subvec(m+r+2*s, m+r+2*s+r*m*m-1) = g_0_in_theta;
  // thetahat.subvec(m+r+2*s+r*m*m, m+r+2*s+r*m*m+s*m*m-1) = g_1_in_theta;
  // thetahat.subvec(m+r+2*s+r*m*m+s*m*m, m+r+2*s+r*m*m+2*s*m*m-1) = g_2_in_theta;
  thetahat.tail(dimbeta) = beta_vec;
  return thetahat;
}

// [[Rcpp::export]]
arma::vec fc_varthetaTOtheta_general(int m, int r, int s, arma::vec varthetahat) {
  // transfer \bm\vartheta to \bm\theta
  int dimtheta = varthetahat.n_elem;
  arma::vec thetahat(dimtheta);
  if(s == 0) {
    thetahat = fc_varthetaTOtheta_sequal0(m, r, s, varthetahat);
  } else if(r == 0) {
    thetahat = fc_varthetaTOtheta_requal0(m, r, s, varthetahat);
  } else {
    thetahat = fc_varthetaTOtheta(m, r, s, varthetahat);
  }
  return thetahat;
}

////////////////////////////////////////////////// ASD //////////////////////////////////////////////////

// arma::mat fc_general_Phi_i(int i, int m, int r, int s, arma::vec kappa_vec) {
//   // r>0 and s>0, \Phi_i(\bm\kappa) = \sum_{k=1}^r \lambda_k^{i-1} G_{0,k} + \sum_{k=1}^s \gamma_k^{i-1} [\cos((i-1)\varphi_k) G_{1,k} + \sin((i-1)\varphi_k) G_{2,k}]
//   arma::vec lambda_vec = kappa_vec.head(r); 
//   arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
//   arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
//   arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
//   arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
//   arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
//   arma::mat sum_lambdaG0(m, m); sum_lambdaG0.fill(0.0);
//   arma::mat sum_gammaG12(m, m); sum_gammaG12.fill(0.0);
//   for(int k = 0; k < r; k++) {
//     double lambda_k = lambda_vec(k);
//     arma::vec g_0k = g_0_mat.col(k);
//     arma::mat G_0k = fc_asmat(g_0k, m, m);
//     sum_lambdaG0 = sum_lambdaG0 + pow(lambda_k, i-1) * G_0k;
//   }
//   for(int k = 0; k < s; k++) {
//     double gamma_k = gamma_vec(k);
//     double varphi_k = varphi_vec(k);
//     arma::vec g_1k = g_1_mat.col(k);
//     arma::vec g_2k = g_2_mat.col(k);
//     arma::mat G_1k = fc_asmat(g_1k, m, m);
//     arma::mat G_2k = fc_asmat(g_2k, m, m);
//     sum_gammaG12 = sum_gammaG12 + pow(gamma_k, i-1) * (cos((i-1) * varphi_k) * G_1k + sin((i-1) * varphi_k) * G_2k);
//   }
//   arma::mat Phi_i = sum_lambdaG0 + sum_gammaG12; 
//   return Phi_i;
// }
// 
// arma::mat fc_general_Phi_i_requal0(int i, int m, int r, int s, arma::vec kappa_vec) {
//   // r=0 and s>0, \Phi_i(\bm\kappa) = \sum_{k=1}^s \gamma_k^{i-1} [\cos((i-1)\varphi_k) G_{1,k} + \sin((i-1)\varphi_k) G_{2,k}]
//   // arma::vec lambda_vec = kappa_vec.head(r); 
//   arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
//   arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
//   // arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
//   arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
//   arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
//   arma::mat sum_lambdaG0(m, m); sum_lambdaG0.fill(0.0);
//   arma::mat sum_gammaG12(m, m); sum_gammaG12.fill(0.0);
//   // for(int k = 0; k < r; k++) {
//   //   double lambda_k = lambda_vec(k);
//   //   arma::vec g_0k = g_0_mat.col(k);
//   //   arma::mat G_0k = fc_asmat(g_0k, m, m);
//   //   sum_lambdaG0 = sum_lambdaG0 + pow(lambda_k, i-1) * G_0k;
//   // }
//   for(int k = 0; k < s; k++) {
//     double gamma_k = gamma_vec(k);
//     double varphi_k = varphi_vec(k);
//     arma::vec g_1k = g_1_mat.col(k);
//     arma::vec g_2k = g_2_mat.col(k);
//     arma::mat G_1k = fc_asmat(g_1k, m, m);
//     arma::mat G_2k = fc_asmat(g_2k, m, m);
//     sum_gammaG12 = sum_gammaG12 + pow(gamma_k, i-1) * (cos((i-1) * varphi_k) * G_1k + sin((i-1) * varphi_k) * G_2k);
//   }
//   arma::mat Phi_i = sum_lambdaG0 + sum_gammaG12; 
//   return Phi_i;
// }
// 
// arma::mat fc_general_Phi_i_sequal0(int i, int m, int r, int s, arma::vec kappa_vec) {
//   // r>0 and s=0, \Phi_i(\bm\kappa) = \sum_{k=1}^r \lambda_k^{i-1} G_{0,k}
//   arma::vec lambda_vec = kappa_vec.head(r); 
//   // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
//   // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
//   arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
//   // arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
//   // arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
//   arma::mat sum_lambdaG0(m, m); sum_lambdaG0.fill(0.0);
//   arma::mat sum_gammaG12(m, m); sum_gammaG12.fill(0.0);
//   for(int k = 0; k < r; k++) {
//     double lambda_k = lambda_vec(k);
//     arma::vec g_0k = g_0_mat.col(k);
//     arma::mat G_0k = fc_asmat(g_0k, m, m);
//     sum_lambdaG0 = sum_lambdaG0 + pow(lambda_k, i-1) * G_0k;
//   }
//   // for(int k = 0; k < s; k++) {
//   //   double gamma_k = gamma_vec(k);
//   //   double varphi_k = varphi_vec(k);
//   //   arma::vec g_1k = g_1_mat.col(k);
//   //   arma::vec g_2k = g_2_mat.col(k);
//   //   arma::mat G_1k = fc_asmat(g_1k, m, m);
//   //   arma::mat G_2k = fc_asmat(g_2k, m, m);
//   //   sum_gammaG12 = sum_gammaG12 + pow(gamma_k, i-1) * (cos((i-1) * varphi_k) * G_1k + sin((i-1) * varphi_k) * G_2k);
//   // }
//   arma::mat Phi_i = sum_lambdaG0 + sum_gammaG12; 
//   return Phi_i;
// }
// 
// arma::mat fc_general_Phi_i_general(int i, int m, int r, int s, arma::vec kappa_vec) {
//   // \Phi_i(\bm\kappa) = \sum_{k=1}^r \lambda_k^{i-1} G_{0,k} + \sum_{k=1}^s \gamma_k^{i-1} [\cos((i-1)\varphi_k) G_{1,k} + \sin((i-1)\varphi_k) G_{2,k}]
//   arma::mat Phi_i(m, m); Phi_i.fill(0.0);
//   if(s == 0) {
//     Phi_i = fc_general_Phi_i_sequal0(i, m, r, s, kappa_vec);
//   } else if(r == 0) {
//     Phi_i = fc_general_Phi_i_requal0(i, m, r, s, kappa_vec);
//   } else {
//     Phi_i = fc_general_Phi_i(i, m, r, s, kappa_vec);
//   }
//   return Phi_i;
// }

// arma::vec fc_general_lnhulinetilde_t(int t, int m, int r, int s, arma::vec omegauline, arma::vec kappa_vec, arma::mat y_m_n) {
//   // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = (\ln\widetilde{h}_{11,t}(\bm\delta), ..., \ln\widetilde{h}_{mm,t}(\bm\delta))'
//   arma::mat Philnyuline_mat(m, t-1); Philnyuline_mat.fill(0.0); // (\Phi_1(\bm\kappa) \ln\underline{\bm{y}}_{t-1}, ..., \Phi_{t-1}(\bm\kappa) \ln\underline{\bm{y}}_{1})
//   for(int i = 1; i < t; i++) {
//     arma::mat Phi_i = fc_general_Phi_i_general(i, m, r, s, kappa_vec); // \Phi_i(\bm\kappa)
//     arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
//     arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln y_{1,t-i}^2, ..., \ln y_{m,t-i}^2)'
//     Philnyuline_mat.col(i-1) = Phi_i * lnyuline_tminusi; // \Phi_i(\bm\kappa) \ln\underline{\bm{y}}_{t-i}
//   }
//   arma::vec lnhulinetilde_t = omegauline + sum(Philnyuline_mat, 1); // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = \underline{\bm{\omega}} + \sum_{i=1}^{t-1} \Phi_i(\bm\kappa) \ln\underline{\bm{y}}_{t-i}
//   return lnhulinetilde_t;
// }

arma::vec fc_general_lnhulinetilde_t(int t, int m, int r, int s, arma::vec omegauline, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_0_m_nminus1_r, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s) {
  // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = (\ln\widetilde{h}_{11,t}(\bm\delta), ..., \ln\widetilde{h}_{mm,t}(\bm\delta))'
  arma::vec lnhulinetilde_t = omegauline; // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = \underline{\bm{\omega}} + \sum_{i=1}^{t-1} \Phi_i(\bm\kappa) \ln\underline{\bm{y}}_{t-i}
    // arma::vec lambda_vec = kappa_vec.head(r);
    // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
    // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
    arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
    arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
    arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
  arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  for(int k = 0; k < r; k++) {
    // double lambda_k = lambda_vec(k);
    arma::vec g_0k = g_0_mat.col(k);
    arma::mat G_0k = fc_asmat(g_0k, m, m);
    arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    lnhulinetilde_t = lnhulinetilde_t + G_0k * sum_lambdaklnyuline_0;
  }
  for(int k = 0; k < s; k++) {
    // double gamma_k = gamma_vec(k);
    // double varphi_k = varphi_vec(k);
    arma::vec g_1k = g_1_mat.col(k);
    arma::vec g_2k = g_2_mat.col(k);
    arma::mat G_1k = fc_asmat(g_1k, m, m);
    arma::mat G_2k = fc_asmat(g_2k, m, m);
    arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    lnhulinetilde_t = lnhulinetilde_t + G_1k * sum_gammaklnyuline_01 + G_2k * sum_gammaklnyuline_02;
  }
  return lnhulinetilde_t;
}

arma::vec fc_general_lnhulinetilde_t_requal0(int t, int m, int r, int s, arma::vec omegauline, arma::vec kappa_vec, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s) {
  // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = (\ln\widetilde{h}_{11,t}(\bm\delta), ..., \ln\widetilde{h}_{mm,t}(\bm\delta))'
  arma::vec lnhulinetilde_t = omegauline; // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = \underline{\bm{\omega}} + \sum_{i=1}^{t-1} \Phi_i(\bm\kappa) \ln\underline{\bm{y}}_{t-i}
  // arma::vec lambda_vec = kappa_vec.head(r);
  // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  // arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
  arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
  arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
  // arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  // for(int k = 0; k < r; k++) {
  //   // double lambda_k = lambda_vec(k);
  //   arma::vec g_0k = g_0_mat.col(k);
  //   arma::mat G_0k = fc_asmat(g_0k, m, m);
  //   arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
  //   lnhulinetilde_t = lnhulinetilde_t + G_0k * sum_lambdaklnyuline_0;
  // }
  for(int k = 0; k < s; k++) {
    // double gamma_k = gamma_vec(k);
    // double varphi_k = varphi_vec(k);
    arma::vec g_1k = g_1_mat.col(k);
    arma::vec g_2k = g_2_mat.col(k);
    arma::mat G_1k = fc_asmat(g_1k, m, m);
    arma::mat G_2k = fc_asmat(g_2k, m, m);
    arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    lnhulinetilde_t = lnhulinetilde_t + G_1k * sum_gammaklnyuline_01 + G_2k * sum_gammaklnyuline_02;
  }
  return lnhulinetilde_t;
}

arma::vec fc_general_lnhulinetilde_t_sequal0(int t, int m, int r, int s, arma::vec omegauline, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_0_m_nminus1_r) {
  // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = (\ln\widetilde{h}_{11,t}(\bm\delta), ..., \ln\widetilde{h}_{mm,t}(\bm\delta))'
  arma::vec lnhulinetilde_t = omegauline; // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = \underline{\bm{\omega}} + \sum_{i=1}^{t-1} \Phi_i(\bm\kappa) \ln\underline{\bm{y}}_{t-i}
  // arma::vec lambda_vec = kappa_vec.head(r);
  // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
  // arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
  // arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
  arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  // arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  for(int k = 0; k < r; k++) {
    // double lambda_k = lambda_vec(k);
    arma::vec g_0k = g_0_mat.col(k);
    arma::mat G_0k = fc_asmat(g_0k, m, m);
    arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    lnhulinetilde_t = lnhulinetilde_t + G_0k * sum_lambdaklnyuline_0;
  }
  // for(int k = 0; k < s; k++) {
  //   // double gamma_k = gamma_vec(k);
  //   // double varphi_k = varphi_vec(k);
  //   arma::vec g_1k = g_1_mat.col(k);
  //   arma::vec g_2k = g_2_mat.col(k);
  //   arma::mat G_1k = fc_asmat(g_1k, m, m);
  //   arma::mat G_2k = fc_asmat(g_2k, m, m);
  //   arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   lnhulinetilde_t = lnhulinetilde_t + G_1k * sum_gammaklnyuline_01 + G_2k * sum_gammaklnyuline_02;
  // }
  return lnhulinetilde_t;
}

arma::mat fc_general_lnhulinetilde_t_general(int t, int m, int r, int s, arma::vec omegauline, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_0_m_nminus1_r, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s) {
  // \ln\widetilde\underline{\bm{h}}_t(\bm\delta) = (\ln\widetilde{h}_{11,t}(\bm\delta), ..., \ln\widetilde{h}_{mm,t}(\bm\delta))'
  arma::vec lnhulinetilde_t(m);
  if(s == 0) {
    lnhulinetilde_t = fc_general_lnhulinetilde_t_sequal0(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r);
  } else if(r == 0) {
    lnhulinetilde_t = fc_general_lnhulinetilde_t_requal0(t, m, r, s, omegauline, kappa_vec, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s);
  } else {
    lnhulinetilde_t = fc_general_lnhulinetilde_t(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s);
  }
  return lnhulinetilde_t;
}

arma::mat fc_general_dlnhulinetilde_ddelta_1_mat(int m, int dimdelta) {
  // the \ell-th column is \partial\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\delta_\ell
  arma::mat dlnhulinetilde_ddelta_1_mat(m, dimdelta); dlnhulinetilde_ddelta_1_mat.fill(0.0);
  arma::mat I_m(m,m); I_m.eye(m,m);
  dlnhulinetilde_ddelta_1_mat.cols(0, m-1) = I_m; // \partial\ln\widetilde\underline{\bm{h}}_1(\bm\delta)/\partial\underline\omega_\ell = \bm{e}_\ell
  return dlnhulinetilde_ddelta_1_mat;
}

arma::mat fc_general_dlnhulinetilde_ddelta_t_mat(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_0_m_nminus1_r, arma::cube sum_lambdaklnyuline_1_m_nminus1_r, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s, arma::cube sum_gammaklnyuline_11_m_nminus1_s, arma::cube sum_gammaklnyuline_12_m_nminus1_s, arma::cube sum_gammaklnyuline_21_m_nminus1_s, arma::cube sum_gammaklnyuline_22_m_nminus1_s) {
  // r>0 and s>0, the \ell-th column is \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
  arma::vec lambda_vec = kappa_vec.head(r); 
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
  arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
  arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat zerosmat_m_1(m, 1); zerosmat_m_1.fill(0.0);
  int Dim1 = m+r+2*s; int Dim2 = Dim1 + r*m*m; int Dim3 = Dim2 + s*m*m;
  // arma::vec y_tminus1 = y_m_n.col(t-2); // \bm{y}_{t-1}
  // arma::vec lnyuline_tminus1 = log(pow(y_tminus1,2)); // \ln\underline{\bm{y}}_{t-1} = (\ln{y_{1,t-1}^2}, ..., \ln{y_{m,t-1}^2})'
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta
  arma::mat dlnhulinetilde_ddelta_t_mat(m, dimdelta); dlnhulinetilde_ddelta_t_mat.fill(0.0);
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\underline\omega
  dlnhulinetilde_ddelta_t_mat.cols(0, m-1) = I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta)/\partial\underline\omega_\ell = \bm{e}_\ell
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\lambda, \partial\bm{g}_{0}'
  arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  arma::mat sum_lambdaklnyuline_1_tminus1_m_r = sum_lambdaklnyuline_1_m_nminus1_r.col(t-2);
  for(int k = 0; k < r; k++) {
    // double lambda_k = lambda_vec(k);
    arma::vec g_0k = g_0_mat.col(k);
    arma::mat G_0k = fc_asmat(g_0k, m, m);
    // arma::vec sum_lambdaklnyuline_0 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_lambdaklnyuline_1 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    // for(int i = 2; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0 + pow(lambda_k, i-1) * lnyuline_tminusi;
    //   sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1 + (i-1) * pow(lambda_k, i-2) * lnyuline_tminusi;
    // }
    arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1_tminus1_m_r.col(k); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    dlnhulinetilde_ddelta_t_mat.col(m+k) = G_0k * sum_lambdaklnyuline_1; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k
    arma::mat sum_lambdaklnyuline_0_mat = fc_asmat(sum_lambdaklnyuline_0, m, 1);
    dlnhulinetilde_ddelta_t_mat.cols(Dim1 +k*m*m,   Dim1 +(k+1)*m*m-1)   = kron(sum_lambdaklnyuline_0_mat.t(), I_m); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{0,k}'
  }
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\gamma, \partial\bm\varphi, \partial\bm{g}_{1}', \partial\bm{g}_{2}'
  arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_11_tminus1_m_s = sum_gammaklnyuline_11_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_12_tminus1_m_s = sum_gammaklnyuline_12_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_21_tminus1_m_s = sum_gammaklnyuline_21_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_22_tminus1_m_s = sum_gammaklnyuline_22_m_nminus1_s.col(t-2);
  for(int k = 0; k < s; k++) {
    // double gamma_k = gamma_vec(k);
    // double varphi_k = varphi_vec(k);
    arma::vec g_1k = g_1_mat.col(k);
    arma::vec g_2k = g_2_mat.col(k);
    arma::mat G_1k = fc_asmat(g_1k, m, m);
    arma::mat G_2k = fc_asmat(g_2k, m, m);
    // arma::vec sum_gammaklnyuline_01 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_02 = zerosmat_m_1.col(0); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_11 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_12 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_21 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_22 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // for(int i = 2; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_gammaklnyuline_01 = sum_gammaklnyuline_01 + pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_02 = sum_gammaklnyuline_02 + pow(gamma_k, i-1) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_11 = sum_gammaklnyuline_11 + (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_12 = sum_gammaklnyuline_12 + (i-1) * pow(gamma_k, i-2) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_21 = sum_gammaklnyuline_21 + (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_22 = sum_gammaklnyuline_22 + (i-1) * pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    // }
    arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_11 = sum_gammaklnyuline_11_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_12 = sum_gammaklnyuline_12_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_21 = sum_gammaklnyuline_21_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_22 = sum_gammaklnyuline_22_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    dlnhulinetilde_ddelta_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_11 + G_2k * sum_gammaklnyuline_12; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k
    dlnhulinetilde_ddelta_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_21 + G_2k * sum_gammaklnyuline_22; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k
    arma::mat sum_gammaklnyuline_01_mat = fc_asmat(sum_gammaklnyuline_01, m, 1);
    arma::mat sum_gammaklnyuline_02_mat = fc_asmat(sum_gammaklnyuline_02, m, 1);
    dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*m*m,     Dim2 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_01_mat.t(), I_m); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*m*m,     Dim3 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_02_mat.t(), I_m); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k}'
  }
  return dlnhulinetilde_ddelta_t_mat;
}

arma::mat fc_general_dlnhulinetilde_ddelta_t_mat_requal0(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s, arma::cube sum_gammaklnyuline_11_m_nminus1_s, arma::cube sum_gammaklnyuline_12_m_nminus1_s, arma::cube sum_gammaklnyuline_21_m_nminus1_s, arma::cube sum_gammaklnyuline_22_m_nminus1_s) {
  // r=0 and s>0, the \ell-th column is \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
  // arma::vec lambda_vec = kappa_vec.head(r); 
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  // arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
  arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
  arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat zerosmat_m_1(m, 1); zerosmat_m_1.fill(0.0);
  int Dim1 = m+r+2*s; int Dim2 = Dim1 + r*m*m; int Dim3 = Dim2 + s*m*m;
  // arma::vec y_tminus1 = y_m_n.col(t-2); // \bm{y}_{t-1}
  // arma::vec lnyuline_tminus1 = log(pow(y_tminus1,2)); // \ln\underline{\bm{y}}_{t-1} = (\ln{y_{1,t-1}^2}, ..., \ln{y_{m,t-1}^2})'
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta
  arma::mat dlnhulinetilde_ddelta_t_mat(m, dimdelta); dlnhulinetilde_ddelta_t_mat.fill(0.0);
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\underline\omega
  dlnhulinetilde_ddelta_t_mat.cols(0, m-1) = I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta)/\partial\underline\omega_\ell = \bm{e}_\ell
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\lambda, \partial\bm{g}_{0}'
  // arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  // arma::mat sum_lambdaklnyuline_1_tminus1_m_r = sum_lambdaklnyuline_1_m_nminus1_r.col(t-2);
  // for(int k = 0; k < r; k++) {
  //   // double lambda_k = lambda_vec(k);
  //   arma::vec g_0k = g_0_mat.col(k);
  //   arma::mat G_0k = fc_asmat(g_0k, m, m);
  //   // arma::vec sum_lambdaklnyuline_0 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_lambdaklnyuline_1 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
  //   // for(int i = 2; i < t; i++) {
  //   //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
  //   //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
  //   //   sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0 + pow(lambda_k, i-1) * lnyuline_tminusi;
  //   //   sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1 + (i-1) * pow(lambda_k, i-2) * lnyuline_tminusi;
  //   // }
  //   arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1_tminus1_m_r.col(k); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
  //   dlnhulinetilde_ddelta_t_mat.col(m+k) = G_0k * sum_lambdaklnyuline_1; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k
  //   arma::mat sum_lambdaklnyuline_0_mat = fc_asmat(sum_lambdaklnyuline_0, m, 1);
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim1 +k*m*m,   Dim1 +(k+1)*m*m-1)   = kron(sum_lambdaklnyuline_0_mat.t(), I_m); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{0,k}'
  // }
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\gamma, \partial\bm\varphi, \partial\bm{g}_{1}', \partial\bm{g}_{2}'
  arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_11_tminus1_m_s = sum_gammaklnyuline_11_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_12_tminus1_m_s = sum_gammaklnyuline_12_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_21_tminus1_m_s = sum_gammaklnyuline_21_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_22_tminus1_m_s = sum_gammaklnyuline_22_m_nminus1_s.col(t-2);
  for(int k = 0; k < s; k++) {
    // double gamma_k = gamma_vec(k);
    // double varphi_k = varphi_vec(k);
    arma::vec g_1k = g_1_mat.col(k);
    arma::vec g_2k = g_2_mat.col(k);
    arma::mat G_1k = fc_asmat(g_1k, m, m);
    arma::mat G_2k = fc_asmat(g_2k, m, m);
    // arma::vec sum_gammaklnyuline_01 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_02 = zerosmat_m_1.col(0); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_11 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_12 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_21 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_22 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // for(int i = 2; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_gammaklnyuline_01 = sum_gammaklnyuline_01 + pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_02 = sum_gammaklnyuline_02 + pow(gamma_k, i-1) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_11 = sum_gammaklnyuline_11 + (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_12 = sum_gammaklnyuline_12 + (i-1) * pow(gamma_k, i-2) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_21 = sum_gammaklnyuline_21 + (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_22 = sum_gammaklnyuline_22 + (i-1) * pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    // }
    arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_11 = sum_gammaklnyuline_11_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_12 = sum_gammaklnyuline_12_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_21 = sum_gammaklnyuline_21_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_22 = sum_gammaklnyuline_22_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    dlnhulinetilde_ddelta_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_11 + G_2k * sum_gammaklnyuline_12; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k
    dlnhulinetilde_ddelta_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_21 + G_2k * sum_gammaklnyuline_22; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k
    arma::mat sum_gammaklnyuline_01_mat = fc_asmat(sum_gammaklnyuline_01, m, 1);
    arma::mat sum_gammaklnyuline_02_mat = fc_asmat(sum_gammaklnyuline_02, m, 1);
    dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*m*m,     Dim2 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_01_mat.t(), I_m); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k}'
    dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*m*m,     Dim3 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_02_mat.t(), I_m); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k}'
  }
  return dlnhulinetilde_ddelta_t_mat;
}

arma::mat fc_general_dlnhulinetilde_ddelta_t_mat_sequal0(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_0_m_nminus1_r, arma::cube sum_lambdaklnyuline_1_m_nminus1_r) {
  // r>0 and s=0, the \ell-th column is \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
  arma::vec lambda_vec = kappa_vec.head(r); 
  // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
  // arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
  // arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat zerosmat_m_1(m, 1); zerosmat_m_1.fill(0.0);
  int Dim1 = m+r+2*s; // int Dim2 = Dim1 + r*m*m; int Dim3 = Dim2 + s*m*m;
  // arma::vec y_tminus1 = y_m_n.col(t-2); // \bm{y}_{t-1}
  // arma::vec lnyuline_tminus1 = log(pow(y_tminus1,2)); // \ln\underline{\bm{y}}_{t-1} = (\ln{y_{1,t-1}^2}, ..., \ln{y_{m,t-1}^2})'
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta
  arma::mat dlnhulinetilde_ddelta_t_mat(m, dimdelta); dlnhulinetilde_ddelta_t_mat.fill(0.0);
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\underline\omega
  dlnhulinetilde_ddelta_t_mat.cols(0, m-1) = I_m; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta)/\partial\underline\omega_\ell = \bm{e}_\ell
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\lambda, \partial\bm{g}_{0}'
  arma::mat sum_lambdaklnyuline_0_tminus1_m_r = sum_lambdaklnyuline_0_m_nminus1_r.col(t-2);
  arma::mat sum_lambdaklnyuline_1_tminus1_m_r = sum_lambdaklnyuline_1_m_nminus1_r.col(t-2);
  for(int k = 0; k < r; k++) {
    // double lambda_k = lambda_vec(k);
    arma::vec g_0k = g_0_mat.col(k);
    arma::mat G_0k = fc_asmat(g_0k, m, m);
    // arma::vec sum_lambdaklnyuline_0 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_lambdaklnyuline_1 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    // for(int i = 2; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0 + pow(lambda_k, i-1) * lnyuline_tminusi;
    //   sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1 + (i-1) * pow(lambda_k, i-2) * lnyuline_tminusi;
    // }
    arma::vec sum_lambdaklnyuline_0 = sum_lambdaklnyuline_0_tminus1_m_r.col(k); // \sum_{i=1}^{t-1} \lambda_k^{i-1} \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1_tminus1_m_r.col(k); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    dlnhulinetilde_ddelta_t_mat.col(m+k) = G_0k * sum_lambdaklnyuline_1; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k
    arma::mat sum_lambdaklnyuline_0_mat = fc_asmat(sum_lambdaklnyuline_0, m, 1);
    dlnhulinetilde_ddelta_t_mat.cols(Dim1 +k*m*m,   Dim1 +(k+1)*m*m-1)   = kron(sum_lambdaklnyuline_0_mat.t(), I_m); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{0,k}'
  }
  // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\gamma, \partial\bm\varphi, \partial\bm{g}_{1}', \partial\bm{g}_{2}'
  // arma::mat sum_gammaklnyuline_01_tminus1_m_s = sum_gammaklnyuline_01_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_02_tminus1_m_s = sum_gammaklnyuline_02_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_11_tminus1_m_s = sum_gammaklnyuline_11_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_12_tminus1_m_s = sum_gammaklnyuline_12_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_21_tminus1_m_s = sum_gammaklnyuline_21_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_22_tminus1_m_s = sum_gammaklnyuline_22_m_nminus1_s.col(t-2);
  // for(int k = 0; k < s; k++) {
  //   // double gamma_k = gamma_vec(k);
  //   // double varphi_k = varphi_vec(k);
  //   arma::vec g_1k = g_1_mat.col(k);
  //   arma::vec g_2k = g_2_mat.col(k);
  //   arma::mat G_1k = fc_asmat(g_1k, m, m);
  //   arma::mat G_2k = fc_asmat(g_2k, m, m);
  //   // arma::vec sum_gammaklnyuline_01 = lnyuline_tminus1; // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_02 = zerosmat_m_1.col(0); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_11 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_12 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_21 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_22 = zerosmat_m_1.col(0); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // for(int i = 2; i < t; i++) {
  //   //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
  //   //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
  //   //   sum_gammaklnyuline_01 = sum_gammaklnyuline_01 + pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_02 = sum_gammaklnyuline_02 + pow(gamma_k, i-1) * sin((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_11 = sum_gammaklnyuline_11 + (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_12 = sum_gammaklnyuline_12 + (i-1) * pow(gamma_k, i-2) * sin((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_21 = sum_gammaklnyuline_21 + (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_22 = sum_gammaklnyuline_22 + (i-1) * pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
  //   // }
  //   arma::vec sum_gammaklnyuline_01 = sum_gammaklnyuline_01_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_02 = sum_gammaklnyuline_02_tminus1_m_s.col(k); // \sum_{i=1}^{t-1} \gamma_k^{i-1} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_11 = sum_gammaklnyuline_11_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_12 = sum_gammaklnyuline_12_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_21 = sum_gammaklnyuline_21_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_22 = sum_gammaklnyuline_22_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   dlnhulinetilde_ddelta_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_11 + G_2k * sum_gammaklnyuline_12; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k
  //   dlnhulinetilde_ddelta_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_21 + G_2k * sum_gammaklnyuline_22; // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k
  //   arma::mat sum_gammaklnyuline_01_mat = fc_asmat(sum_gammaklnyuline_01, m, 1);
  //   arma::mat sum_gammaklnyuline_02_mat = fc_asmat(sum_gammaklnyuline_02, m, 1);
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim2 +k*m*m,     Dim2 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_01_mat.t(), I_m); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{1,k}'
  //   dlnhulinetilde_ddelta_t_mat.cols(Dim3 +k*m*m,     Dim3 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_02_mat.t(), I_m); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm{g}_{2,k}'
  // }
  return dlnhulinetilde_ddelta_t_mat;
}

arma::mat fc_general_dlnhulinetilde_ddelta_t_mat_general(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_0_m_nminus1_r, arma::cube sum_lambdaklnyuline_1_m_nminus1_r, arma::cube sum_gammaklnyuline_01_m_nminus1_s, arma::cube sum_gammaklnyuline_02_m_nminus1_s, arma::cube sum_gammaklnyuline_11_m_nminus1_s, arma::cube sum_gammaklnyuline_12_m_nminus1_s, arma::cube sum_gammaklnyuline_21_m_nminus1_s, arma::cube sum_gammaklnyuline_22_m_nminus1_s) {
  // the \ell-th column is \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
  arma::mat dlnhulinetilde_ddelta_t_mat(m, dimdelta); dlnhulinetilde_ddelta_t_mat.fill(0.0);
  if(s == 0) {
    dlnhulinetilde_ddelta_t_mat = fc_general_dlnhulinetilde_ddelta_t_mat_sequal0(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_lambdaklnyuline_1_m_nminus1_r);
  } else if(r == 0) {
    dlnhulinetilde_ddelta_t_mat = fc_general_dlnhulinetilde_ddelta_t_mat_requal0(t, m, r, s, dimdelta, kappa_vec, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s);
  } else {
    dlnhulinetilde_ddelta_t_mat = fc_general_dlnhulinetilde_ddelta_t_mat(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_lambdaklnyuline_1_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s);
  }
  return dlnhulinetilde_ddelta_t_mat;
}

arma::cube fc_dDtilde_ddelta_t_cube(int m, int dimdelta, arma::mat Dtilde_t, arma::mat dlnhulinetilde_ddelta_t_mat) {
  // the \ell-th slice is \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
  arma::cube dDtilde_ddelta_t_cube(m, m, dimdelta); dDtilde_ddelta_t_cube.fill(0.0);
  for(int l = 0; l < dimdelta; l++) {
    arma::vec dlnhulinetilde_ddelta_t_l = dlnhulinetilde_ddelta_t_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
    dDtilde_ddelta_t_cube.slice(l) = fc_dDtilde_ddelta_t_l(Dtilde_t, dlnhulinetilde_ddelta_t_l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
  }
  return dDtilde_ddelta_t_cube;
}

arma::cube fc_dPsitilde_ddelta_tminus1_cube(int m, int dimdelta, arma::mat varepsilontilde_tminusBkTOtminus1_mat, arma::mat Psitilde_tminus1, arma::cube dlnhulinetilde_ddelta_tminusBkTOtminus1_cube) {
  // the \ell-th slice is \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
  arma::cube dPsitilde_ddelta_tminus1_cube(m, m, dimdelta); dPsitilde_ddelta_tminus1_cube.fill(0.0);
  for(int l = 0; l < dimdelta; l++) {
    arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.col(l); // the k-th column is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-k)}(\bm\delta) / \partial\delta_\ell
    dPsitilde_ddelta_tminus1_cube.slice(l) = fc_dPsitilde_ddelta_tminus1_l(m, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
  }
  return dPsitilde_ddelta_tminus1_cube;
}

arma::cube fc2_dRtilde_ddelta_t_cube(int m, int dimdelta, double beta_1, double beta_2, arma::cube dPsitilde_ddelta_tminus1_cube, arma::cube dRtilde_ddelta_tminus1_cube) {
  // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
  arma::cube dRtilde_ddelta_t_cube(m, m, dimdelta); dRtilde_ddelta_t_cube.fill(0.0);
  for(int l = 0; l < dimdelta; l++) {
    arma::mat dPsitilde_ddelta_tminus1_l = dPsitilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
    arma::mat dRtilde_ddelta_tminus1_l = dRtilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
    dRtilde_ddelta_t_cube.slice(l) = beta_1 * dPsitilde_ddelta_tminus1_l + beta_2 * dRtilde_ddelta_tminus1_l;
  }
  return dRtilde_ddelta_t_cube;
}

arma::cube fc_dHtilde_ddelta_t_cube(int m, int dimdelta, arma::mat Dtilde_t, arma::mat Rtilde_t, arma::cube dDtilde_ddelta_t_cube, arma::cube dRtilde_ddelta_t_cube) {
  // the \ell-th slice is \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
  arma::cube dHtilde_ddelta_t_cube(m, m, dimdelta); dHtilde_ddelta_t_cube.fill(0.0);
  for(int l = 0; l < dimdelta; l++) {
    arma::mat dDtilde_ddelta_t_l = dDtilde_ddelta_t_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
    arma::mat dRtilde_ddelta_t_l = dRtilde_ddelta_t_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
    dHtilde_ddelta_t_cube.slice(l) = fc_dHtilde_ddelta_t_l(Dtilde_t, Rtilde_t, dDtilde_ddelta_t_l, dRtilde_ddelta_t_l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
  }
  return dHtilde_ddelta_t_cube;
}

arma::cube fc_dHtilde_dbeta_t_cube(int m, int dimbeta, arma::mat Dtilde_t, arma::cube dRtilde_dbeta_t_cube) {
  // the \ell-th slice is \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
  arma::cube dHtilde_dbeta_t_cube(m, m, dimbeta); dHtilde_dbeta_t_cube.fill(0.0);
  for(int l = 0; l < dimbeta; l++) {
    arma::mat dRtilde_dbeta_t_l = dRtilde_dbeta_t_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_\ell
    dHtilde_dbeta_t_cube.slice(l) = fc_dHtilde_dbeta_t_l(Dtilde_t, dRtilde_dbeta_t_l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
  }
  return dHtilde_dbeta_t_cube;
}

arma::vec fc2_dltilde_dtheta_t(int m, int dimdelta, int dimbeta, int dimtheta, arma::vec y_t, arma::mat Htilde_t, arma::cube dHtilde_ddelta_t_cube, arma::cube dHtilde_dbeta_t_cube) {
  // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
  arma::vec dltilde_dtheta_t(dimtheta); // \partial\widetilde{\ell}_t(\bm\theta) / \partial\theta_\ell = \tr[(I_m - \widetilde{H}_t^{-1}(\bm\theta) \bm{y}_t \bm{y}_t') \widetilde{H}_t^{-1}(\bm\theta) \partial\widetilde{H}_t(\bm\theta)/\partial\theta_\ell]
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat y_t_mat = fc_asmat(y_t, m, 1);
  arma::mat inverse_Htilde_t = Htilde_t.i();
  arma::mat left_mat = (I_m - inverse_Htilde_t * y_t * y_t_mat.t()) * inverse_Htilde_t;
  for(int l = 0; l < dimdelta; l++) {
    arma::mat dHtilde_ddelta_t_l = dHtilde_ddelta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
    dltilde_dtheta_t(l) = trace(left_mat * dHtilde_ddelta_t_l); // \partial\widetilde{\ell}_t(\bm\theta) / \partial\delta_\ell
  }
  for(int l = 0; l < dimbeta; l++) {
    arma::mat dHtilde_dbeta_t_l = dHtilde_dbeta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
    dltilde_dtheta_t(dimdelta+l) = trace(left_mat * dHtilde_dbeta_t_l); // \partial\widetilde{\ell}_t(\bm\theta) / \partial\beta_\ell
  }
  return dltilde_dtheta_t;
}

arma::cube fc_general_ddlnhulinetilde_ddeltaddeltaprime_1_cube(int m, int dimdelta) {
  // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\delta_k\partial\delta_\ell = \bm{0}_m
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_1_cube(m, dimdelta, dimdelta); ddlnhulinetilde_ddeltaddeltaprime_1_cube.fill(0.0);
  return ddlnhulinetilde_ddeltaddeltaprime_1_cube;
}

arma::cube fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_1_m_nminus1_r, arma::cube sum_lambdaklnyuline_2_m_nminus1_r, arma::cube sum_gammaklnyuline_11_m_nminus1_s, arma::cube sum_gammaklnyuline_12_m_nminus1_s, arma::cube sum_gammaklnyuline_21_m_nminus1_s, arma::cube sum_gammaklnyuline_22_m_nminus1_s, arma::cube sum_gammaklnyuline_31_m_nminus1_s, arma::cube sum_gammaklnyuline_32_m_nminus1_s, arma::cube sum_gammaklnyuline_41_m_nminus1_s, arma::cube sum_gammaklnyuline_42_m_nminus1_s, arma::cube sum_gammaklnyuline_51_m_nminus1_s, arma::cube sum_gammaklnyuline_52_m_nminus1_s) {
  // r>0 and s>0, the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_t_cube(m, dimdelta, dimdelta); ddlnhulinetilde_ddeltaddeltaprime_t_cube.fill(0.0); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::vec lambda_vec = kappa_vec.head(r); 
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
  arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
  arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat zerosmat_m_1(m, 1); zerosmat_m_1.fill(0.0);
  int Dim1 = m+r+2*s; int Dim2 = Dim1 + r*m*m; int Dim3 = Dim2 + s*m*m;
  // arma::vec lnyuline_tminus2(m);
  // if(t == 2) {
  //   lnyuline_tminus2 = zerosmat_m_1.col(0);
  // } else {
  //   arma::vec y_tminus2 = y_m_n.col(t-3); // \bm{y}_{t-2}
  //   lnyuline_tminus2 = log(pow(y_tminus2, 2)); // \ln\underline{\bm{y}}_{t-2} = (\ln{y_{1,t-2}^2}, ..., \ln{y_{m,t-2}^2})'
  // }
  arma::mat sum_lambdaklnyuline_1_tminus1_m_r = sum_lambdaklnyuline_1_m_nminus1_r.col(t-2);
  arma::mat sum_lambdaklnyuline_2_tminus1_m_r = sum_lambdaklnyuline_2_m_nminus1_r.col(t-2);
  for(int k = 0; k < r; k++) {
    // double lambda_k = lambda_vec(k);
    arma::vec g_0k = g_0_mat.col(k);
    arma::mat G_0k = fc_asmat(g_0k, m, m);
    // arma::vec sum_lambdaklnyuline_1 = lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_lambdaklnyuline_2 = zerosmat_m_1.col(0); // \sum_{i=3}^{t-1} (i-1) (i-2) \lambda_k^{i-3} \ln\underline{\bm{y}}_{t-i}
    // for(int i = 3; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi, 2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1 + (i-1) * pow(lambda_k, i-2) * lnyuline_tminusi;
    //   sum_lambdaklnyuline_2 = sum_lambdaklnyuline_2 + (i-1) * (i-2) * pow(lambda_k, i-3) * lnyuline_tminusi;
    // }
    arma::vec sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1_tminus1_m_r.col(k); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_lambdaklnyuline_2 = sum_lambdaklnyuline_2_tminus1_m_r.col(k); // \sum_{i=3}^{t-1} (i-1) (i-2) \lambda_k^{i-3} \ln\underline{\bm{y}}_{t-i}
    // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\bm\delta'
    arma::mat ddlnhulinetilde_dlambdakddeltaprime_t_mat(m, dimdelta); ddlnhulinetilde_dlambdakddeltaprime_t_mat.fill(0.0); // the \ell-th column is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\delta_\ell
    ddlnhulinetilde_dlambdakddeltaprime_t_mat.col(m+k) = G_0k * sum_lambdaklnyuline_2; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\lambda_k
    arma::mat sum_lambdaklnyuline_1_mat = fc_asmat(sum_lambdaklnyuline_1, m, 1);
    ddlnhulinetilde_dlambdakddeltaprime_t_mat.cols(Dim1 +k*m*m,   Dim1 +(k+1)*m*m-1) = kron(sum_lambdaklnyuline_1_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\bm{g}_{0,k}'
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(m+k) = ddlnhulinetilde_dlambdakddeltaprime_t_mat;
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.col(m+k) = ddlnhulinetilde_dlambdakddeltaprime_t_mat;
  }
  arma::mat sum_gammaklnyuline_11_tminus1_m_s = sum_gammaklnyuline_11_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_12_tminus1_m_s = sum_gammaklnyuline_12_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_21_tminus1_m_s = sum_gammaklnyuline_21_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_22_tminus1_m_s = sum_gammaklnyuline_22_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_31_tminus1_m_s = sum_gammaklnyuline_31_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_32_tminus1_m_s = sum_gammaklnyuline_32_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_41_tminus1_m_s = sum_gammaklnyuline_41_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_42_tminus1_m_s = sum_gammaklnyuline_42_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_51_tminus1_m_s = sum_gammaklnyuline_51_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_52_tminus1_m_s = sum_gammaklnyuline_52_m_nminus1_s.col(t-2);
  for(int k = 0; k < s; k++) {
    // double gamma_k = gamma_vec(k);
    // double varphi_k = varphi_vec(k);
    arma::vec g_1k = g_1_mat.col(k);
    arma::vec g_2k = g_2_mat.col(k);
    arma::mat G_1k = fc_asmat(g_1k, m, m);
    arma::mat G_2k = fc_asmat(g_2k, m, m);
    // arma::vec sum_gammaklnyuline_11 = cos(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_12 = sin(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_21 = gamma_k * (-sin(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_22 = gamma_k * cos(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_31 = (-sin(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_32 = cos(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_41 = gamma_k * (-cos(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\cos((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_42 = gamma_k * (-sin(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_51 = zerosmat_m_1.col(0); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_52 = zerosmat_m_1.col(0); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // for(int i = 3; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_gammaklnyuline_11 = sum_gammaklnyuline_11 + (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_12 = sum_gammaklnyuline_12 + (i-1) * pow(gamma_k, i-2) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_21 = sum_gammaklnyuline_21 + (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_22 = sum_gammaklnyuline_22 + (i-1) * pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_31 = sum_gammaklnyuline_31 + (i-1) * (i-1) * pow(gamma_k, i-2) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_32 = sum_gammaklnyuline_32 + (i-1) * (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_41 = sum_gammaklnyuline_41 + (i-1) * (i-1) * pow(gamma_k, i-1) * (-cos((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_42 = sum_gammaklnyuline_42 + (i-1) * (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_51 = sum_gammaklnyuline_51 + (i-1) * (i-2) * pow(gamma_k, i-3) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_52 = sum_gammaklnyuline_52 + (i-1) * (i-2) * pow(gamma_k, i-3) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    // }
    arma::vec sum_gammaklnyuline_11 = sum_gammaklnyuline_11_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_12 = sum_gammaklnyuline_12_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_21 = sum_gammaklnyuline_21_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_22 = sum_gammaklnyuline_22_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_31 = sum_gammaklnyuline_31_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_32 = sum_gammaklnyuline_32_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_41 = sum_gammaklnyuline_41_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\cos((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_42 = sum_gammaklnyuline_42_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_51 = sum_gammaklnyuline_51_tminus1_m_s.col(k); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_52 = sum_gammaklnyuline_52_tminus1_m_s.col(k); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\bm\delta'
    arma::mat ddlnhulinetilde_dgammakddeltaprime_t_mat(m, dimdelta); ddlnhulinetilde_dgammakddeltaprime_t_mat.fill(0.0); // the \ell-th column is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\delta_\ell
    ddlnhulinetilde_dgammakddeltaprime_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_51 + G_2k * sum_gammaklnyuline_52; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\gamma_k
    ddlnhulinetilde_dgammakddeltaprime_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_31 + G_2k * sum_gammaklnyuline_32; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\varphi_k
    arma::mat sum_gammaklnyuline_11_mat = fc_asmat(sum_gammaklnyuline_11, m, 1);
    arma::mat sum_gammaklnyuline_12_mat = fc_asmat(sum_gammaklnyuline_12, m, 1);
    ddlnhulinetilde_dgammakddeltaprime_t_mat.cols(Dim2 +k*m*m,     Dim2 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_11_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\bm{g}_{1,k}'
    ddlnhulinetilde_dgammakddeltaprime_t_mat.cols(Dim3 +k*m*m,     Dim3 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_12_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\bm{g}_{2,k}'
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(m+r+k) = ddlnhulinetilde_dgammakddeltaprime_t_mat;
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.col(m+r+k) = ddlnhulinetilde_dgammakddeltaprime_t_mat;
    // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\bm\delta'
    arma::mat ddlnhulinetilde_dvarphikddeltaprime_t_mat(m, dimdelta); ddlnhulinetilde_dvarphikddeltaprime_t_mat.fill(0.0); // the \ell-th column is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\delta_\ell
    ddlnhulinetilde_dvarphikddeltaprime_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_31 + G_2k * sum_gammaklnyuline_32; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\gamma_k
    ddlnhulinetilde_dvarphikddeltaprime_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_41 + G_2k * sum_gammaklnyuline_42; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\varphi_k
    arma::mat sum_gammaklnyuline_21_mat = fc_asmat(sum_gammaklnyuline_21, m, 1);
    arma::mat sum_gammaklnyuline_22_mat = fc_asmat(sum_gammaklnyuline_22, m, 1);
    ddlnhulinetilde_dvarphikddeltaprime_t_mat.cols(Dim2 +k*m*m,     Dim2 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_21_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\bm{g}_{1,k}'
    ddlnhulinetilde_dvarphikddeltaprime_t_mat.cols(Dim3 +k*m*m,     Dim3 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_22_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\bm{g}_{2,k}'
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(m+r+s+k) = ddlnhulinetilde_dvarphikddeltaprime_t_mat;
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.col(m+r+s+k) = ddlnhulinetilde_dvarphikddeltaprime_t_mat;
  }
  return ddlnhulinetilde_ddeltaddeltaprime_t_cube;
}

arma::cube fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube_requal0(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::cube sum_gammaklnyuline_11_m_nminus1_s, arma::cube sum_gammaklnyuline_12_m_nminus1_s, arma::cube sum_gammaklnyuline_21_m_nminus1_s, arma::cube sum_gammaklnyuline_22_m_nminus1_s, arma::cube sum_gammaklnyuline_31_m_nminus1_s, arma::cube sum_gammaklnyuline_32_m_nminus1_s, arma::cube sum_gammaklnyuline_41_m_nminus1_s, arma::cube sum_gammaklnyuline_42_m_nminus1_s, arma::cube sum_gammaklnyuline_51_m_nminus1_s, arma::cube sum_gammaklnyuline_52_m_nminus1_s) {
  // r=0 and s>0, the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_t_cube(m, dimdelta, dimdelta); ddlnhulinetilde_ddeltaddeltaprime_t_cube.fill(0.0); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
  // arma::vec lambda_vec = kappa_vec.head(r); 
  arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  // arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
  arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
  arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat zerosmat_m_1(m, 1); zerosmat_m_1.fill(0.0);
  int Dim1 = m+r+2*s; int Dim2 = Dim1 + r*m*m; int Dim3 = Dim2 + s*m*m;
  // arma::vec lnyuline_tminus2(m);
  // if(t == 2) {
  //   lnyuline_tminus2 = zerosmat_m_1.col(0);
  // } else {
  //   arma::vec y_tminus2 = y_m_n.col(t-3); // \bm{y}_{t-2}
  //   lnyuline_tminus2 = log(pow(y_tminus2, 2)); // \ln\underline{\bm{y}}_{t-2} = (\ln{y_{1,t-2}^2}, ..., \ln{y_{m,t-2}^2})'
  // }
  // arma::mat sum_lambdaklnyuline_1_tminus1_m_r = sum_lambdaklnyuline_1_m_nminus1_r.col(t-2);
  // arma::mat sum_lambdaklnyuline_2_tminus1_m_r = sum_lambdaklnyuline_2_m_nminus1_r.col(t-2);
  // for(int k = 0; k < r; k++) {
  //   // double lambda_k = lambda_vec(k);
  //   arma::vec g_0k = g_0_mat.col(k);
  //   arma::mat G_0k = fc_asmat(g_0k, m, m);
  //   // arma::vec sum_lambdaklnyuline_1 = lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_lambdaklnyuline_2 = zerosmat_m_1.col(0); // \sum_{i=3}^{t-1} (i-1) (i-2) \lambda_k^{i-3} \ln\underline{\bm{y}}_{t-i}
  //   // for(int i = 3; i < t; i++) {
  //   //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
  //   //   arma::vec lnyuline_tminusi = log(pow(y_tminusi, 2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
  //   //   sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1 + (i-1) * pow(lambda_k, i-2) * lnyuline_tminusi;
  //   //   sum_lambdaklnyuline_2 = sum_lambdaklnyuline_2 + (i-1) * (i-2) * pow(lambda_k, i-3) * lnyuline_tminusi;
  //   // }
  //   arma::vec sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1_tminus1_m_r.col(k); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_lambdaklnyuline_2 = sum_lambdaklnyuline_2_tminus1_m_r.col(k); // \sum_{i=3}^{t-1} (i-1) (i-2) \lambda_k^{i-3} \ln\underline{\bm{y}}_{t-i}
  //   // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\bm\delta'
  //   arma::mat ddlnhulinetilde_dlambdakddeltaprime_t_mat(m, dimdelta); ddlnhulinetilde_dlambdakddeltaprime_t_mat.fill(0.0); // the \ell-th column is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\delta_\ell
  //   ddlnhulinetilde_dlambdakddeltaprime_t_mat.col(m+k) = G_0k * sum_lambdaklnyuline_2; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\lambda_k
  //   arma::mat sum_lambdaklnyuline_1_mat = fc_asmat(sum_lambdaklnyuline_1, m, 1);
  //   ddlnhulinetilde_dlambdakddeltaprime_t_mat.cols(Dim1 +k*m*m,   Dim1 +(k+1)*m*m-1) = kron(sum_lambdaklnyuline_1_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\bm{g}_{0,k}'
  //   ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(m+k) = ddlnhulinetilde_dlambdakddeltaprime_t_mat;
  //   ddlnhulinetilde_ddeltaddeltaprime_t_cube.col(m+k) = ddlnhulinetilde_dlambdakddeltaprime_t_mat;
  // }
  arma::mat sum_gammaklnyuline_11_tminus1_m_s = sum_gammaklnyuline_11_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_12_tminus1_m_s = sum_gammaklnyuline_12_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_21_tminus1_m_s = sum_gammaklnyuline_21_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_22_tminus1_m_s = sum_gammaklnyuline_22_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_31_tminus1_m_s = sum_gammaklnyuline_31_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_32_tminus1_m_s = sum_gammaklnyuline_32_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_41_tminus1_m_s = sum_gammaklnyuline_41_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_42_tminus1_m_s = sum_gammaklnyuline_42_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_51_tminus1_m_s = sum_gammaklnyuline_51_m_nminus1_s.col(t-2);
  arma::mat sum_gammaklnyuline_52_tminus1_m_s = sum_gammaklnyuline_52_m_nminus1_s.col(t-2);
  for(int k = 0; k < s; k++) {
    // double gamma_k = gamma_vec(k);
    // double varphi_k = varphi_vec(k);
    arma::vec g_1k = g_1_mat.col(k);
    arma::vec g_2k = g_2_mat.col(k);
    arma::mat G_1k = fc_asmat(g_1k, m, m);
    arma::mat G_2k = fc_asmat(g_2k, m, m);
    // arma::vec sum_gammaklnyuline_11 = cos(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_12 = sin(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_21 = gamma_k * (-sin(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_22 = gamma_k * cos(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_31 = (-sin(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_32 = cos(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_41 = gamma_k * (-cos(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\cos((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_42 = gamma_k * (-sin(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_51 = zerosmat_m_1.col(0); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_gammaklnyuline_52 = zerosmat_m_1.col(0); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // for(int i = 3; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_gammaklnyuline_11 = sum_gammaklnyuline_11 + (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_12 = sum_gammaklnyuline_12 + (i-1) * pow(gamma_k, i-2) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_21 = sum_gammaklnyuline_21 + (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_22 = sum_gammaklnyuline_22 + (i-1) * pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_31 = sum_gammaklnyuline_31 + (i-1) * (i-1) * pow(gamma_k, i-2) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_32 = sum_gammaklnyuline_32 + (i-1) * (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_41 = sum_gammaklnyuline_41 + (i-1) * (i-1) * pow(gamma_k, i-1) * (-cos((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_42 = sum_gammaklnyuline_42 + (i-1) * (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
    //   sum_gammaklnyuline_51 = sum_gammaklnyuline_51 + (i-1) * (i-2) * pow(gamma_k, i-3) * cos((i-1) * varphi_k) * lnyuline_tminusi;
    //   sum_gammaklnyuline_52 = sum_gammaklnyuline_52 + (i-1) * (i-2) * pow(gamma_k, i-3) * sin((i-1) * varphi_k) * lnyuline_tminusi;
    // }
    arma::vec sum_gammaklnyuline_11 = sum_gammaklnyuline_11_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_12 = sum_gammaklnyuline_12_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_21 = sum_gammaklnyuline_21_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_22 = sum_gammaklnyuline_22_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_31 = sum_gammaklnyuline_31_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_32 = sum_gammaklnyuline_32_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_41 = sum_gammaklnyuline_41_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\cos((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_42 = sum_gammaklnyuline_42_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_51 = sum_gammaklnyuline_51_tminus1_m_s.col(k); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_gammaklnyuline_52 = sum_gammaklnyuline_52_tminus1_m_s.col(k); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
    // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\bm\delta'
    arma::mat ddlnhulinetilde_dgammakddeltaprime_t_mat(m, dimdelta); ddlnhulinetilde_dgammakddeltaprime_t_mat.fill(0.0); // the \ell-th column is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\delta_\ell
    ddlnhulinetilde_dgammakddeltaprime_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_51 + G_2k * sum_gammaklnyuline_52; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\gamma_k
    ddlnhulinetilde_dgammakddeltaprime_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_31 + G_2k * sum_gammaklnyuline_32; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\varphi_k
    arma::mat sum_gammaklnyuline_11_mat = fc_asmat(sum_gammaklnyuline_11, m, 1);
    arma::mat sum_gammaklnyuline_12_mat = fc_asmat(sum_gammaklnyuline_12, m, 1);
    ddlnhulinetilde_dgammakddeltaprime_t_mat.cols(Dim2 +k*m*m,     Dim2 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_11_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\bm{g}_{1,k}'
    ddlnhulinetilde_dgammakddeltaprime_t_mat.cols(Dim3 +k*m*m,     Dim3 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_12_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\bm{g}_{2,k}'
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(m+r+k) = ddlnhulinetilde_dgammakddeltaprime_t_mat;
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.col(m+r+k) = ddlnhulinetilde_dgammakddeltaprime_t_mat;
    // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\bm\delta'
    arma::mat ddlnhulinetilde_dvarphikddeltaprime_t_mat(m, dimdelta); ddlnhulinetilde_dvarphikddeltaprime_t_mat.fill(0.0); // the \ell-th column is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\delta_\ell
    ddlnhulinetilde_dvarphikddeltaprime_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_31 + G_2k * sum_gammaklnyuline_32; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\gamma_k
    ddlnhulinetilde_dvarphikddeltaprime_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_41 + G_2k * sum_gammaklnyuline_42; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\varphi_k
    arma::mat sum_gammaklnyuline_21_mat = fc_asmat(sum_gammaklnyuline_21, m, 1);
    arma::mat sum_gammaklnyuline_22_mat = fc_asmat(sum_gammaklnyuline_22, m, 1);
    ddlnhulinetilde_dvarphikddeltaprime_t_mat.cols(Dim2 +k*m*m,     Dim2 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_21_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\bm{g}_{1,k}'
    ddlnhulinetilde_dvarphikddeltaprime_t_mat.cols(Dim3 +k*m*m,     Dim3 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_22_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\bm{g}_{2,k}'
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(m+r+s+k) = ddlnhulinetilde_dvarphikddeltaprime_t_mat;
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.col(m+r+s+k) = ddlnhulinetilde_dvarphikddeltaprime_t_mat;
  }
  return ddlnhulinetilde_ddeltaddeltaprime_t_cube;
}

arma::cube fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube_sequal0(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_1_m_nminus1_r, arma::cube sum_lambdaklnyuline_2_m_nminus1_r) {
  // r>0 and s=0, the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_t_cube(m, dimdelta, dimdelta); ddlnhulinetilde_ddeltaddeltaprime_t_cube.fill(0.0); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::vec lambda_vec = kappa_vec.head(r); 
  // arma::vec gamma_vec = kappa_vec.subvec(r, r+s-1);
  // arma::vec varphi_vec = kappa_vec.subvec(r+s, r+2*s-1);
  arma::vec g_0_vec = kappa_vec.subvec(r+2*s, r+2*s+r*m*m-1); arma::mat g_0_mat = fc_asmat(g_0_vec, m*m, r);
  // arma::vec g_1_vec = kappa_vec.subvec(r+2*s+r*m*m, r+2*s+r*m*m+s*m*m-1); arma::mat g_1_mat = fc_asmat(g_1_vec, m*m, s);
  // arma::vec g_2_vec = kappa_vec.tail(s*m*m); arma::mat g_2_mat = fc_asmat(g_2_vec, m*m, s);
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat zerosmat_m_1(m, 1); zerosmat_m_1.fill(0.0);
  int Dim1 = m+r+2*s; // int Dim2 = Dim1 + r*m*m; int Dim3 = Dim2 + s*m*m;
  // arma::vec lnyuline_tminus2(m);
  // if(t == 2) {
  //   lnyuline_tminus2 = zerosmat_m_1.col(0);
  // } else {
  //   arma::vec y_tminus2 = y_m_n.col(t-3); // \bm{y}_{t-2}
  //   lnyuline_tminus2 = log(pow(y_tminus2, 2)); // \ln\underline{\bm{y}}_{t-2} = (\ln{y_{1,t-2}^2}, ..., \ln{y_{m,t-2}^2})'
  // }
  arma::mat sum_lambdaklnyuline_1_tminus1_m_r = sum_lambdaklnyuline_1_m_nminus1_r.col(t-2);
  arma::mat sum_lambdaklnyuline_2_tminus1_m_r = sum_lambdaklnyuline_2_m_nminus1_r.col(t-2);
  for(int k = 0; k < r; k++) {
    // double lambda_k = lambda_vec(k);
    arma::vec g_0k = g_0_mat.col(k);
    arma::mat G_0k = fc_asmat(g_0k, m, m);
    // arma::vec sum_lambdaklnyuline_1 = lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    // arma::vec sum_lambdaklnyuline_2 = zerosmat_m_1.col(0); // \sum_{i=3}^{t-1} (i-1) (i-2) \lambda_k^{i-3} \ln\underline{\bm{y}}_{t-i}
    // for(int i = 3; i < t; i++) {
    //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
    //   arma::vec lnyuline_tminusi = log(pow(y_tminusi, 2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
    //   sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1 + (i-1) * pow(lambda_k, i-2) * lnyuline_tminusi;
    //   sum_lambdaklnyuline_2 = sum_lambdaklnyuline_2 + (i-1) * (i-2) * pow(lambda_k, i-3) * lnyuline_tminusi;
    // }
    arma::vec sum_lambdaklnyuline_1 = sum_lambdaklnyuline_1_tminus1_m_r.col(k); // \sum_{i=2}^{t-1} (i-1) \lambda_k^{i-2} \ln\underline{\bm{y}}_{t-i}
    arma::vec sum_lambdaklnyuline_2 = sum_lambdaklnyuline_2_tminus1_m_r.col(k); // \sum_{i=3}^{t-1} (i-1) (i-2) \lambda_k^{i-3} \ln\underline{\bm{y}}_{t-i}
    // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\bm\delta'
    arma::mat ddlnhulinetilde_dlambdakddeltaprime_t_mat(m, dimdelta); ddlnhulinetilde_dlambdakddeltaprime_t_mat.fill(0.0); // the \ell-th column is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\delta_\ell
    ddlnhulinetilde_dlambdakddeltaprime_t_mat.col(m+k) = G_0k * sum_lambdaklnyuline_2; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\lambda_k
    arma::mat sum_lambdaklnyuline_1_mat = fc_asmat(sum_lambdaklnyuline_1, m, 1);
    ddlnhulinetilde_dlambdakddeltaprime_t_mat.cols(Dim1 +k*m*m,   Dim1 +(k+1)*m*m-1) = kron(sum_lambdaklnyuline_1_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\lambda_k\partial\bm{g}_{0,k}'
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(m+k) = ddlnhulinetilde_dlambdakddeltaprime_t_mat;
    ddlnhulinetilde_ddeltaddeltaprime_t_cube.col(m+k) = ddlnhulinetilde_dlambdakddeltaprime_t_mat;
  }
  // arma::mat sum_gammaklnyuline_11_tminus1_m_s = sum_gammaklnyuline_11_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_12_tminus1_m_s = sum_gammaklnyuline_12_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_21_tminus1_m_s = sum_gammaklnyuline_21_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_22_tminus1_m_s = sum_gammaklnyuline_22_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_31_tminus1_m_s = sum_gammaklnyuline_31_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_32_tminus1_m_s = sum_gammaklnyuline_32_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_41_tminus1_m_s = sum_gammaklnyuline_41_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_42_tminus1_m_s = sum_gammaklnyuline_42_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_51_tminus1_m_s = sum_gammaklnyuline_51_m_nminus1_s.col(t-2);
  // arma::mat sum_gammaklnyuline_52_tminus1_m_s = sum_gammaklnyuline_52_m_nminus1_s.col(t-2);
  // for(int k = 0; k < s; k++) {
  //   // double gamma_k = gamma_vec(k);
  //   // double varphi_k = varphi_vec(k);
  //   arma::vec g_1k = g_1_mat.col(k);
  //   arma::vec g_2k = g_2_mat.col(k);
  //   arma::mat G_1k = fc_asmat(g_1k, m, m);
  //   arma::mat G_2k = fc_asmat(g_2k, m, m);
  //   // arma::vec sum_gammaklnyuline_11 = cos(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_12 = sin(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_21 = gamma_k * (-sin(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_22 = gamma_k * cos(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_31 = (-sin(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_32 = cos(varphi_k) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_41 = gamma_k * (-cos(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\cos((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_42 = gamma_k * (-sin(varphi_k)) * lnyuline_tminus2; // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_51 = zerosmat_m_1.col(0); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // arma::vec sum_gammaklnyuline_52 = zerosmat_m_1.col(0); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // for(int i = 3; i < t; i++) {
  //   //   arma::vec y_tminusi = y_m_n.col(t-i-1); // \bm{y}_{t-i}
  //   //   arma::vec lnyuline_tminusi = log(pow(y_tminusi,2)); // \ln\underline{\bm{y}}_{t-i} = (\ln{y_{1,t-i}^2}, ..., \ln{y_{m,t-i}^2})'
  //   //   sum_gammaklnyuline_11 = sum_gammaklnyuline_11 + (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_12 = sum_gammaklnyuline_12 + (i-1) * pow(gamma_k, i-2) * sin((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_21 = sum_gammaklnyuline_21 + (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_22 = sum_gammaklnyuline_22 + (i-1) * pow(gamma_k, i-1) * cos((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_31 = sum_gammaklnyuline_31 + (i-1) * (i-1) * pow(gamma_k, i-2) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_32 = sum_gammaklnyuline_32 + (i-1) * (i-1) * pow(gamma_k, i-2) * cos((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_41 = sum_gammaklnyuline_41 + (i-1) * (i-1) * pow(gamma_k, i-1) * (-cos((i-1) * varphi_k)) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_42 = sum_gammaklnyuline_42 + (i-1) * (i-1) * pow(gamma_k, i-1) * (-sin((i-1) * varphi_k)) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_51 = sum_gammaklnyuline_51 + (i-1) * (i-2) * pow(gamma_k, i-3) * cos((i-1) * varphi_k) * lnyuline_tminusi;
  //   //   sum_gammaklnyuline_52 = sum_gammaklnyuline_52 + (i-1) * (i-2) * pow(gamma_k, i-3) * sin((i-1) * varphi_k) * lnyuline_tminusi;
  //   // }
  //   arma::vec sum_gammaklnyuline_11 = sum_gammaklnyuline_11_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_12 = sum_gammaklnyuline_12_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-2} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_21 = sum_gammaklnyuline_21_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_22 = sum_gammaklnyuline_22_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1) \gamma_k^{i-1} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_31 = sum_gammaklnyuline_31_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_32 = sum_gammaklnyuline_32_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-2} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_41 = sum_gammaklnyuline_41_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\cos((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_42 = sum_gammaklnyuline_42_tminus1_m_s.col(k); // \sum_{i=2}^{t-1} (i-1)^2 \gamma_k^{i-1} (-\sin((i-1)\varphi_k)) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_51 = sum_gammaklnyuline_51_tminus1_m_s.col(k); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \cos((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   arma::vec sum_gammaklnyuline_52 = sum_gammaklnyuline_52_tminus1_m_s.col(k); // \sum_{i=3}^{t-1} (i-1) (i-2) \gamma_k^{i-3} \sin((i-1)\varphi_k) \ln\underline{\bm{y}}_{t-i}
  //   // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\bm\delta'
  //   arma::mat ddlnhulinetilde_dgammakddeltaprime_t_mat(m, dimdelta); ddlnhulinetilde_dgammakddeltaprime_t_mat.fill(0.0); // the \ell-th column is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\delta_\ell
  //   ddlnhulinetilde_dgammakddeltaprime_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_51 + G_2k * sum_gammaklnyuline_52; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\gamma_k
  //   ddlnhulinetilde_dgammakddeltaprime_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_31 + G_2k * sum_gammaklnyuline_32; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\varphi_k
  //   arma::mat sum_gammaklnyuline_11_mat = fc_asmat(sum_gammaklnyuline_11, m, 1);
  //   arma::mat sum_gammaklnyuline_12_mat = fc_asmat(sum_gammaklnyuline_12, m, 1);
  //   ddlnhulinetilde_dgammakddeltaprime_t_mat.cols(Dim2 +k*m*m,     Dim2 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_11_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\bm{g}_{1,k}'
  //   ddlnhulinetilde_dgammakddeltaprime_t_mat.cols(Dim3 +k*m*m,     Dim3 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_12_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\gamma_k\partial\bm{g}_{2,k}'
  //   ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(m+r+k) = ddlnhulinetilde_dgammakddeltaprime_t_mat;
  //   ddlnhulinetilde_ddeltaddeltaprime_t_cube.col(m+r+k) = ddlnhulinetilde_dgammakddeltaprime_t_mat;
  //   // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\bm\delta'
  //   arma::mat ddlnhulinetilde_dvarphikddeltaprime_t_mat(m, dimdelta); ddlnhulinetilde_dvarphikddeltaprime_t_mat.fill(0.0); // the \ell-th column is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\delta_\ell
  //   ddlnhulinetilde_dvarphikddeltaprime_t_mat.col(m+r+k) = G_1k * sum_gammaklnyuline_31 + G_2k * sum_gammaklnyuline_32; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\gamma_k
  //   ddlnhulinetilde_dvarphikddeltaprime_t_mat.col(m+r+s+k) = G_1k * sum_gammaklnyuline_41 + G_2k * sum_gammaklnyuline_42; // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\varphi_k
  //   arma::mat sum_gammaklnyuline_21_mat = fc_asmat(sum_gammaklnyuline_21, m, 1);
  //   arma::mat sum_gammaklnyuline_22_mat = fc_asmat(sum_gammaklnyuline_22, m, 1);
  //   ddlnhulinetilde_dvarphikddeltaprime_t_mat.cols(Dim2 +k*m*m,     Dim2 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_21_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\bm{g}_{1,k}'
  //   ddlnhulinetilde_dvarphikddeltaprime_t_mat.cols(Dim3 +k*m*m,     Dim3 +(k+1)*m*m-1)   = kron(sum_gammaklnyuline_22_mat.t(), I_m); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\varphi_k\partial\bm{g}_{2,k}'
  //   ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(m+r+s+k) = ddlnhulinetilde_dvarphikddeltaprime_t_mat;
  //   ddlnhulinetilde_ddeltaddeltaprime_t_cube.col(m+r+s+k) = ddlnhulinetilde_dvarphikddeltaprime_t_mat;
  // }
  return ddlnhulinetilde_ddeltaddeltaprime_t_cube;
}

arma::cube fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube_general(int t, int m, int r, int s, int dimdelta, arma::vec kappa_vec, arma::cube sum_lambdaklnyuline_1_m_nminus1_r, arma::cube sum_lambdaklnyuline_2_m_nminus1_r, arma::cube sum_gammaklnyuline_11_m_nminus1_s, arma::cube sum_gammaklnyuline_12_m_nminus1_s, arma::cube sum_gammaklnyuline_21_m_nminus1_s, arma::cube sum_gammaklnyuline_22_m_nminus1_s, arma::cube sum_gammaklnyuline_31_m_nminus1_s, arma::cube sum_gammaklnyuline_32_m_nminus1_s, arma::cube sum_gammaklnyuline_41_m_nminus1_s, arma::cube sum_gammaklnyuline_42_m_nminus1_s, arma::cube sum_gammaklnyuline_51_m_nminus1_s, arma::cube sum_gammaklnyuline_52_m_nminus1_s) {
  // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_t_cube(m, dimdelta, dimdelta); ddlnhulinetilde_ddeltaddeltaprime_t_cube.fill(0.0); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
  if(s == 0) {
    ddlnhulinetilde_ddeltaddeltaprime_t_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube_sequal0(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_1_m_nminus1_r, sum_lambdaklnyuline_2_m_nminus1_r);
  } else if(r == 0) {
    ddlnhulinetilde_ddeltaddeltaprime_t_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube_requal0(t, m, r, s, dimdelta, kappa_vec, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s, sum_gammaklnyuline_31_m_nminus1_s, sum_gammaklnyuline_32_m_nminus1_s, sum_gammaklnyuline_41_m_nminus1_s, sum_gammaklnyuline_42_m_nminus1_s, sum_gammaklnyuline_51_m_nminus1_s, sum_gammaklnyuline_52_m_nminus1_s);
  } else {
    ddlnhulinetilde_ddeltaddeltaprime_t_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_1_m_nminus1_r, sum_lambdaklnyuline_2_m_nminus1_r, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s, sum_gammaklnyuline_31_m_nminus1_s, sum_gammaklnyuline_32_m_nminus1_s, sum_gammaklnyuline_41_m_nminus1_s, sum_gammaklnyuline_42_m_nminus1_s, sum_gammaklnyuline_51_m_nminus1_s, sum_gammaklnyuline_52_m_nminus1_s);
  }
  return ddlnhulinetilde_ddeltaddeltaprime_t_cube;
}

arma::mat fc_ddDtilde_ddeltaddeltaprime_t_kl(arma::vec dlnhulinetilde_ddelta_t_l, arma::vec dlnhulinetilde_ddelta_t_k, arma::vec ddlnhulinetilde_ddeltaddeltaprime_t_kl, arma::mat Dtilde_t) {
  // \partial^2\widetilde{D}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::vec ddDtilde_ddeltaddeltaprime_t_kl_diag = 0.25 * Dtilde_t.diag(0) % dlnhulinetilde_ddelta_t_k % dlnhulinetilde_ddelta_t_l + 0.5 * Dtilde_t.diag(0) % ddlnhulinetilde_ddeltaddeltaprime_t_kl;
  arma::mat ddDtilde_ddeltaddeltaprime_t_kl = diagmat(ddDtilde_ddeltaddeltaprime_t_kl_diag, 0);
  return ddDtilde_ddeltaddeltaprime_t_kl;
}

arma::mat fc_dvarepsilontilde_ddelta_tminusBkTOtminus1_l_mat(arma::mat varepsilontilde_tminusBkTOtminus1_mat, arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat) {
  // the k-th column is \partial\widetilde\varepsilon_{t-(\Bbbk+1-k)}/\partial\delta_ell
  arma::mat dvarepsilontilde_ddelta_tminusBkTOtminus1_l_mat = -0.5 * dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat % varepsilontilde_tminusBkTOtminus1_mat; 
  return dvarepsilontilde_ddelta_tminusBkTOtminus1_l_mat;
}

arma::mat fc_ddvarepsilontilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat(arma::mat varepsilontilde_tminusBkTOtminus1_mat, arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat, arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat, arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat) {
  // the i-th column is \partial^2\widetilde\varepsilon_{t-(\Bbbk+1-i)}/\partial\delta_k\partial\delta_ell
  arma::mat ddvarepsilontilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat = 0.25 * dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat % dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat % varepsilontilde_tminusBkTOtminus1_mat - 0.5 * ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat % varepsilontilde_tminusBkTOtminus1_mat;
  return ddvarepsilontilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat;
}

arma::mat fc_ddPsitilde_ddeltaddeltaprime_tminus1_kl(int m, arma::mat Psitilde_tminus1, arma::mat dPsitilde_ddelta_tminus1_k, arma::mat varepsilontilde_tminusBkTOtminus1_mat, arma::mat dvarepsilontilde_ddelta_tminusBkTOtminus1_l_mat, arma::mat dvarepsilontilde_ddelta_tminusBkTOtminus1_k_mat, arma::mat ddvarepsilontilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat) {
  // \partial^2\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::mat ddPsitilde_ddeltaddeltaprime_tminus1_kl(m, m); ddPsitilde_ddeltaddeltaprime_tminus1_kl.fill(0.0);
  for(int i = 1; i < m; i++) { // {2,1}; {3,1}, {3,2}; ...; {m,1}, ..., {m,m-1}
    for(int j = 0; j < i; j++) {
      double Psitilde_tminus1_ij = Psitilde_tminus1(i,j);
      double dPsitilde_ddelta_tminus1_k_ij = dPsitilde_ddelta_tminus1_k(i,j);
      arma::rowvec varepsilon_ith_Bk = varepsilontilde_tminusBkTOtminus1_mat.row(i); // (\widetilde\varepsilon_{i,t-\Bbbk}(\bm\delta), ..., \widetilde\varepsilon_{i,t-1}(\bm\delta))'
      arma::rowvec varepsilon_jth_Bk = varepsilontilde_tminusBkTOtminus1_mat.row(j); 
      arma::rowvec dvarepsilon_ddelta_l_ith_Bk = dvarepsilontilde_ddelta_tminusBkTOtminus1_l_mat.row(i); // (\partial\widetilde\varepsilon_{i,t-\Bbbk}(\bm\delta)/\partial\delta_\ell, ..., \partial\widetilde\varepsilon_{i,t-1}(\bm\delta)/\partial\delta_\ell)'
      arma::rowvec dvarepsilon_ddelta_l_jth_Bk = dvarepsilontilde_ddelta_tminusBkTOtminus1_l_mat.row(j);
      arma::rowvec dvarepsilon_ddelta_k_ith_Bk = dvarepsilontilde_ddelta_tminusBkTOtminus1_k_mat.row(i);
      arma::rowvec dvarepsilon_ddelta_k_jth_Bk = dvarepsilontilde_ddelta_tminusBkTOtminus1_k_mat.row(j);
      arma::rowvec ddvarepsilon_ddelta_kl_ith_Bk = ddvarepsilontilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.row(i); // (\partial^2\widetilde\varepsilon_{i,t-\Bbbk}(\bm\delta)/\partial\delta_k\partial\delta_\ell, ..., \partial^2\widetilde\varepsilon_{i,t-1}(\bm\delta)/\partial\delta_k\partial\delta_\ell)'
      arma::rowvec ddvarepsilon_ddelta_kl_jth_Bk = ddvarepsilontilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.row(j); // (\partial^2\widetilde\varepsilon_{i,t-\Bbbk}(\bm\delta)/\partial\delta_k\partial\delta_\ell, ..., \partial^2\widetilde\varepsilon_{i,t-1}(\bm\delta)/\partial\delta_k\partial\delta_\ell)'
      double d0d1_i_il = dot(varepsilon_ith_Bk, dvarepsilon_ddelta_l_ith_Bk) / dot(varepsilon_ith_Bk, varepsilon_ith_Bk);
      double d0d1_i_ik = dot(varepsilon_ith_Bk, dvarepsilon_ddelta_k_ith_Bk) / dot(varepsilon_ith_Bk, varepsilon_ith_Bk);
      double d0d1_j_jl = dot(varepsilon_jth_Bk, dvarepsilon_ddelta_l_jth_Bk) / dot(varepsilon_jth_Bk, varepsilon_jth_Bk);
      double d0d1_j_jk = dot(varepsilon_jth_Bk, dvarepsilon_ddelta_k_jth_Bk) / dot(varepsilon_jth_Bk, varepsilon_jth_Bk);
      double d0d1_i_jl = dot(varepsilon_ith_Bk, dvarepsilon_ddelta_l_jth_Bk) / sqrt(dot(varepsilon_ith_Bk, varepsilon_ith_Bk) * dot(varepsilon_jth_Bk, varepsilon_jth_Bk));
      double d0d1_il_j = dot(dvarepsilon_ddelta_l_ith_Bk, varepsilon_jth_Bk) / sqrt(dot(varepsilon_ith_Bk, varepsilon_ith_Bk) * dot(varepsilon_jth_Bk, varepsilon_jth_Bk));
      double d1d1_ik_il = dot(dvarepsilon_ddelta_k_ith_Bk, dvarepsilon_ddelta_l_ith_Bk) / dot(varepsilon_ith_Bk, varepsilon_ith_Bk);
      double d1d1_jk_jl = dot(dvarepsilon_ddelta_k_jth_Bk, dvarepsilon_ddelta_l_jth_Bk) / dot(varepsilon_jth_Bk, varepsilon_jth_Bk);
      double d1d1_il_jk = dot(dvarepsilon_ddelta_l_ith_Bk, dvarepsilon_ddelta_k_jth_Bk) / sqrt(dot(varepsilon_ith_Bk, varepsilon_ith_Bk) * dot(varepsilon_jth_Bk, varepsilon_jth_Bk));
      double d1d1_ik_jl = dot(dvarepsilon_ddelta_k_ith_Bk, dvarepsilon_ddelta_l_jth_Bk) / sqrt(dot(varepsilon_ith_Bk, varepsilon_ith_Bk) * dot(varepsilon_jth_Bk, varepsilon_jth_Bk));
      double d0d2_i_ikl = dot(varepsilon_ith_Bk, ddvarepsilon_ddelta_kl_ith_Bk) / dot(varepsilon_ith_Bk, varepsilon_ith_Bk);
      double d0d2_j_jkl = dot(varepsilon_jth_Bk, ddvarepsilon_ddelta_kl_jth_Bk) / dot(varepsilon_jth_Bk, varepsilon_jth_Bk);
      double d0d2_ikl_j = dot(ddvarepsilon_ddelta_kl_ith_Bk, varepsilon_jth_Bk) / sqrt(dot(varepsilon_ith_Bk, varepsilon_ith_Bk) * dot(varepsilon_jth_Bk, varepsilon_jth_Bk));
      double d0d2_i_jkl = dot(varepsilon_ith_Bk, ddvarepsilon_ddelta_kl_jth_Bk) / sqrt(dot(varepsilon_ith_Bk, varepsilon_ith_Bk) * dot(varepsilon_jth_Bk, varepsilon_jth_Bk));
      ddPsitilde_ddeltaddeltaprime_tminus1_kl(i,j) = ddPsitilde_ddeltaddeltaprime_tminus1_kl(j,i) = (d1d1_il_jk + d1d1_ik_jl + d0d2_ikl_j + d0d2_i_jkl) - 
                                                                                                    (d0d1_il_j + d0d1_i_jl) * (d0d1_i_ik + d0d1_j_jk) - 
                                                                                                    dPsitilde_ddelta_tminus1_k_ij * (d0d1_i_il + d0d1_j_jl) - 
                                                                                                    Psitilde_tminus1_ij * (d1d1_ik_il + d0d2_i_ikl + d1d1_jk_jl + d0d2_j_jkl) + 
                                                                                                    2.0 * Psitilde_tminus1_ij * (d0d1_i_ik * d0d1_i_il + d0d1_j_jk * d0d1_j_jl);
    }
  }
  return ddPsitilde_ddeltaddeltaprime_tminus1_kl;
}

arma::field<arma::mat> fc_ddRtilde_ddeltaddeltaprime_1_field(int m, int dimdelta) {
  // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\delta_k\partial\delta_\ell = 0_m
  arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_1_field(dimdelta, dimdelta);
  arma::mat zerosmat_m_m(m, m); zerosmat_m_m.fill(0.0);
  for(int k = 0; k < dimdelta; k++) {
    for(int l = 0; l < (k+1); l++) {
      ddRtilde_ddeltaddeltaprime_1_field(k,l) = ddRtilde_ddeltaddeltaprime_1_field(l,k) = zerosmat_m_m;
    }
  }
  return ddRtilde_ddeltaddeltaprime_1_field;
}

arma::mat fc_ddRtilde_ddeltaddeltaprime_t_kl(int m, double beta_1, double beta_2, arma::mat varepsilontilde_tminusBkTOtminus1_mat, arma::mat Psitilde_tminus1, arma::mat dPsitilde_ddelta_tminus1_k, arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat, arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat, arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat, arma::mat ddRtilde_ddeltaddeltaprime_tminus1_kl) {
  // \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
  arma::mat dvarepsilontilde_ddelta_tminusBkTOtminus1_l_mat = fc_dvarepsilontilde_ddelta_tminusBkTOtminus1_l_mat(varepsilontilde_tminusBkTOtminus1_mat, dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat); // (\partial\widetilde\varepsilon_{t-\Bbbk}/\partial\delta_ell, ..., \partial\widetilde\varepsilon_{t-1}/\partial\delta_ell)
  arma::mat dvarepsilontilde_ddelta_tminusBkTOtminus1_k_mat = fc_dvarepsilontilde_ddelta_tminusBkTOtminus1_l_mat(varepsilontilde_tminusBkTOtminus1_mat, dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat); // (\partial\widetilde\varepsilon_{t-\Bbbk}/\partial\delta_k, ..., \partial\widetilde\varepsilon_{t-1}/\partial\delta_k)
  arma::mat ddvarepsilontilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat = fc_ddvarepsilontilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat(varepsilontilde_tminusBkTOtminus1_mat, dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat, dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat, ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat); // (\partial^2\widetilde\varepsilon_{t-\Bbbk}/\partial\delta_k\partial\delta_ell, ..., \partial^2\widetilde\varepsilon_{t-1}/\partial\delta_k\partial\delta_ell)
  arma::mat ddPsitilde_ddeltaddeltaprime_tminus1_kl = fc_ddPsitilde_ddeltaddeltaprime_tminus1_kl(m, Psitilde_tminus1, dPsitilde_ddelta_tminus1_k, varepsilontilde_tminusBkTOtminus1_mat, dvarepsilontilde_ddelta_tminusBkTOtminus1_l_mat, dvarepsilontilde_ddelta_tminusBkTOtminus1_k_mat, ddvarepsilontilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat); // \partial^2\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::mat ddRtilde_ddeltaddeltaprime_t_kl = beta_1 * ddPsitilde_ddeltaddeltaprime_tminus1_kl + beta_2 * ddRtilde_ddeltaddeltaprime_tminus1_kl;
  return ddRtilde_ddeltaddeltaprime_t_kl;
}

arma::field<arma::mat> fc_ddRtilde_dbetadbetaprime_1_field(int m, int dimbeta, double beta_1, double beta_2, arma::mat Ruline, arma::cube dRuline_dbeta_cube) {
  // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::field<arma::mat> ddRtilde_dbetadbetaprime_1_field(dimbeta, dimbeta);
  arma::mat zerosmat_m_m(m, m); zerosmat_m_m.fill(0.0);
  arma::mat onesmat_m_m(m, m); onesmat_m_m.fill(1.0);
  ddRtilde_dbetadbetaprime_1_field(0,0) = zerosmat_m_m;
  ddRtilde_dbetadbetaprime_1_field(1,1) = - 2.0 * beta_1 / pow(1.0 - beta_2, 3) * (Ruline - onesmat_m_m);
  ddRtilde_dbetadbetaprime_1_field(1,0) = ddRtilde_dbetadbetaprime_1_field(0,1) = - 1.0 / pow(1 - beta_2, 2) * (Ruline - onesmat_m_m);
  for(int k = 2; k < dimbeta; k++) {
    ddRtilde_dbetadbetaprime_1_field(k,0) = ddRtilde_dbetadbetaprime_1_field(0,k) = - 1.0 / (1.0 - beta_2) * dRuline_dbeta_cube.slice(k);
    ddRtilde_dbetadbetaprime_1_field(k,1) = ddRtilde_dbetadbetaprime_1_field(1,k) = - beta_1 / pow(1.0 - beta_2, 2) * dRuline_dbeta_cube.slice(k);
    for(int l = 2; l < (k+1); l++) {
      ddRtilde_dbetadbetaprime_1_field(k,l) = ddRtilde_dbetadbetaprime_1_field(l,k) = zerosmat_m_m;
    }
  }
  return ddRtilde_dbetadbetaprime_1_field;
}

arma::mat fc_ddRtilde_dbetadbetaprime_t_kl(int l, int k, int m, double beta_2, arma::mat dRtilde_dbeta_tminus1_l, arma::mat ddRtilde_dbetadbetaprime_tminus1_kl) {
  // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::mat ddRtilde_dbetadbetaprime_t_kl(m, m); ddRtilde_dbetadbetaprime_t_kl.fill(0.0);
  arma::mat zerosmat_m_m(m, m); zerosmat_m_m.fill(0.0);
  if(k == 0 && l == 0) {
    ddRtilde_dbetadbetaprime_t_kl = zerosmat_m_m;
  } else if(k == 1 && l == 0) {
    ddRtilde_dbetadbetaprime_t_kl = dRtilde_dbeta_tminus1_l + beta_2 * ddRtilde_dbetadbetaprime_tminus1_kl;
  } else if(k == 1 && l == 1) {
    ddRtilde_dbetadbetaprime_t_kl = 2.0 * dRtilde_dbeta_tminus1_l + beta_2 * ddRtilde_dbetadbetaprime_tminus1_kl;
  } else if(k > 1 && l < 2) {
    ddRtilde_dbetadbetaprime_t_kl = ddRtilde_dbetadbetaprime_tminus1_kl;
  } else {
    ddRtilde_dbetadbetaprime_t_kl = zerosmat_m_m;
  }
  return ddRtilde_dbetadbetaprime_t_kl;
}

arma::field<arma::mat> fc_ddRtilde_dbetaddeltaprime_1_field(int m, int dimdelta, int dimbeta) {
  // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\beta_k\partial\delta_\ell = 0_m
  arma::field<arma::mat> ddRtilde_dbetaddeltaprime_1_field(dimbeta, dimdelta);
  arma::mat zerosmat_m_m(m, m); zerosmat_m_m.fill(0.0);
  for(int k = 0; k < dimbeta; k++) {
    for(int l = 0; l < dimdelta; l++) {
      ddRtilde_dbetaddeltaprime_1_field(k,l) = zerosmat_m_m;
    }
  }
  return ddRtilde_dbetaddeltaprime_1_field;
}

arma::mat fc_ddRtilde_dbetaddeltaprime_t_kl(int l, int k, int m, double beta_2, arma::mat dPsitilde_ddelta_tminus1_l, arma::mat dRtilde_ddelta_tminus1_l, arma::mat ddRtilde_dbetaddeltaprime_tminus1_kl) {
  // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
  arma::mat ddRtilde_dbetaddeltaprime_t_kl(m, m); ddRtilde_dbetaddeltaprime_t_kl.fill(0.0);
  arma::mat zerosmat_m_m(m, m); zerosmat_m_m.fill(0.0);
  if(k == 0) {
    ddRtilde_dbetaddeltaprime_t_kl = dPsitilde_ddelta_tminus1_l + beta_2 * ddRtilde_dbetaddeltaprime_tminus1_kl;
  } else if(k == 1) {
    ddRtilde_dbetaddeltaprime_t_kl = dRtilde_ddelta_tminus1_l + beta_2 * ddRtilde_dbetaddeltaprime_tminus1_kl;
  } else {
    ddRtilde_dbetaddeltaprime_t_kl = zerosmat_m_m;
  }
  return ddRtilde_dbetaddeltaprime_t_kl;
}

arma::mat fc_ddHtilde_ddeltaddeltaprime_t_kl(arma::mat Dtilde_t, arma::mat Rtilde_t, arma::mat dDtilde_ddelta_t_l, arma::mat dDtilde_ddelta_t_k, arma::mat dRtilde_ddelta_t_l, arma::mat dRtilde_ddelta_t_k, arma::mat ddDtilde_ddeltaddeltaprime_t_kl, arma::mat ddRtilde_ddeltaddeltaprime_t_kl) {
  // \partial^2\widetilde{H}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
  arma::mat ddHtilde_ddeltaddeltaprime_t_kl = ddDtilde_ddeltaddeltaprime_t_kl * Rtilde_t * Dtilde_t +
                                              dDtilde_ddelta_t_l * dRtilde_ddelta_t_k * Dtilde_t +
                                              dDtilde_ddelta_t_l * Rtilde_t * dDtilde_ddelta_t_k +
                                              dDtilde_ddelta_t_k * dRtilde_ddelta_t_l * Dtilde_t +
                                              Dtilde_t * ddRtilde_ddeltaddeltaprime_t_kl * Dtilde_t +
                                              Dtilde_t * dRtilde_ddelta_t_l * dDtilde_ddelta_t_k +
                                              dDtilde_ddelta_t_k * Rtilde_t * dDtilde_ddelta_t_l +
                                              Dtilde_t * dRtilde_ddelta_t_k * dDtilde_ddelta_t_l +
                                              Dtilde_t * Rtilde_t * ddDtilde_ddeltaddeltaprime_t_kl;
  return ddHtilde_ddeltaddeltaprime_t_kl;
}

arma::mat fc_ddHtilde_dbetadbetaprime_t_kl(arma::mat Dtilde_t, arma::mat ddRtilde_dbetadbetaprime_t_kl) {
  // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::mat ddHtilde_dbetadbetaprime_t_kl = Dtilde_t * ddRtilde_dbetadbetaprime_t_kl * Dtilde_t;
  return ddHtilde_dbetadbetaprime_t_kl;
}

arma::mat fc_ddHtilde_dbetaddeltaprime_t_kl(arma::mat Dtilde_t, arma::mat dDtilde_ddelta_t_l, arma::mat dRtilde_dbeta_t_k, arma::mat ddRtilde_dbetaddeltaprime_t_kl) {
  // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
  arma::mat ddHtilde_dbetaddeltaprime_t_kl = dDtilde_ddelta_t_l * dRtilde_dbeta_t_k * Dtilde_t +
                                              Dtilde_t * ddRtilde_dbetaddeltaprime_t_kl * Dtilde_t +
                                              Dtilde_t * dRtilde_dbeta_t_k * dDtilde_ddelta_t_l;
  return ddHtilde_dbetaddeltaprime_t_kl;
}

double fc_ddltilde_dthetadthetaprime_t_kl(int m, arma::vec y_t, arma::mat Htilde_t, arma::mat dHtilde_dtheta_t_l, arma::mat dHtilde_dtheta_t_k, arma::mat ddHtilde_dthetadthetaprime_t_kl) {
  // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\theta_k\partial\theta_\ell
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat y_t_mat = fc_asmat(y_t, m, 1);
  arma::mat inverse_Htilde_t = Htilde_t.i();
  double trace1 = trace((I_m - 2.0 * inverse_Htilde_t * y_t * y_t_mat.t()) * inverse_Htilde_t * dHtilde_dtheta_t_k * inverse_Htilde_t * dHtilde_dtheta_t_l);
  double trace2 = trace((I_m - inverse_Htilde_t * y_t * y_t_mat.t()) * inverse_Htilde_t * ddHtilde_dthetadthetaprime_t_kl);
  double ddltilde_dthetadthetaprime_t_kl = - trace1 + trace2;
  return ddltilde_dthetadthetaprime_t_kl;
}

// [[Rcpp::export]]
arma::mat fc_derivative_of_transformation_function(int m, int r, int s, arma::vec vartheta_vec) {
  // the derivative matrix of transformation function, that is \partial\bm\theta / \partial\bm\vartheta^\prime
  int Dim1_theta = m+r+2*s; int Dim2_theta = Dim1_theta + r*m*m; int Dim3_theta = Dim2_theta + s*m*m; int Dim4_theta = Dim3_theta + s*m*m;
  int Dim1_vartheta = m+r+2*s; int Dim2_vartheta = Dim1_vartheta + 2*r*m; int Dim3_vartheta = Dim2_vartheta + 4*s*m; int Dim4_vartheta = Dim3_vartheta + 4*s*m;
  int dimtheta = Dim4_theta + (2 + m*(m-1)/2);
  int dimvartheta = Dim4_vartheta + (2 + m*(m-1)/2);
  arma::mat Delta_theta_vartheta(dimtheta, dimvartheta); Delta_theta_vartheta.fill(0.0);
  arma::mat I_front(m+r+2*s, m+r+2*s); I_front.eye(m+r+2*s, m+r+2*s);
  arma::mat I_last(2+m*(m-1)/2, 2+m*(m-1)/2); I_last.eye(2+m*(m-1)/2, 2+m*(m-1)/2);
  Delta_theta_vartheta.submat(0, 0, Dim1_theta-1, Dim1_vartheta-1) = I_front;
  Delta_theta_vartheta.submat(Dim4_theta, Dim4_vartheta, dimtheta-1, dimvartheta-1) = I_last;
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::vec g_in_vartheta = vartheta_vec.subvec(Dim1_vartheta, Dim4_vartheta-1); // (\bm{g}_0', \bm{g}_1', \bm{g}_2')' in \bm\vartheta
  if(s == 0) {
    arma::vec g_0_vec_in_vartheta = g_in_vartheta.subvec(0, 2*r*m-1); arma::mat g_0_mat_in_vartheta = fc_asmat(g_0_vec_in_vartheta, 2*m, r);
    for(int k = 0; k < r; k++) {
      int dimk0_theta = Dim1_theta + k*m*m; 
      int dimk0_vartheta = Dim1_vartheta + k*2*m; 
      arma::vec g_0k = g_0_mat_in_vartheta.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2); 
      arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
      arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
      Delta_theta_vartheta.submat(dimk0_theta, dimk0_vartheta, dimk0_theta+m*m-1, dimk0_vartheta+2*m-1) = join_rows(kron(g_0k2, I_m), kron(I_m, g_0k1));
    }
  } else if(r == 0) {
    arma::vec g_1_vec_in_vartheta = g_in_vartheta.subvec(2*r*m, 2*r*m+4*s*m-1); arma::mat g_1_mat_in_vartheta = fc_asmat(g_1_vec_in_vartheta, 4*m, s);
    arma::vec g_2_vec_in_vartheta = g_in_vartheta.subvec(2*r*m+4*s*m, 2*r*m+8*s*m-1); arma::mat g_2_mat_in_vartheta = fc_asmat(g_2_vec_in_vartheta, 4*m, s);
    for(int k = 0; k < s; k++) {
      int dimk1_theta = Dim2_theta + k*m*m; 
      int dimk1_vartheta = Dim2_vartheta + k*4*m; 
      int dimk2_theta = Dim3_theta + k*m*m; 
      int dimk2_vartheta = Dim3_vartheta + k*4*m; 
      arma::vec g_1k = g_1_mat_in_vartheta.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4); 
      arma::vec g_2k = g_2_mat_in_vartheta.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4); 
      arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
      arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
      arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
      arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
      arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
      arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
      arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
      arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
      Delta_theta_vartheta.submat(dimk1_theta, dimk1_vartheta, dimk1_theta+m*m-1, dimk1_vartheta+4*m-1) = join_rows(kron(g_1k2, I_m), kron(I_m, g_1k1), kron(g_1k4, I_m), kron(I_m, g_1k3));
      Delta_theta_vartheta.submat(dimk2_theta, dimk2_vartheta, dimk2_theta+m*m-1, dimk2_vartheta+4*m-1) = join_rows(kron(g_2k2, I_m), kron(I_m, g_2k1), kron(g_2k4, I_m), kron(I_m, g_2k3));
    }
  } else {
    arma::vec g_0_vec_in_vartheta = g_in_vartheta.subvec(0, 2*r*m-1); arma::mat g_0_mat_in_vartheta = fc_asmat(g_0_vec_in_vartheta, 2*m, r);
    arma::vec g_1_vec_in_vartheta = g_in_vartheta.subvec(2*r*m, 2*r*m+4*s*m-1); arma::mat g_1_mat_in_vartheta = fc_asmat(g_1_vec_in_vartheta, 4*m, s);
    arma::vec g_2_vec_in_vartheta = g_in_vartheta.subvec(2*r*m+4*s*m, 2*r*m+8*s*m-1); arma::mat g_2_mat_in_vartheta = fc_asmat(g_2_vec_in_vartheta, 4*m, s);
    for(int k = 0; k < r; k++) {
      int dimk0_theta = Dim1_theta + k*m*m; 
      int dimk0_vartheta = Dim1_vartheta + k*2*m; 
      arma::vec g_0k = g_0_mat_in_vartheta.col(k); arma::mat g_0k_mat = fc_asmat(g_0k, m, 2); 
      arma::mat g_0k1(m, 1); g_0k1.col(0) = g_0k_mat.col(0);
      arma::mat g_0k2(m, 1); g_0k2.col(0) = g_0k_mat.col(1);
      Delta_theta_vartheta.submat(dimk0_theta, dimk0_vartheta, dimk0_theta+m*m-1, dimk0_vartheta+2*m-1) = join_rows(kron(g_0k2, I_m), kron(I_m, g_0k1));
    }
    for(int k = 0; k < s; k++) {
      int dimk1_theta = Dim2_theta + k*m*m; 
      int dimk1_vartheta = Dim2_vartheta + k*4*m; 
      int dimk2_theta = Dim3_theta + k*m*m; 
      int dimk2_vartheta = Dim3_vartheta + k*4*m; 
      arma::vec g_1k = g_1_mat_in_vartheta.col(k); arma::mat g_1k_mat = fc_asmat(g_1k, m, 4); 
      arma::vec g_2k = g_2_mat_in_vartheta.col(k); arma::mat g_2k_mat = fc_asmat(g_2k, m, 4); 
      arma::mat g_1k1(m, 1); g_1k1.col(0) = g_1k_mat.col(0);
      arma::mat g_1k2(m, 1); g_1k2.col(0) = g_1k_mat.col(1);
      arma::mat g_1k3(m, 1); g_1k3.col(0) = g_1k_mat.col(2);
      arma::mat g_1k4(m, 1); g_1k4.col(0) = g_1k_mat.col(3);
      arma::mat g_2k1(m, 1); g_2k1.col(0) = g_2k_mat.col(0);
      arma::mat g_2k2(m, 1); g_2k2.col(0) = g_2k_mat.col(1);
      arma::mat g_2k3(m, 1); g_2k3.col(0) = g_2k_mat.col(2);
      arma::mat g_2k4(m, 1); g_2k4.col(0) = g_2k_mat.col(3);
      Delta_theta_vartheta.submat(dimk1_theta, dimk1_vartheta, dimk1_theta+m*m-1, dimk1_vartheta+4*m-1) = join_rows(kron(g_1k2, I_m), kron(I_m, g_1k1), kron(g_1k4, I_m), kron(I_m, g_1k3));
      Delta_theta_vartheta.submat(dimk2_theta, dimk2_vartheta, dimk2_theta+m*m-1, dimk2_vartheta+4*m-1) = join_rows(kron(g_2k2, I_m), kron(I_m, g_2k1), kron(g_2k4, I_m), kron(I_m, g_2k3));
    }
  }
  return Delta_theta_vartheta; 
}

// [[Rcpp::export]]
arma::vec fc_lowrank_ASD_thetahat(int n, int m, int r, int s, int Bk, arma::vec thetahat, arma::mat y_m_n, arma::mat Deltahat_theta_vartheta) {
  // ASD of \widehat\bm\theta
  int dimkappa = r + 2*s + r*m*m + 2*s*m*m;
  int dimbeta = 2 + m*(m-1)/2;
  int dimdelta = m + dimkappa;
  int dimtheta = dimdelta + dimbeta;
  arma::vec omegauline = thetahat.subvec(0, m-1); // \underline{\bm\omega}
  arma::vec kappa_vec = thetahat.subvec(m, m+dimkappa-1); // \bm\kappa
  arma::vec beta_vec = thetahat.subvec(m+dimkappa, m+dimkappa+dimbeta-1); // \bm\beta
  double beta_1 = beta_vec(0); double beta_2 = beta_vec(1); // \beta_1 and \beta_2
  arma::vec ruline = beta_vec.tail(dimbeta-2); // \underline\bm{r}
  arma::mat Ruline = fc_Ruline(m, ruline); // \underline{R}
  // summations using FFT algorithm
  int max_r_1 = r; if(max_r_1 == 0) {max_r_1 = 1;}
  int max_s_1 = s; if(max_s_1 == 0) {max_s_1 = 1;}
  arma::cube sum_lambdaklnyuline_0_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_0_m_nminus1_r.fill(0.0);
  arma::cube sum_lambdaklnyuline_1_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_1_m_nminus1_r.fill(0.0);
  arma::cube sum_lambdaklnyuline_2_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_2_m_nminus1_r.fill(0.0);
  arma::cube sum_gammaklnyuline_01_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_01_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_02_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_02_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_11_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_11_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_12_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_12_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_21_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_21_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_22_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_22_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_31_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_31_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_32_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_32_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_41_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_41_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_42_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_42_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_51_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_51_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_52_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_52_m_nminus1_s.fill(0.0);
  if(s == 0) {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_2_m_nminus1_r = fc_sum_lambdaklnyuline_2_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
  } else if(r == 0) {
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_31_m_nminus1_s = fc_sum_gammaklnyuline_31_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_32_m_nminus1_s = fc_sum_gammaklnyuline_32_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_41_m_nminus1_s = fc_sum_gammaklnyuline_41_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_42_m_nminus1_s = fc_sum_gammaklnyuline_42_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_51_m_nminus1_s = fc_sum_gammaklnyuline_51_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_52_m_nminus1_s = fc_sum_gammaklnyuline_52_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  } else {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_2_m_nminus1_r = fc_sum_lambdaklnyuline_2_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_31_m_nminus1_s = fc_sum_gammaklnyuline_31_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_32_m_nminus1_s = fc_sum_gammaklnyuline_32_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_41_m_nminus1_s = fc_sum_gammaklnyuline_41_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_42_m_nminus1_s = fc_sum_gammaklnyuline_42_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_51_m_nminus1_s = fc_sum_gammaklnyuline_51_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_52_m_nminus1_s = fc_sum_gammaklnyuline_52_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  }
  // n \widehat{\Sigma}
  arma::mat nSigmahat(dimtheta, dimtheta); nSigmahat.fill(0.0); // n \widehat{\Sigma} = \sum_{t=1}^{n} \widehat{\Sigma}_t = \sum_{t=1}^{n} \partial\widetilde{\ell}_t(\widehat\bm\theta)/\partial\bm\theta \partial\widetilde{\ell}_t(\widehat\bm\theta)/\partial\bm\theta'
  // n \widehat{\Sigma}_{*}
  arma::mat nSigmaStarhat(dimtheta, dimtheta); nSigmaStarhat.fill(0.0); // n \widehat{\Sigma}_{*} = \sum_{t=1}^{n} \widehat{\Sigma}_{*, t} = \sum_{t=1}^{n} \partial^2\widetilde{\ell}_t(\widehat\bm\theta) / \partial\bm\theta\partial\bm\theta'
  // t=1
  arma::vec y_1 = y_m_n.col(0); // \bm{y}_{1}
  arma::mat Dtilde_1 = fc_Dtilde_1(omegauline); // \widetilde{D}_{1}(\bm\delta)
  arma::mat inverse_Dtilde_1 = fc_inverse_Dtilde_1(omegauline); // \widetilde{D}_{1}^{-1}(\bm\delta)
  arma::mat Rtilde_1 = fc_Rtilde_1(m, beta_1, beta_2, Ruline); // \widetilde{R}_{1}(\bm\theta)
  arma::mat Htilde_1 = fc_Htilde_t(Dtilde_1, Rtilde_1); // \widetilde{H}_{1}(\bm\theta)
  arma::vec varepsilontilde_1 = inverse_Dtilde_1 * y_1; // \widetilde{\varepsilon}_{1}(\bm\delta) = \widetilde{D}_{1}^{-1}(\bm\delta) \bm{y}_{1}
  // first derivatives
  arma::mat dlnhulinetilde_ddelta_1_mat = fc_general_dlnhulinetilde_ddelta_1_mat(m, dimdelta); // \partial\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\bm\delta'_\ell, which is equal to \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'_\ell for t <= 0
  arma::cube dDtilde_ddelta_1_cube = fc_dDtilde_ddelta_t_cube(m, dimdelta, Dtilde_1, dlnhulinetilde_ddelta_1_mat); // the \ell-th slice is \partial\widetilde{D}_1(\bm\delta) / \partial\delta_\ell
  arma::cube dRtilde_ddelta_1_cube = fc_dRtilde_ddelta_1_cube(m, dimdelta); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\delta_\ell
  arma::cube dRuline_dbeta_cube = fc_dRuline_dbeta_cube(m, dimbeta); // the \ell-th slice is \partial\underline{R} / \partial\beta_\ell
  arma::cube dRtilde_dbeta_1_cube = fc_dRtilde_dbeta_1_cube(m, dimbeta, beta_1, beta_2, Ruline, dRuline_dbeta_cube); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\beta_\ell
  arma::cube dHtilde_ddelta_1_cube = fc_dHtilde_ddelta_t_cube(m, dimdelta, Dtilde_1, Rtilde_1, dDtilde_ddelta_1_cube, dRtilde_ddelta_1_cube); // the \ell-th slice is \partial\widetilde{H}_1(\bm\theta) / \partial\delta_\ell
  arma::cube dHtilde_dbeta_1_cube = fc_dHtilde_dbeta_t_cube(m, dimbeta, Dtilde_1, dRtilde_dbeta_1_cube); // the \ell-th slice is \partial\widetilde{H}_1(\bm\theta) / \partial\beta_\ell
  // \widehat{\Sigma}_1
  arma::vec dltilde_dtheta_1 = fc2_dltilde_dtheta_t(m, dimdelta, dimbeta, dimtheta, y_1, Htilde_1, dHtilde_ddelta_1_cube, dHtilde_dbeta_1_cube); // \partial\widetilde{\ell}_1(\bm\theta) / \partial\bm\theta
  arma::mat dltilde_dtheta_1_mat = fc_asmat(dltilde_dtheta_1, dimtheta, 1);
  arma::mat Sigmahat_1 = dltilde_dtheta_1 * dltilde_dtheta_1_mat.t(); // \widehat{\Sigma}_1
  nSigmahat = nSigmahat + Sigmahat_1;
  // second derivatives
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_1_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_1_cube(m, dimdelta); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\delta_k\partial\delta_\ell, which is equal to \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell for t <= 0
  arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_1_field = fc_ddRtilde_ddeltaddeltaprime_1_field(m, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddRtilde_dbetadbetaprime_1_field = fc_ddRtilde_dbetadbetaprime_1_field(m, dimbeta, beta_1, beta_2, Ruline, dRuline_dbeta_cube); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::field<arma::mat> ddRtilde_dbetaddeltaprime_1_field = fc_ddRtilde_dbetaddeltaprime_1_field(m, dimdelta, dimbeta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\beta_k\partial\delta_\ell
  // \partial^2\ln\widetilde\underline{\bm{h}}_{1-\Bbbk}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{0}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix, of which the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{-\Bbbk+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  // \partial^2\widetilde{\ell}_{1}(\bm\theta) / \partial\bm\theta \partial\bm\theta'
  arma::mat ddltilde_dthetadthetaprime_1(dimtheta, dimtheta); ddltilde_dthetadthetaprime_1.fill(0.0);
  for(int k = 0; k < dimdelta; k++) {
    for(int l = 0; l < (k+1); l++) {
      arma::vec dlnhulinetilde_ddelta_1_l = dlnhulinetilde_ddelta_1_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
      arma::vec dlnhulinetilde_ddelta_1_k = dlnhulinetilde_ddelta_1_mat.col(k); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k
      arma::mat ddlnhulinetilde_ddeltaddeltaprime_1_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_1_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\bm\delta'
      arma::vec ddlnhulinetilde_ddeltaddeltaprime_1_kl = ddlnhulinetilde_ddeltaddeltaprime_1_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat(m, Bk); ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat.fill(0.0); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{-\Bbbk+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
      ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat.each_col() = ddlnhulinetilde_ddeltaddeltaprime_1_kl;
      ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(k, l) = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(l, k) = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat;
      arma::mat dDtilde_ddelta_1_l = dDtilde_ddelta_1_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
      arma::mat dDtilde_ddelta_1_k = dDtilde_ddelta_1_cube.slice(k); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_k
      arma::mat dRtilde_ddelta_1_l = dRtilde_ddelta_1_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dRtilde_ddelta_1_k = dRtilde_ddelta_1_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_k
      arma::mat dHtilde_ddelta_1_l = dHtilde_ddelta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dHtilde_ddelta_1_k = dHtilde_ddelta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_k
      arma::mat ddDtilde_ddeltaddeltaprime_1_kl = fc_ddDtilde_ddeltaddeltaprime_t_kl(dlnhulinetilde_ddelta_1_l, dlnhulinetilde_ddelta_1_k, ddlnhulinetilde_ddeltaddeltaprime_1_kl, Dtilde_1); // \partial^2\widetilde{D}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddRtilde_ddeltaddeltaprime_1_kl = ddRtilde_ddeltaddeltaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddHtilde_ddeltaddeltaprime_1_kl = fc_ddHtilde_ddeltaddeltaprime_t_kl(Dtilde_1, Rtilde_1, dDtilde_ddelta_1_l, dDtilde_ddelta_1_k, dRtilde_ddelta_1_l, dRtilde_ddelta_1_k, ddDtilde_ddeltaddeltaprime_1_kl, ddRtilde_ddeltaddeltaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
      ddltilde_dthetadthetaprime_1(k,l) = ddltilde_dthetadthetaprime_1(l,k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_ddelta_1_l, dHtilde_ddelta_1_k, ddHtilde_ddeltaddeltaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\delta_k\partial\delta_\ell
    }
  }
  for(int k = 0; k < dimbeta; k++) {
    for(int l = 0; l < (k+1); l++) {
      arma::mat dHtilde_dbeta_1_l = dHtilde_dbeta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
      arma::mat dHtilde_dbeta_1_k = dHtilde_dbeta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
      arma::mat ddRtilde_dbetadbetaprime_1_kl = ddRtilde_dbetadbetaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
      arma::mat ddHtilde_dbetadbetaprime_1_kl = fc_ddHtilde_dbetadbetaprime_t_kl(Dtilde_1, ddRtilde_dbetadbetaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
      ddltilde_dthetadthetaprime_1(dimdelta+k,dimdelta+l) = ddltilde_dthetadthetaprime_1(dimdelta+l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_dbeta_1_l, dHtilde_dbeta_1_k, ddHtilde_dbetadbetaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\beta_\ell
    }
  }
  for(int k = 0; k < dimbeta; k++) {
    for(int l = 0; l < dimdelta; l++) {
      arma::mat dDtilde_ddelta_1_l = dDtilde_ddelta_1_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
      arma::mat dRtilde_dbeta_1_k = dRtilde_dbeta_1_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_k
      arma::mat dHtilde_ddelta_1_l = dHtilde_ddelta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dHtilde_dbeta_1_k = dHtilde_dbeta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
      arma::mat ddRtilde_dbetaddeltaprime_1_kl = ddRtilde_dbetaddeltaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
      arma::mat ddHtilde_dbetaddeltaprime_1_kl = fc_ddHtilde_dbetaddeltaprime_t_kl(Dtilde_1, dDtilde_ddelta_1_l, dRtilde_dbeta_1_k, ddRtilde_dbetaddeltaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
      ddltilde_dthetadthetaprime_1(dimdelta+k,l) = ddltilde_dthetadthetaprime_1(l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_ddelta_1_l, dHtilde_dbeta_1_k, ddHtilde_dbetaddeltaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\delta_\ell
    }
  }
  // \widehat{\Sigma}_{*, 1}
  arma::mat SigmaStarhat_1 = ddltilde_dthetadthetaprime_1; // \widehat{\Sigma}_{*, 1}
  nSigmaStarhat = nSigmaStarhat + SigmaStarhat_1;
  // t=2,...,n
  // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
  arma::mat varepsilontilde_tminusBkminus1TOtminus2_mat(m, Bk); varepsilontilde_tminusBkminus1TOtminus2_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta))
  varepsilontilde_tminusBkminus1TOtminus2_mat.each_col() = exp(-0.5*omegauline); // for s <= 0, \widetilde{\varepsilon}_s(\bm\delta) = (\exp{-1/2*\underliner{\omega}_1}, ..., \exp{-1/2*\underliner{\omega}_m})'
  arma::vec varepsilontilde_tminus1 = varepsilontilde_1; // \widetilde{\varepsilon}_{t-1}(\bm\delta)
  // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
  arma::cube dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.fill(0.0); // the i-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-2+i}(\bm\delta) / \partial\bm\delta'_\ell
  dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.each_slice() = dlnhulinetilde_ddelta_1_mat;
  arma::mat dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_1_mat;
  // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_k\partial\delta_\ell, and \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field; // the (k, \ell, i)-th element is the vector \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-2+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube = ddlnhulinetilde_ddeltaddeltaprime_1_cube;
  // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell, \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell, and \partial^2\widetilde{R}_{t-1}(\bm\theta)/\partial\theta_\k\partial\theta_\ell
  arma::mat Rtilde_tminus1 = Rtilde_1; // \widetilde{R}_{t-1}(\bm\theta)
  arma::cube dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
  arma::cube dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_\ell
  arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_tminus1_field = ddRtilde_ddeltaddeltaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddRtilde_dbetadbetaprime_tminus1_field = ddRtilde_dbetadbetaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::field<arma::mat> ddRtilde_dbetaddeltaprime_tminus1_field = ddRtilde_dbetaddeltaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\delta_\ell
  for(int t = 2; t < (n+1); t++) {
    arma::vec y_t = y_m_n.col(t-1); // \bm{y}_{t}
    arma::vec lnhulinetilde_t = fc_general_lnhulinetilde_t_general(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s); // \ln\widetilde\underline{\bm{h}}_t(\bm\delta)
    arma::mat Dtilde_t = fc_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}(\bm\delta)
    arma::mat inverse_Dtilde_t = fc_inverse_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}^{-1}(\bm\delta)
    arma::vec varepsilontilde_t = inverse_Dtilde_t * y_t; // \widetilde{\varepsilon}_{t}(\bm\delta) = \widetilde{D}_{t}^{-1}(\bm\delta) \bm{y}_{t}
    arma::mat varepsilontilde_tminusBkTOtminus1_mat(m, Bk); varepsilontilde_tminusBkTOtminus1_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk}(\bm\delta), ..., \widetilde{\varepsilon}_{t-1}(\bm\delta))
    varepsilontilde_tminusBkTOtminus1_mat.cols(0, Bk-2) = varepsilontilde_tminusBkminus1TOtminus2_mat.cols(1, Bk-1);
    varepsilontilde_tminusBkTOtminus1_mat.col(Bk-1) = varepsilontilde_tminus1;
    arma::mat Psitilde_tminus1 = fc_Psi(m, varepsilontilde_tminusBkTOtminus1_mat); // \widetilde{\Psi}_{t-1}(\bm\delta)
    arma::mat Rtilde_t = fc_Rtilde_t(m, Bk, beta_1, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1); // \widetilde{R}_{t}(\bm\theta)
    arma::mat Htilde_t = fc_Htilde_t(Dtilde_t, Rtilde_t); // \widetilde{H}_{t}(\bm\theta)
    // first derivatives
    arma::mat dlnhulinetilde_ddelta_t_mat = fc_general_dlnhulinetilde_ddelta_t_mat_general(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_lambdaklnyuline_1_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'
    arma::cube dlnhulinetilde_ddelta_tminusBkTOtminus1_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.fill(0.0); // the k-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-k)}(\bm\delta) / \partial\bm\delta'
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slices(0, Bk-2) = dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.slices(1, Bk-1);
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slice(Bk-1) = dlnhulinetilde_ddelta_tminus1_mat;
    arma::cube dDtilde_ddelta_t_cube = fc_dDtilde_ddelta_t_cube(m, dimdelta, Dtilde_t, dlnhulinetilde_ddelta_t_mat); // the \ell-th slice is \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
    arma::cube dPsitilde_ddelta_tminus1_cube = fc_dPsitilde_ddelta_tminus1_cube(m, dimdelta, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dlnhulinetilde_ddelta_tminusBkTOtminus1_cube); // the \ell-th slice is \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
    arma::cube dRtilde_ddelta_t_cube = fc2_dRtilde_ddelta_t_cube(m, dimdelta, beta_1, beta_2, dPsitilde_ddelta_tminus1_cube, dRtilde_ddelta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
    arma::cube dRtilde_dbeta_t_cube = fc_dRtilde_dbeta_t_cube(m, dimbeta, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1, dRtilde_dbeta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\beta_\ell
    arma::cube dHtilde_ddelta_t_cube = fc_dHtilde_ddelta_t_cube(m, dimdelta, Dtilde_t, Rtilde_t, dDtilde_ddelta_t_cube, dRtilde_ddelta_t_cube); // the \ell-th slice is \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
    arma::cube dHtilde_dbeta_t_cube = fc_dHtilde_dbeta_t_cube(m, dimbeta, Dtilde_t, dRtilde_dbeta_t_cube); // the \ell-th slice is \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
    // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
    arma::vec dltilde_dtheta_t = fc2_dltilde_dtheta_t(m, dimdelta, dimbeta, dimtheta, y_t, Htilde_t, dHtilde_ddelta_t_cube, dHtilde_dbeta_t_cube); // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
    // \widehat{\Sigma}_t
    arma::mat dltilde_dtheta_t_mat = fc_asmat(dltilde_dtheta_t, dimtheta, 1);
    arma::mat Sigmahat_t = dltilde_dtheta_t * dltilde_dtheta_t_mat.t(); // \widehat{\Sigma}_t
    nSigmahat = nSigmahat + Sigmahat_t;
    // second derivatives
    arma::cube ddlnhulinetilde_ddeltaddeltaprime_t_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube_general(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_1_m_nminus1_r, sum_lambdaklnyuline_2_m_nminus1_r, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s, sum_gammaklnyuline_31_m_nminus1_s, sum_gammaklnyuline_32_m_nminus1_s, sum_gammaklnyuline_41_m_nminus1_s, sum_gammaklnyuline_42_m_nminus1_s, sum_gammaklnyuline_51_m_nminus1_s, sum_gammaklnyuline_52_m_nminus1_s); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_t_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddRtilde_dbetadbetaprime_t_field(dimbeta, dimbeta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
    arma::field<arma::mat> ddRtilde_dbetaddeltaprime_t_field(dimbeta, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
    // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix, of which the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta \partial\bm\theta'
    arma::mat ddltilde_dthetadthetaprime_t(dimtheta, dimtheta); ddltilde_dthetadthetaprime_t.fill(0.0);
    for(int k = 0; k < dimdelta; k++) {
      for(int l = 0; l < (k+1); l++) {
        arma::vec dlnhulinetilde_ddelta_t_l = dlnhulinetilde_ddelta_t_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
        arma::vec dlnhulinetilde_ddelta_t_k = dlnhulinetilde_ddelta_t_mat.col(k); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k
        arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.col(l); // the i-th column is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_\ell
        arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.col(k); // the i-th column is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_k
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_t_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\bm\delta'
        arma::vec ddlnhulinetilde_ddeltaddeltaprime_t_kl = ddlnhulinetilde_ddeltaddeltaprime_t_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field(k, l); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\bm\delta'
        arma::vec ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl = ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat(m, Bk); ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.fill(0.0); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.cols(0, Bk-2) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_kl_mat.cols(1, Bk-1);
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.col(Bk-1) = ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl;
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(k, l) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(l, k) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat;
        arma::mat dDtilde_ddelta_t_l = dDtilde_ddelta_t_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
        arma::mat dDtilde_ddelta_t_k = dDtilde_ddelta_t_cube.slice(k); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_k
        arma::mat dPsitilde_ddelta_tminus1_k = dPsitilde_ddelta_tminus1_cube.slice(k); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_k
        arma::mat dRtilde_ddelta_t_l = dRtilde_ddelta_t_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dRtilde_ddelta_t_k = dRtilde_ddelta_t_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_k
        arma::mat dHtilde_ddelta_t_l = dHtilde_ddelta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_ddelta_t_k = dHtilde_ddelta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_k
        arma::mat ddDtilde_ddeltaddeltaprime_t_kl = fc_ddDtilde_ddeltaddeltaprime_t_kl(dlnhulinetilde_ddelta_t_l, dlnhulinetilde_ddelta_t_k, ddlnhulinetilde_ddeltaddeltaprime_t_kl, Dtilde_t); // \partial^2\widetilde{D}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddRtilde_ddeltaddeltaprime_tminus1_kl = ddRtilde_ddeltaddeltaprime_tminus1_field(k,l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddRtilde_ddeltaddeltaprime_t_kl = fc_ddRtilde_ddeltaddeltaprime_t_kl(m, beta_1, beta_2, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dPsitilde_ddelta_tminus1_k, dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat, dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat, ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat, ddRtilde_ddeltaddeltaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
        ddRtilde_ddeltaddeltaprime_t_field(k, l) = ddRtilde_ddeltaddeltaprime_t_field(l, k) = ddRtilde_ddeltaddeltaprime_t_kl;
        arma::mat ddHtilde_ddeltaddeltaprime_t_kl = fc_ddHtilde_ddeltaddeltaprime_t_kl(Dtilde_t, Rtilde_t, dDtilde_ddelta_t_l, dDtilde_ddelta_t_k, dRtilde_ddelta_t_l, dRtilde_ddelta_t_k, ddDtilde_ddeltaddeltaprime_t_kl, ddRtilde_ddeltaddeltaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
        ddltilde_dthetadthetaprime_t(k,l) = ddltilde_dthetadthetaprime_t(l,k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_ddelta_t_l, dHtilde_ddelta_t_k, ddHtilde_ddeltaddeltaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\delta_k\partial\delta_\ell
      }
    }
    for(int k = 0; k < dimbeta; k++) {
      for(int l = 0; l < (k+1); l++) {
        arma::mat dRtilde_dbeta_tminus1_l = dRtilde_dbeta_tminus1_cube.slice(l); // \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_\ell
        arma::mat dHtilde_dbeta_t_l = dHtilde_dbeta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
        arma::mat dHtilde_dbeta_t_k = dHtilde_dbeta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
        arma::mat ddRtilde_dbetadbetaprime_tminus1_kl = ddRtilde_dbetadbetaprime_tminus1_field(k, l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\beta_\ell
        arma::mat ddRtilde_dbetadbetaprime_t_kl = fc_ddRtilde_dbetadbetaprime_t_kl(l, k, m, beta_2, dRtilde_dbeta_tminus1_l, ddRtilde_dbetadbetaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
        ddRtilde_dbetadbetaprime_t_field(k, l) = ddRtilde_dbetadbetaprime_t_field(l, k) = ddRtilde_dbetadbetaprime_t_kl;
        arma::mat ddHtilde_dbetadbetaprime_t_kl = fc_ddHtilde_dbetadbetaprime_t_kl(Dtilde_t, ddRtilde_dbetadbetaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
        ddltilde_dthetadthetaprime_t(dimdelta+k,dimdelta+l) = ddltilde_dthetadthetaprime_t(dimdelta+l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_dbeta_t_l, dHtilde_dbeta_t_k, ddHtilde_dbetadbetaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\beta_\ell
      }
    }
    for(int k = 0; k < dimbeta; k++) {
      for(int l = 0; l < dimdelta; l++) {
        arma::mat dDtilde_ddelta_t_l = dDtilde_ddelta_t_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
        arma::mat dPsitilde_ddelta_tminus1_l = dPsitilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
        arma::mat dRtilde_dbeta_t_k = dRtilde_dbeta_t_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_k
        arma::mat dRtilde_ddelta_tminus1_l = dRtilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_ddelta_t_l = dHtilde_ddelta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_dbeta_t_k = dHtilde_dbeta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
        arma::mat ddRtilde_dbetaddeltaprime_tminus1_kl = ddRtilde_dbetaddeltaprime_tminus1_field(k, l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\delta_\ell
        arma::mat ddRtilde_dbetaddeltaprime_t_kl = fc_ddRtilde_dbetaddeltaprime_t_kl(l, k, m, beta_2, dPsitilde_ddelta_tminus1_l, dRtilde_ddelta_tminus1_l, ddRtilde_dbetaddeltaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
        ddRtilde_dbetaddeltaprime_t_field(k, l) = ddRtilde_dbetaddeltaprime_t_kl;
        arma::mat ddHtilde_dbetaddeltaprime_t_kl = fc_ddHtilde_dbetaddeltaprime_t_kl(Dtilde_t, dDtilde_ddelta_t_l, dRtilde_dbeta_t_k, ddRtilde_dbetaddeltaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
        ddltilde_dthetadthetaprime_t(dimdelta+k,l) = ddltilde_dthetadthetaprime_t(l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_ddelta_t_l, dHtilde_dbeta_t_k, ddHtilde_dbetaddeltaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\delta_\ell
      }
    }
    // \widehat{\Sigma}_{*, t}
    arma::mat SigmaStarhat_t = ddltilde_dthetadthetaprime_t; // \widehat{\Sigma}_{*, t}
    nSigmaStarhat = nSigmaStarhat + SigmaStarhat_t;
    // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
    varepsilontilde_tminusBkminus1TOtminus2_mat = varepsilontilde_tminusBkTOtminus1_mat;
    varepsilontilde_tminus1 = varepsilontilde_t;
    // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
    dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube;
    dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_t_mat;
    // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_k\partial\delta_\ell and \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field;
    ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube = ddlnhulinetilde_ddeltaddeltaprime_t_cube;
    // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell, \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell, and \partial^2\widetilde{R}_{t-1}(\bm\theta)/\partial\theta_\k\partial\theta_\ell
    Rtilde_tminus1 = Rtilde_t;
    dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_t_cube;
    dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_t_cube;
    ddRtilde_ddeltaddeltaprime_tminus1_field = ddRtilde_ddeltaddeltaprime_t_field;
    ddRtilde_dbetadbetaprime_tminus1_field = ddRtilde_dbetadbetaprime_t_field;
    ddRtilde_dbetaddeltaprime_tminus1_field = ddRtilde_dbetaddeltaprime_t_field;
  }
  arma::mat Sigmahat = nSigmahat / (n * 1.0); // \widehat{\Sigma} = 1/n \sum_{t=1}^{n} \widehat{\Sigma}_t
  arma::mat SigmaStarhat = nSigmaStarhat / (n * 1.0); // \widehat{\Sigma}_{*} = 1/n \sum_{t=1}^{n} \widehat{\Sigma}_{*, t}
  arma::mat Sigmahat_transf = Deltahat_theta_vartheta.t() * Sigmahat * Deltahat_theta_vartheta; // \Delta_\theta' \widehat{\Sigma} \Delta_\theta
  arma::mat SigmaStarhat_transf = Deltahat_theta_vartheta.t() * SigmaStarhat * Deltahat_theta_vartheta; // \Delta_\theta' \widehat{\Sigma}_{*} \Delta_\theta
  arma::mat CovMat = Deltahat_theta_vartheta * pinv(SigmaStarhat_transf) * Sigmahat_transf * pinv(SigmaStarhat_transf) * Deltahat_theta_vartheta.t(); // \Delta_\theta (\Delta_\theta' \widehat{\Sigma}_{*} \Delta_\theta)^{-1} (\Delta_\theta' \widehat{\Sigma} \Delta_\theta) (\Delta_\theta' \widehat{\Sigma}_{*} \Delta_\theta)^{-1} \Delta_\theta'
  arma::vec var_sqrtnthetahat = CovMat.diag(0); // variances of \sqrt{n} \widehat\bm\theta
  arma::vec var_thetahat = var_sqrtnthetahat / (n * 1.0); // variances of \widehat\bm\theta
  arma::vec ASD_thetahat = sqrt(var_thetahat); // ASD of \widehat\bm\theta
  return ASD_thetahat;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////// Estimation \widehat\bm\theta //////////////////////////////////////////////////
////// In general,
////// \bm\theta = (\bm\delta', \bm\beta')'
//////           = (\underline\bm\omega'(m*1), \bm\kappa', \beta_1, \beta_2, \underline\bmr'(m(m-1)/2*1))'
//////           = (\underline\bm\omega'(m*1), \bm\lambda'(r*1), \bm\gamma'(s*1), \bm\varphi'(s*1), \bm{g}_0'(rm^2*1), \bm{g}_1'(sm^2*1), \bm{g}_2'(sm^2*1), \beta_1, \beta_2, \underline\bm{r}'(m(m-1)/2*1))'
////// \bm{g}_0 = (\bm{g}_{01}', ..., \bm{g}_{0r}')', \bm{g}_1 = (\bm{g}_{11}', ..., \bm{g}_{1s}')', \bm{g}_2 = (\bm{g}_{21}', ..., \bm{g}_{2s}')' with 
////// \bm{g}_{0k} = vec(G_{0k}),
////// \bm{g}_{1k} = vec(G_{1k}),
////// \bm{g}_{2k} = vec(G_{2k}).
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////// constraints of parameters //////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat fc_general_linear_constr_ui(int m, int r, int s, int dimtheta) {
  // linear inequality constraints, ui %*% theta - ci >=0
  int dimruline = m*(m-1)/2;
  int dimchi = dimtheta - dimruline;
  arma::vec dim_inequality_lambda_choices(2); dim_inequality_lambda_choices(0) = r-1; dim_inequality_lambda_choices(1) = 0;
  int dim_inequality_lambda = max(dim_inequality_lambda_choices);
  int dim_constr_lambda = 2*r + 4*s + dim_inequality_lambda;
  int dim_constr_beta12 = 5;
  int dim_constr_ruline = 2*dimruline;
  arma::mat ui(dim_constr_lambda+dim_constr_beta12+dim_constr_ruline, dimtheta); ui.fill(0.0);
  if(s == 0) {
    for(int k = 0; k < r; k++) {
      ui(2*k, m+k) = -1.0;
      ui(2*k+1, m+k) = 1.0;
    }
    for(int k = 1; k < r; k++) {
      ui(2*r+4*s+k-1, m+k-1) = 1.0; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
      ui(2*r+4*s+k-1, m+k) = -1.0; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
    }
  } else if(r == 0) {
    for(int k = 0; k < s; k++) {
      ui(2*r+2*k, m+r+k) = -1.0;
      ui(2*r+2*k+1, m+r+k) = 1.0;
      ui(2*r+2*s+2*k, m+r+s+k) = -1.0;
      ui(2*r+2*s+2*k+1, m+r+s+k) = 1.0;
    }
  } else {
    for(int k = 0; k < r; k++) {
      ui(2*k, m+k) = -1.0;
      ui(2*k+1, m+k) = 1.0;
    }
    for(int k = 0; k < s; k++) {
      ui(2*r+2*k, m+r+k) = -1.0;
      ui(2*r+2*k+1, m+r+k) = 1.0;
      ui(2*r+2*s+2*k, m+r+s+k) = -1.0;
      ui(2*r+2*s+2*k+1, m+r+s+k) = 1.0;
    }
    for(int k = 1; k < r; k++) {
      ui(2*r+4*s+k-1, m+k-1) = 1.0; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
      ui(2*r+4*s+k-1, m+k) = -1.0; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
    }
  }
  ui(dim_constr_lambda, dimchi-2) = -1.0;
  ui(dim_constr_lambda+1, dimchi-2) = 1.0;
  ui(dim_constr_lambda+2, dimchi-1) = -1.0;
  ui(dim_constr_lambda+3, dimchi-1) = 1.0;
  ui(dim_constr_lambda+4, dimchi-2) = -1.0; // \beta_1 + \beta_2 <= 0.999
  ui(dim_constr_lambda+4, dimchi-1) = -1.0; // \beta_1 + \beta_2 <= 0.999
  for(int k = 0; k < dimruline; k++) {
    ui(dim_constr_lambda+5+2*k, dimchi+k) = -1.0;
    ui(dim_constr_lambda+5+2*k+1, dimchi+k) = 1.0;
  }
  return ui;
}

// [[Rcpp::export]]
arma::vec fc_general_linear_constr_ci(int m, int r, int s) {
  // linear inequality constraints, ui %*% theta - ci >=0
  int dimruline = m*(m-1)/2;
  arma::vec dim_inequality_lambda_choices(2); dim_inequality_lambda_choices(0) = r-1; dim_inequality_lambda_choices(1) = 0;
  int dim_inequality_lambda = max(dim_inequality_lambda_choices);
  int dim_constr_lambda = 2*r + 4*s + dim_inequality_lambda;
  int dim_constr_beta12 = 5;
  int dim_constr_ruline = 2*dimruline;
  arma::mat ci_mat(dim_constr_lambda+dim_constr_beta12+dim_constr_ruline, 1); ci_mat.fill(0.0);
  arma::vec ci = ci_mat.col(0);
  if(s == 0) {
    for(int k = 0; k < r; k++) {
      ci(2*k) = -0.999;
      ci(2*k+1) = -0.999;
    }
    for(int k = 1; k < r; k++) {
      ci(2*r+4*s+k-1) = 0.001; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
    }
  } else if(r == 0) {
    for(int k = 0; k < s; k++) {
      ci(2*r+2*k) = -0.999;
      ci(2*r+2*k+1) = 0.001;
      ci(2*r+2*s+2*k) = -3.141;
      ci(2*r+2*s+2*k+1) = 0.001;
    }
  } else {
    for(int k = 0; k < r; k++) {
      ci(2*k) = -0.999;
      ci(2*k+1) = -0.999;
    }
    for(int k = 0; k < s; k++) {
      ci(2*r+2*k) = -0.999;
      ci(2*r+2*k+1) = 0.001;
      ci(2*r+2*s+2*k) = -3.141;
      ci(2*r+2*s+2*k+1) = 0.001;
    }
    for(int k = 1; k < r; k++) {
      ci(2*r+4*s+k-1) = 0.001; // lambda_{k} - \lambda_{k+1} - 0.001 >= 0
    }
  }
  ci(dim_constr_lambda) = -0.999;
  ci(dim_constr_lambda+1) = 0.001;
  ci(dim_constr_lambda+2) = -0.999;
  ci(dim_constr_lambda+3) = 0.001;
  ci(dim_constr_lambda+4) = -0.999; // \beta_1 + \beta_2 <= 0.999
  for(int k = 0; k < dimruline; k++) {
    ci(dim_constr_lambda+5+2*k) = -0.999;
    ci(dim_constr_lambda+5+2*k+1) = -0.999;
  }
  return ci;
}

////////////////////////////////////////////////// loss function //////////////////////////////////////////////////

// [[Rcpp::export]]
double fc_general_mathcalLtilde_n(int n, int m, int r, int s, int Bk, arma::vec theta_vec, arma::mat y_m_n) {
  // \widetilde{\mathcal{L}}_{n}(\bm\theta), under initial values \widetilde{y}_s = 1_m for s <= 0
  int dimkappa = r + 2*s + r*m*m + 2*s*m*m;
  int dimbeta = 2 + m*(m-1)/2;
  arma::vec omegauline = theta_vec.subvec(0, m-1); // \underline{\bm\omega}
  arma::vec kappa_vec = theta_vec.subvec(m, m+dimkappa-1); // \bm\kappa
  arma::vec beta_vec = theta_vec.subvec(m+dimkappa, m+dimkappa+dimbeta-1); // \bm\beta
  double beta_1 = beta_vec(0); double beta_2 = beta_vec(1); // \beta_1 and \beta_2
  arma::vec ruline = beta_vec.tail(dimbeta-2); // \underline\bm{r}
  arma::mat Ruline = fc_Ruline(m, ruline); // \underline{R}
  // summations using FFT algorithm
  int max_r_1 = r; if(max_r_1 == 0) {max_r_1 = 1;}
  int max_s_1 = s; if(max_s_1 == 0) {max_s_1 = 1;}
  arma::cube sum_lambdaklnyuline_0_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_0_m_nminus1_r.fill(0.0);
  arma::cube sum_gammaklnyuline_01_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_01_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_02_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_02_m_nminus1_s.fill(0.0);
  if(s == 0) {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
  } else if(r == 0) {
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  } else {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  }
  // \widetilde{\ell}_{t}(\bm\theta) and \widetilde\mathcal{L}_n(\bm\theta)
  double sum_ltilde = 0.0; // \sum_{t=1}^{n} \widetilde{\ell}_{t}(\bm\theta)
  // t=1
  arma::vec y_1 = y_m_n.col(0); // \bm{y}_{1}
  arma::mat Dtilde_1 = fc_Dtilde_1(omegauline); // \widetilde{D}_{1}(\bm\delta)
  arma::mat inverse_Dtilde_1 = fc_inverse_Dtilde_1(omegauline); // \widetilde{D}_{1}^{-1}(\bm\delta)
  arma::mat Rtilde_1 = fc_Rtilde_1(m, beta_1, beta_2, Ruline); // \widetilde{R}_{1}(\bm\theta)
  arma::mat Htilde_1 = fc_Htilde_t(Dtilde_1, Rtilde_1); // \widetilde{H}_{1}(\bm\theta)
  arma::vec varepsilontilde_1 = inverse_Dtilde_1 * y_1; // \widetilde{\varepsilon}_{1}(\bm\delta) = \widetilde{D}_{1}^{-1}(\bm\delta) \bm{y}_{1}
  // \widetilde{\ell}_1(\bm\theta)
  double ltilde_1 = fc_ltilde_t(y_1, Htilde_1); // \widetilde{\ell}_{1}(\bm\theta)
  sum_ltilde = sum_ltilde + ltilde_1;
  // t=2,...,n
  // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
  arma::mat varepsilontilde_tminusBkminus1TOtminus2_mat(m, Bk); varepsilontilde_tminusBkminus1TOtminus2_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta))
  varepsilontilde_tminusBkminus1TOtminus2_mat.each_col() = exp(-0.5*omegauline); // for s <= 0, \widetilde{\varepsilon}_s(\bm\delta) = (\exp{-1/2*\underliner{\omega}_1}, ..., \exp{-1/2*\underliner{\omega}_m})'
  arma::vec varepsilontilde_tminus1 = varepsilontilde_1; // \widetilde{\varepsilon}_{t-1}(\bm\delta)
  // \widetilde{R}_{t-1}(\bm\theta)
  arma::mat Rtilde_tminus1 = Rtilde_1; // \widetilde{R}_{t-1}(\bm\theta)
  for(int t = 2; t < (n+1); t++) {
    // \bm{y}_2, ..., \bm{y}_n
    arma::vec y_t = y_m_n.col(t-1); // \bm{y}_{t}
    // \widetilde{D}_2(\bm\delta), ..., \widetilde{D}_n(\bm\delta)
    arma::vec lnhulinetilde_t = fc_general_lnhulinetilde_t_general(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s); // \ln\widetilde\underline{\bm{h}}_t(\bm\delta)
    arma::mat Dtilde_t = fc_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}(\bm\delta)
    arma::mat inverse_Dtilde_t = fc_inverse_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}^{-1}(\bm\delta)
    // \widetilde{R}_2(\bm\theta), ..., \widetilde{R}_n(\bm\theta)
    arma::vec varepsilontilde_t = inverse_Dtilde_t * y_t; // \widetilde{\varepsilon}_{t}(\bm\delta) = \widetilde{D}_{t}^{-1}(\bm\delta) \bm{y}_{t}
    arma::mat varepsilontilde_tminusBkTOtminus1_mat(m, Bk); varepsilontilde_tminusBkTOtminus1_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk}(\bm\delta), ..., \widetilde{\varepsilon}_{t-1}(\bm\delta))
    varepsilontilde_tminusBkTOtminus1_mat.cols(0, Bk-2) = varepsilontilde_tminusBkminus1TOtminus2_mat.cols(1, Bk-1);
    varepsilontilde_tminusBkTOtminus1_mat.col(Bk-1) = varepsilontilde_tminus1;
    arma::mat Psitilde_tminus1 = fc_Psi(m, varepsilontilde_tminusBkTOtminus1_mat); // \widetilde{\Psi}_{t-1}(\bm\delta)
    arma::mat Rtilde_t = fc_Rtilde_t(m, Bk, beta_1, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1); // \widetilde{R}_{t}(\bm\theta)
    // \widetilde{H}_2(\bm\theta), ..., \widetilde{H}_n(\bm\theta)
    arma::mat Htilde_t = fc_Htilde_t(Dtilde_t, Rtilde_t); // \widetilde{H}_{t}(\bm\theta)
    // \widetilde{\ell}_2(\bm\theta), ..., \widetilde{\ell}_n(\bm\theta)
    double ltilde_t = fc_ltilde_t(y_t, Htilde_t); // \widetilde{\ell}_{t}(\bm\theta)
    sum_ltilde = sum_ltilde + ltilde_t;
    // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
    varepsilontilde_tminusBkminus1TOtminus2_mat = varepsilontilde_tminusBkTOtminus1_mat;
    varepsilontilde_tminus1 = varepsilontilde_t;
    // \widetilde{R}_{t-1}(\bm\theta)
    Rtilde_tminus1 = Rtilde_t;
  }
  double mathcalLtilde_n = sum_ltilde / (n * 1.0); // \widetilde\mathcal{L}_n(\bm\theta)
  return mathcalLtilde_n;
}

////////////////////////////////////////////////// gradient //////////////////////////////////////////////////

arma::vec fc_dltilde_dtheta_t(int m, int dimdelta, int dimbeta, int dimtheta, arma::vec y_t, arma::mat Htilde_t, arma::mat Dtilde_t, arma::mat Rtilde_t, arma::mat dlnhulinetilde_ddelta_t_mat, arma::cube dRtilde_ddelta_t_cube, arma::cube dRtilde_dbeta_t_cube) {
  // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
  arma::vec dltilde_dtheta_t(dimtheta); // \partial\widetilde{\ell}_t(\bm\theta) / \partial\theta_\ell = \tr[(I_m - \widetilde{H}_t^{-1}(\bm\theta) \bm{y}_t \bm{y}_t') \widetilde{H}_t^{-1}(\bm\theta) \partial\widetilde{H}_t(\bm\theta)/\partial\theta_\ell]
  arma::mat I_m(m,m); I_m.eye(m,m);
  arma::mat y_t_mat = fc_asmat(y_t, m, 1);
  arma::mat inverse_Htilde_t = Htilde_t.i();
  arma::mat left_mat = (I_m - inverse_Htilde_t * y_t * y_t_mat.t()) * inverse_Htilde_t;
  for(int l = 0; l < dimdelta; l++) {
    arma::vec dlnhulinetilde_ddelta_t_l = dlnhulinetilde_ddelta_t_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
    arma::mat dDtilde_ddelta_t_l = fc_dDtilde_ddelta_t_l(Dtilde_t, dlnhulinetilde_ddelta_t_l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
    arma::mat dRtilde_ddelta_t_l = dRtilde_ddelta_t_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
    arma::mat dHtilde_ddelta_t_l = fc_dHtilde_ddelta_t_l(Dtilde_t, Rtilde_t, dDtilde_ddelta_t_l, dRtilde_ddelta_t_l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
    dltilde_dtheta_t(l) = trace(left_mat * dHtilde_ddelta_t_l); // \partial\widetilde{\ell}_t(\bm\theta) / \partial\delta_\ell
  }
  for(int l = 0; l < dimbeta; l++) {
    arma::mat dRtilde_dbeta_t_l = dRtilde_dbeta_t_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_\ell
    arma::mat dHtilde_dbeta_t_l = fc_dHtilde_dbeta_t_l(Dtilde_t, dRtilde_dbeta_t_l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
    dltilde_dtheta_t(dimdelta+l) = trace(left_mat * dHtilde_dbeta_t_l); // \partial\widetilde{\ell}_t(\bm\theta) / \partial\beta_\ell
  }
  return dltilde_dtheta_t;
}

// [[Rcpp::export]]
arma::vec fc_dmathcalLtilde_dtheta(int n, int m, int r, int s, int Bk, arma::vec theta_vec, arma::mat y_m_n) { 
  // \partial\widetilde{\mathcal{L}}_{n}(\bm\theta) / \partial\bm\theta, under initial values \widetilde{y}_s = 1_m for s <= 0
  int dimkappa = r + 2*s + r*m*m + 2*s*m*m;
  int dimbeta = 2 + m*(m-1)/2;
  int dimdelta = m + dimkappa;
  int dimtheta = dimdelta + dimbeta;
  arma::vec omegauline = theta_vec.subvec(0, m-1); // \underline{\bm\omega}
  arma::vec kappa_vec = theta_vec.subvec(m, m+dimkappa-1); // \bm\kappa
  arma::vec beta_vec = theta_vec.subvec(m+dimkappa, m+dimkappa+dimbeta-1); // \bm\beta
  double beta_1 = beta_vec(0); double beta_2 = beta_vec(1); // \beta_1 and \beta_2
  arma::vec ruline = beta_vec.tail(dimbeta-2); // \underline\bm{r}
  arma::mat Ruline = fc_Ruline(m, ruline); // \underline{R}
  // summations using FFT algorithm
  int max_r_1 = r; if(max_r_1 == 0) {max_r_1 = 1;}
  int max_s_1 = s; if(max_s_1 == 0) {max_s_1 = 1;}
  arma::cube sum_lambdaklnyuline_0_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_0_m_nminus1_r.fill(0.0);
  arma::cube sum_lambdaklnyuline_1_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_1_m_nminus1_r.fill(0.0);
  arma::cube sum_gammaklnyuline_01_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_01_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_02_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_02_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_11_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_11_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_12_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_12_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_21_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_21_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_22_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_22_m_nminus1_s.fill(0.0);
  if(s == 0) {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
  } else if(r == 0) {
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  } else {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  }
  // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
  arma::mat sum_dltilde_dtheta_mat(dimtheta, 1); sum_dltilde_dtheta_mat.fill(0.0);
  arma::vec sum_dltilde_dtheta = sum_dltilde_dtheta_mat.col(0); // \sum_{t=1}^{n} \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
  // t=1
  arma::vec y_1 = y_m_n.col(0); // \bm{y}_{1}
  arma::mat Dtilde_1 = fc_Dtilde_1(omegauline); // \widetilde{D}_{1}(\bm\delta)
  arma::mat inverse_Dtilde_1 = fc_inverse_Dtilde_1(omegauline); // \widetilde{D}_{1}^{-1}(\bm\delta)
  arma::mat Rtilde_1 = fc_Rtilde_1(m, beta_1, beta_2, Ruline); // \widetilde{R}_{1}(\bm\theta)
  arma::mat Htilde_1 = fc_Htilde_t(Dtilde_1, Rtilde_1); // \widetilde{H}_{1}(\bm\theta)
  arma::vec varepsilontilde_1 = inverse_Dtilde_1 * y_1; // \widetilde{\varepsilon}_{1}(\bm\delta) = \widetilde{D}_{1}^{-1}(\bm\delta) \bm{y}_{1}
  // derivatives
  arma::mat dlnhulinetilde_ddelta_1_mat = fc_general_dlnhulinetilde_ddelta_1_mat(m, dimdelta); // \partial\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\bm\delta'_\ell, which is equal to \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'_\ell for t <= 0
  arma::cube dRtilde_ddelta_1_cube = fc_dRtilde_ddelta_1_cube(m, dimdelta); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\delta_\ell
  arma::cube dRuline_dbeta_cube = fc_dRuline_dbeta_cube(m, dimbeta); // the \ell-th slice is \partial\underline{R} / \partial\beta_\ell
  arma::cube dRtilde_dbeta_1_cube = fc_dRtilde_dbeta_1_cube(m, dimbeta, beta_1, beta_2, Ruline, dRuline_dbeta_cube); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\beta_\ell
  // \partial\widetilde{\ell}_1(\bm\theta)/\partial\bm\theta
  arma::vec dltilde_dtheta_1 = fc_dltilde_dtheta_t(m, dimdelta, dimbeta, dimtheta, y_1, Htilde_1, Dtilde_1, Rtilde_1, dlnhulinetilde_ddelta_1_mat, dRtilde_ddelta_1_cube, dRtilde_dbeta_1_cube); // \partial\widetilde{\ell}_{1}(\bm\theta) / \partial\bm\theta
  sum_dltilde_dtheta = sum_dltilde_dtheta + dltilde_dtheta_1;
  // t=2,...,n
  // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
  arma::mat varepsilontilde_tminusBkminus1TOtminus2_mat(m, Bk); varepsilontilde_tminusBkminus1TOtminus2_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta))
  varepsilontilde_tminusBkminus1TOtminus2_mat.each_col() = exp(-0.5*omegauline); // for s <= 0, \widetilde{\varepsilon}_s(\bm\delta) = (\exp{-1/2*\underliner{\omega}_1}, ..., \exp{-1/2*\underliner{\omega}_m})'
  arma::vec varepsilontilde_tminus1 = varepsilontilde_1; // \widetilde{\varepsilon}_{t-1}(\bm\delta)
  // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
  arma::cube dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.fill(0.0); // the i-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-2+i}(\bm\delta) / \partial\bm\delta'_\ell
  dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.each_slice() = dlnhulinetilde_ddelta_1_mat;
  arma::mat dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_1_mat;
  // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell and \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell
  arma::mat Rtilde_tminus1 = Rtilde_1; // \widetilde{R}_{t-1}(\bm\theta)
  arma::cube dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
  arma::cube dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_\ell
  for(int t = 2; t < (n+1); t++) {
    arma::vec y_t = y_m_n.col(t-1); // \bm{y}_{t}
    arma::vec lnhulinetilde_t = fc_general_lnhulinetilde_t_general(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s); // \ln\widetilde\underline{\bm{h}}_t(\bm\delta)
    arma::mat Dtilde_t = fc_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}(\bm\delta)
    arma::mat inverse_Dtilde_t = fc_inverse_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}^{-1}(\bm\delta)
    arma::vec varepsilontilde_t = inverse_Dtilde_t * y_t; // \widetilde{\varepsilon}_{t}(\bm\delta) = \widetilde{D}_{t}^{-1}(\bm\delta) \bm{y}_{t}
    arma::mat varepsilontilde_tminusBkTOtminus1_mat(m, Bk); varepsilontilde_tminusBkTOtminus1_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk}(\bm\delta), ..., \widetilde{\varepsilon}_{t-1}(\bm\delta))
    varepsilontilde_tminusBkTOtminus1_mat.cols(0, Bk-2) = varepsilontilde_tminusBkminus1TOtminus2_mat.cols(1, Bk-1);
    varepsilontilde_tminusBkTOtminus1_mat.col(Bk-1) = varepsilontilde_tminus1;
    arma::mat Psitilde_tminus1 = fc_Psi(m, varepsilontilde_tminusBkTOtminus1_mat); // \widetilde{\Psi}_{t-1}(\bm\delta)
    arma::mat Rtilde_t = fc_Rtilde_t(m, Bk, beta_1, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1); // \widetilde{R}_{t}(\bm\theta)
    arma::mat Htilde_t = fc_Htilde_t(Dtilde_t, Rtilde_t); // \widetilde{H}_{t}(\bm\theta)
    // derivatives
    arma::mat dlnhulinetilde_ddelta_t_mat = fc_general_dlnhulinetilde_ddelta_t_mat_general(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_lambdaklnyuline_1_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'
    arma::cube dlnhulinetilde_ddelta_tminusBkTOtminus1_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.fill(0.0); // the k-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-k)}(\bm\delta) / \partial\bm\delta'
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slices(0, Bk-2) = dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.slices(1, Bk-1);
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slice(Bk-1) = dlnhulinetilde_ddelta_tminus1_mat;
    arma::cube dRtilde_ddelta_t_cube = fc_dRtilde_ddelta_t_cube(m, dimdelta, beta_1, beta_2, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dlnhulinetilde_ddelta_tminusBkTOtminus1_cube, dRtilde_ddelta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
    arma::cube dRtilde_dbeta_t_cube = fc_dRtilde_dbeta_t_cube(m, dimbeta, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1, dRtilde_dbeta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\beta_\ell
    // \partial\widetilde{\ell}_2(\bm\theta)/\partial\bm\theta, ..., \partial\widetilde{\ell}_n(\bm\theta)/\partial\bm\theta
    arma::vec dltilde_dtheta_t = fc_dltilde_dtheta_t(m, dimdelta, dimbeta, dimtheta, y_t, Htilde_t, Dtilde_t, Rtilde_t, dlnhulinetilde_ddelta_t_mat, dRtilde_ddelta_t_cube, dRtilde_dbeta_t_cube); // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
    sum_dltilde_dtheta = sum_dltilde_dtheta + dltilde_dtheta_t;
    // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
    varepsilontilde_tminusBkminus1TOtminus2_mat = varepsilontilde_tminusBkTOtminus1_mat;
    varepsilontilde_tminus1 = varepsilontilde_t;
    // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
    dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube; 
    dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_t_mat; 
    // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell and \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell
    Rtilde_tminus1 = Rtilde_t;
    dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_t_cube;
    dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_t_cube;
  }
  arma::vec dmathcalLtilde_dtheta = sum_dltilde_dtheta / (n * 1.0); // \partial\widetilde\mathcal{L}_n(\bm\theta) / \partial\bm\theta
  return dmathcalLtilde_dtheta;
}

////////////////////////////////////////////////// ASD //////////////////////////////////////////////////

// [[Rcpp::export]]
arma::vec fc_general_ASD_thetahat(int n, int m, int r, int s, int Bk, arma::vec thetahat, arma::mat y_m_n) {
  // ASD of \widehat\bm\theta
  int dimkappa = r + 2*s + r*m*m + 2*s*m*m;
  int dimbeta = 2 + m*(m-1)/2;
  int dimdelta = m + dimkappa;
  int dimtheta = dimdelta + dimbeta;
  arma::vec omegauline = thetahat.subvec(0, m-1); // \underline{\bm\omega}
  arma::vec kappa_vec = thetahat.subvec(m, m+dimkappa-1); // \bm\kappa
  arma::vec beta_vec = thetahat.subvec(m+dimkappa, m+dimkappa+dimbeta-1); // \bm\beta
  double beta_1 = beta_vec(0); double beta_2 = beta_vec(1); // \beta_1 and \beta_2
  arma::vec ruline = beta_vec.tail(dimbeta-2); // \underline\bm{r}
  arma::mat Ruline = fc_Ruline(m, ruline); // \underline{R}
  // summations using FFT algorithm
  int max_r_1 = r; if(max_r_1 == 0) {max_r_1 = 1;}
  int max_s_1 = s; if(max_s_1 == 0) {max_s_1 = 1;}
  arma::cube sum_lambdaklnyuline_0_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_0_m_nminus1_r.fill(0.0);
  arma::cube sum_lambdaklnyuline_1_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_1_m_nminus1_r.fill(0.0);
  arma::cube sum_lambdaklnyuline_2_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_2_m_nminus1_r.fill(0.0);
  arma::cube sum_gammaklnyuline_01_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_01_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_02_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_02_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_11_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_11_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_12_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_12_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_21_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_21_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_22_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_22_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_31_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_31_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_32_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_32_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_41_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_41_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_42_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_42_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_51_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_51_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_52_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_52_m_nminus1_s.fill(0.0);
  if(s == 0) {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_2_m_nminus1_r = fc_sum_lambdaklnyuline_2_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
  } else if(r == 0) {
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_31_m_nminus1_s = fc_sum_gammaklnyuline_31_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_32_m_nminus1_s = fc_sum_gammaklnyuline_32_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_41_m_nminus1_s = fc_sum_gammaklnyuline_41_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_42_m_nminus1_s = fc_sum_gammaklnyuline_42_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_51_m_nminus1_s = fc_sum_gammaklnyuline_51_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_52_m_nminus1_s = fc_sum_gammaklnyuline_52_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  } else {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_2_m_nminus1_r = fc_sum_lambdaklnyuline_2_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_31_m_nminus1_s = fc_sum_gammaklnyuline_31_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_32_m_nminus1_s = fc_sum_gammaklnyuline_32_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_41_m_nminus1_s = fc_sum_gammaklnyuline_41_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_42_m_nminus1_s = fc_sum_gammaklnyuline_42_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_51_m_nminus1_s = fc_sum_gammaklnyuline_51_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_52_m_nminus1_s = fc_sum_gammaklnyuline_52_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  }
  // n \widehat{\Sigma}
  arma::mat nSigmahat(dimtheta, dimtheta); nSigmahat.fill(0.0); // n \widehat{\Sigma} = \sum_{t=1}^{n} \widehat{\Sigma}_t = \sum_{t=1}^{n} \partial\widetilde{\ell}_t(\widehat\bm\theta)/\partial\bm\theta \partial\widetilde{\ell}_t(\widehat\bm\theta)/\partial\bm\theta'
  // n \widehat{\Sigma}_{*}
  arma::mat nSigmaStarhat(dimtheta, dimtheta); nSigmaStarhat.fill(0.0); // n \widehat{\Sigma}_{*} = \sum_{t=1}^{n} \widehat{\Sigma}_{*, t} = \sum_{t=1}^{n} \partial^2\widetilde{\ell}_t(\widehat\bm\theta) / \partial\bm\theta\partial\bm\theta'
  // t=1
  arma::vec y_1 = y_m_n.col(0); // \bm{y}_{1}
  arma::mat Dtilde_1 = fc_Dtilde_1(omegauline); // \widetilde{D}_{1}(\bm\delta)
  arma::mat inverse_Dtilde_1 = fc_inverse_Dtilde_1(omegauline); // \widetilde{D}_{1}^{-1}(\bm\delta)
  arma::mat Rtilde_1 = fc_Rtilde_1(m, beta_1, beta_2, Ruline); // \widetilde{R}_{1}(\bm\theta)
  arma::mat Htilde_1 = fc_Htilde_t(Dtilde_1, Rtilde_1); // \widetilde{H}_{1}(\bm\theta)
  arma::vec varepsilontilde_1 = inverse_Dtilde_1 * y_1; // \widetilde{\varepsilon}_{1}(\bm\delta) = \widetilde{D}_{1}^{-1}(\bm\delta) \bm{y}_{1}
  // first derivatives
  arma::mat dlnhulinetilde_ddelta_1_mat = fc_general_dlnhulinetilde_ddelta_1_mat(m, dimdelta); // \partial\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\bm\delta'_\ell, which is equal to \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'_\ell for t <= 0
  arma::cube dDtilde_ddelta_1_cube = fc_dDtilde_ddelta_t_cube(m, dimdelta, Dtilde_1, dlnhulinetilde_ddelta_1_mat); // the \ell-th slice is \partial\widetilde{D}_1(\bm\delta) / \partial\delta_\ell
  arma::cube dRtilde_ddelta_1_cube = fc_dRtilde_ddelta_1_cube(m, dimdelta); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\delta_\ell
  arma::cube dRuline_dbeta_cube = fc_dRuline_dbeta_cube(m, dimbeta); // the \ell-th slice is \partial\underline{R} / \partial\beta_\ell
  arma::cube dRtilde_dbeta_1_cube = fc_dRtilde_dbeta_1_cube(m, dimbeta, beta_1, beta_2, Ruline, dRuline_dbeta_cube); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\beta_\ell
  arma::cube dHtilde_ddelta_1_cube = fc_dHtilde_ddelta_t_cube(m, dimdelta, Dtilde_1, Rtilde_1, dDtilde_ddelta_1_cube, dRtilde_ddelta_1_cube); // the \ell-th slice is \partial\widetilde{H}_1(\bm\theta) / \partial\delta_\ell
  arma::cube dHtilde_dbeta_1_cube = fc_dHtilde_dbeta_t_cube(m, dimbeta, Dtilde_1, dRtilde_dbeta_1_cube); // the \ell-th slice is \partial\widetilde{H}_1(\bm\theta) / \partial\beta_\ell
  // \widehat{\Sigma}_1
  arma::vec dltilde_dtheta_1 = fc2_dltilde_dtheta_t(m, dimdelta, dimbeta, dimtheta, y_1, Htilde_1, dHtilde_ddelta_1_cube, dHtilde_dbeta_1_cube); // \partial\widetilde{\ell}_1(\bm\theta) / \partial\bm\theta
  arma::mat dltilde_dtheta_1_mat = fc_asmat(dltilde_dtheta_1, dimtheta, 1);
  arma::mat Sigmahat_1 = dltilde_dtheta_1 * dltilde_dtheta_1_mat.t(); // \widehat{\Sigma}_1
  nSigmahat = nSigmahat + Sigmahat_1;
  // second derivatives
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_1_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_1_cube(m, dimdelta); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\delta_k\partial\delta_\ell, which is equal to \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell for t <= 0
  arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_1_field = fc_ddRtilde_ddeltaddeltaprime_1_field(m, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddRtilde_dbetadbetaprime_1_field = fc_ddRtilde_dbetadbetaprime_1_field(m, dimbeta, beta_1, beta_2, Ruline, dRuline_dbeta_cube); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::field<arma::mat> ddRtilde_dbetaddeltaprime_1_field = fc_ddRtilde_dbetaddeltaprime_1_field(m, dimdelta, dimbeta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\beta_k\partial\delta_\ell
  // \partial^2\ln\widetilde\underline{\bm{h}}_{1-\Bbbk}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{0}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix, of which the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{-\Bbbk+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  // \partial^2\widetilde{\ell}_{1}(\bm\theta) / \partial\bm\theta \partial\bm\theta'
  arma::mat ddltilde_dthetadthetaprime_1(dimtheta, dimtheta); ddltilde_dthetadthetaprime_1.fill(0.0);
  for(int k = 0; k < dimdelta; k++) {
    for(int l = 0; l < (k+1); l++) {
      arma::vec dlnhulinetilde_ddelta_1_l = dlnhulinetilde_ddelta_1_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
      arma::vec dlnhulinetilde_ddelta_1_k = dlnhulinetilde_ddelta_1_mat.col(k); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k
      arma::mat ddlnhulinetilde_ddeltaddeltaprime_1_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_1_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\bm\delta'
      arma::vec ddlnhulinetilde_ddeltaddeltaprime_1_kl = ddlnhulinetilde_ddeltaddeltaprime_1_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat(m, Bk); ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat.fill(0.0); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{-\Bbbk+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
      ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat.each_col() = ddlnhulinetilde_ddeltaddeltaprime_1_kl;
      ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(k, l) = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(l, k) = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat;
      arma::mat dDtilde_ddelta_1_l = dDtilde_ddelta_1_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
      arma::mat dDtilde_ddelta_1_k = dDtilde_ddelta_1_cube.slice(k); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_k
      arma::mat dRtilde_ddelta_1_l = dRtilde_ddelta_1_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dRtilde_ddelta_1_k = dRtilde_ddelta_1_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_k
      arma::mat dHtilde_ddelta_1_l = dHtilde_ddelta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dHtilde_ddelta_1_k = dHtilde_ddelta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_k
      arma::mat ddDtilde_ddeltaddeltaprime_1_kl = fc_ddDtilde_ddeltaddeltaprime_t_kl(dlnhulinetilde_ddelta_1_l, dlnhulinetilde_ddelta_1_k, ddlnhulinetilde_ddeltaddeltaprime_1_kl, Dtilde_1); // \partial^2\widetilde{D}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddRtilde_ddeltaddeltaprime_1_kl = ddRtilde_ddeltaddeltaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddHtilde_ddeltaddeltaprime_1_kl = fc_ddHtilde_ddeltaddeltaprime_t_kl(Dtilde_1, Rtilde_1, dDtilde_ddelta_1_l, dDtilde_ddelta_1_k, dRtilde_ddelta_1_l, dRtilde_ddelta_1_k, ddDtilde_ddeltaddeltaprime_1_kl, ddRtilde_ddeltaddeltaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
      ddltilde_dthetadthetaprime_1(k,l) = ddltilde_dthetadthetaprime_1(l,k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_ddelta_1_l, dHtilde_ddelta_1_k, ddHtilde_ddeltaddeltaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\delta_k\partial\delta_\ell
    }
  }
  for(int k = 0; k < dimbeta; k++) {
    for(int l = 0; l < (k+1); l++) {
      arma::mat dHtilde_dbeta_1_l = dHtilde_dbeta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
      arma::mat dHtilde_dbeta_1_k = dHtilde_dbeta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
      arma::mat ddRtilde_dbetadbetaprime_1_kl = ddRtilde_dbetadbetaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
      arma::mat ddHtilde_dbetadbetaprime_1_kl = fc_ddHtilde_dbetadbetaprime_t_kl(Dtilde_1, ddRtilde_dbetadbetaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
      ddltilde_dthetadthetaprime_1(dimdelta+k,dimdelta+l) = ddltilde_dthetadthetaprime_1(dimdelta+l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_dbeta_1_l, dHtilde_dbeta_1_k, ddHtilde_dbetadbetaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\beta_\ell
    }
  }
  for(int k = 0; k < dimbeta; k++) {
    for(int l = 0; l < dimdelta; l++) {
      arma::mat dDtilde_ddelta_1_l = dDtilde_ddelta_1_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
      arma::mat dRtilde_dbeta_1_k = dRtilde_dbeta_1_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_k
      arma::mat dHtilde_ddelta_1_l = dHtilde_ddelta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dHtilde_dbeta_1_k = dHtilde_dbeta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
      arma::mat ddRtilde_dbetaddeltaprime_1_kl = ddRtilde_dbetaddeltaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
      arma::mat ddHtilde_dbetaddeltaprime_1_kl = fc_ddHtilde_dbetaddeltaprime_t_kl(Dtilde_1, dDtilde_ddelta_1_l, dRtilde_dbeta_1_k, ddRtilde_dbetaddeltaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
      ddltilde_dthetadthetaprime_1(dimdelta+k,l) = ddltilde_dthetadthetaprime_1(l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_ddelta_1_l, dHtilde_dbeta_1_k, ddHtilde_dbetaddeltaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\delta_\ell
    }
  }
  // \widehat{\Sigma}_{*, 1}
  arma::mat SigmaStarhat_1 = ddltilde_dthetadthetaprime_1; // \widehat{\Sigma}_{*, 1}
  nSigmaStarhat = nSigmaStarhat + SigmaStarhat_1;
  // t=2,...,n
  // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
  arma::mat varepsilontilde_tminusBkminus1TOtminus2_mat(m, Bk); varepsilontilde_tminusBkminus1TOtminus2_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta))
  varepsilontilde_tminusBkminus1TOtminus2_mat.each_col() = exp(-0.5*omegauline); // for s <= 0, \widetilde{\varepsilon}_s(\bm\delta) = (\exp{-1/2*\underliner{\omega}_1}, ..., \exp{-1/2*\underliner{\omega}_m})'
  arma::vec varepsilontilde_tminus1 = varepsilontilde_1; // \widetilde{\varepsilon}_{t-1}(\bm\delta)
  // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
  arma::cube dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.fill(0.0); // the i-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-2+i}(\bm\delta) / \partial\bm\delta'_\ell
  dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.each_slice() = dlnhulinetilde_ddelta_1_mat;
  arma::mat dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_1_mat;
  // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_k\partial\delta_\ell, and \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field; // the (k, \ell, i)-th element is the vector \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-2+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube = ddlnhulinetilde_ddeltaddeltaprime_1_cube;
  // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell, \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell, and \partial^2\widetilde{R}_{t-1}(\bm\theta)/\partial\theta_\k\partial\theta_\ell
  arma::mat Rtilde_tminus1 = Rtilde_1; // \widetilde{R}_{t-1}(\bm\theta)
  arma::cube dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
  arma::cube dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_\ell
  arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_tminus1_field = ddRtilde_ddeltaddeltaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddRtilde_dbetadbetaprime_tminus1_field = ddRtilde_dbetadbetaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::field<arma::mat> ddRtilde_dbetaddeltaprime_tminus1_field = ddRtilde_dbetaddeltaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\delta_\ell
  for(int t = 2; t < (n+1); t++) {
    arma::vec y_t = y_m_n.col(t-1); // \bm{y}_{t}
    arma::vec lnhulinetilde_t = fc_general_lnhulinetilde_t_general(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s); // \ln\widetilde\underline{\bm{h}}_t(\bm\delta)
    arma::mat Dtilde_t = fc_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}(\bm\delta)
    arma::mat inverse_Dtilde_t = fc_inverse_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}^{-1}(\bm\delta)
    arma::vec varepsilontilde_t = inverse_Dtilde_t * y_t; // \widetilde{\varepsilon}_{t}(\bm\delta) = \widetilde{D}_{t}^{-1}(\bm\delta) \bm{y}_{t}
    arma::mat varepsilontilde_tminusBkTOtminus1_mat(m, Bk); varepsilontilde_tminusBkTOtminus1_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk}(\bm\delta), ..., \widetilde{\varepsilon}_{t-1}(\bm\delta))
    varepsilontilde_tminusBkTOtminus1_mat.cols(0, Bk-2) = varepsilontilde_tminusBkminus1TOtminus2_mat.cols(1, Bk-1);
    varepsilontilde_tminusBkTOtminus1_mat.col(Bk-1) = varepsilontilde_tminus1;
    arma::mat Psitilde_tminus1 = fc_Psi(m, varepsilontilde_tminusBkTOtminus1_mat); // \widetilde{\Psi}_{t-1}(\bm\delta)
    arma::mat Rtilde_t = fc_Rtilde_t(m, Bk, beta_1, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1); // \widetilde{R}_{t}(\bm\theta)
    arma::mat Htilde_t = fc_Htilde_t(Dtilde_t, Rtilde_t); // \widetilde{H}_{t}(\bm\theta)
    // first derivatives
    arma::mat dlnhulinetilde_ddelta_t_mat = fc_general_dlnhulinetilde_ddelta_t_mat_general(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_lambdaklnyuline_1_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'
    arma::cube dlnhulinetilde_ddelta_tminusBkTOtminus1_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.fill(0.0); // the k-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-k)}(\bm\delta) / \partial\bm\delta'
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slices(0, Bk-2) = dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.slices(1, Bk-1);
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slice(Bk-1) = dlnhulinetilde_ddelta_tminus1_mat;
    arma::cube dDtilde_ddelta_t_cube = fc_dDtilde_ddelta_t_cube(m, dimdelta, Dtilde_t, dlnhulinetilde_ddelta_t_mat); // the \ell-th slice is \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
    arma::cube dPsitilde_ddelta_tminus1_cube = fc_dPsitilde_ddelta_tminus1_cube(m, dimdelta, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dlnhulinetilde_ddelta_tminusBkTOtminus1_cube); // the \ell-th slice is \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
    arma::cube dRtilde_ddelta_t_cube = fc2_dRtilde_ddelta_t_cube(m, dimdelta, beta_1, beta_2, dPsitilde_ddelta_tminus1_cube, dRtilde_ddelta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
    arma::cube dRtilde_dbeta_t_cube = fc_dRtilde_dbeta_t_cube(m, dimbeta, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1, dRtilde_dbeta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\beta_\ell
    arma::cube dHtilde_ddelta_t_cube = fc_dHtilde_ddelta_t_cube(m, dimdelta, Dtilde_t, Rtilde_t, dDtilde_ddelta_t_cube, dRtilde_ddelta_t_cube); // the \ell-th slice is \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
    arma::cube dHtilde_dbeta_t_cube = fc_dHtilde_dbeta_t_cube(m, dimbeta, Dtilde_t, dRtilde_dbeta_t_cube); // the \ell-th slice is \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
    // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
    arma::vec dltilde_dtheta_t = fc2_dltilde_dtheta_t(m, dimdelta, dimbeta, dimtheta, y_t, Htilde_t, dHtilde_ddelta_t_cube, dHtilde_dbeta_t_cube); // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
    // \widehat{\Sigma}_t
    arma::mat dltilde_dtheta_t_mat = fc_asmat(dltilde_dtheta_t, dimtheta, 1);
    arma::mat Sigmahat_t = dltilde_dtheta_t * dltilde_dtheta_t_mat.t(); // \widehat{\Sigma}_t
    nSigmahat = nSigmahat + Sigmahat_t;
    // second derivatives
    arma::cube ddlnhulinetilde_ddeltaddeltaprime_t_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube_general(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_1_m_nminus1_r, sum_lambdaklnyuline_2_m_nminus1_r, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s, sum_gammaklnyuline_31_m_nminus1_s, sum_gammaklnyuline_32_m_nminus1_s, sum_gammaklnyuline_41_m_nminus1_s, sum_gammaklnyuline_42_m_nminus1_s, sum_gammaklnyuline_51_m_nminus1_s, sum_gammaklnyuline_52_m_nminus1_s); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_t_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddRtilde_dbetadbetaprime_t_field(dimbeta, dimbeta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
    arma::field<arma::mat> ddRtilde_dbetaddeltaprime_t_field(dimbeta, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
    // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix, of which the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta \partial\bm\theta'
    arma::mat ddltilde_dthetadthetaprime_t(dimtheta, dimtheta); ddltilde_dthetadthetaprime_t.fill(0.0);
    for(int k = 0; k < dimdelta; k++) {
      for(int l = 0; l < (k+1); l++) {
        arma::vec dlnhulinetilde_ddelta_t_l = dlnhulinetilde_ddelta_t_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
        arma::vec dlnhulinetilde_ddelta_t_k = dlnhulinetilde_ddelta_t_mat.col(k); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k
        arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.col(l); // the i-th column is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_\ell
        arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.col(k); // the i-th column is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_k
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_t_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\bm\delta'
        arma::vec ddlnhulinetilde_ddeltaddeltaprime_t_kl = ddlnhulinetilde_ddeltaddeltaprime_t_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field(k, l); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\bm\delta'
        arma::vec ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl = ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat(m, Bk); ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.fill(0.0); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.cols(0, Bk-2) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_kl_mat.cols(1, Bk-1);
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.col(Bk-1) = ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl;
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(k, l) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(l, k) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat;
        arma::mat dDtilde_ddelta_t_l = dDtilde_ddelta_t_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
        arma::mat dDtilde_ddelta_t_k = dDtilde_ddelta_t_cube.slice(k); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_k
        arma::mat dPsitilde_ddelta_tminus1_k = dPsitilde_ddelta_tminus1_cube.slice(k); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_k
        arma::mat dRtilde_ddelta_t_l = dRtilde_ddelta_t_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dRtilde_ddelta_t_k = dRtilde_ddelta_t_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_k
        arma::mat dHtilde_ddelta_t_l = dHtilde_ddelta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_ddelta_t_k = dHtilde_ddelta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_k
        arma::mat ddDtilde_ddeltaddeltaprime_t_kl = fc_ddDtilde_ddeltaddeltaprime_t_kl(dlnhulinetilde_ddelta_t_l, dlnhulinetilde_ddelta_t_k, ddlnhulinetilde_ddeltaddeltaprime_t_kl, Dtilde_t); // \partial^2\widetilde{D}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddRtilde_ddeltaddeltaprime_tminus1_kl = ddRtilde_ddeltaddeltaprime_tminus1_field(k,l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddRtilde_ddeltaddeltaprime_t_kl = fc_ddRtilde_ddeltaddeltaprime_t_kl(m, beta_1, beta_2, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dPsitilde_ddelta_tminus1_k, dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat, dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat, ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat, ddRtilde_ddeltaddeltaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
        ddRtilde_ddeltaddeltaprime_t_field(k, l) = ddRtilde_ddeltaddeltaprime_t_field(l, k) = ddRtilde_ddeltaddeltaprime_t_kl;
        arma::mat ddHtilde_ddeltaddeltaprime_t_kl = fc_ddHtilde_ddeltaddeltaprime_t_kl(Dtilde_t, Rtilde_t, dDtilde_ddelta_t_l, dDtilde_ddelta_t_k, dRtilde_ddelta_t_l, dRtilde_ddelta_t_k, ddDtilde_ddeltaddeltaprime_t_kl, ddRtilde_ddeltaddeltaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
        ddltilde_dthetadthetaprime_t(k,l) = ddltilde_dthetadthetaprime_t(l,k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_ddelta_t_l, dHtilde_ddelta_t_k, ddHtilde_ddeltaddeltaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\delta_k\partial\delta_\ell
      }
    }
    for(int k = 0; k < dimbeta; k++) {
      for(int l = 0; l < (k+1); l++) {
        arma::mat dRtilde_dbeta_tminus1_l = dRtilde_dbeta_tminus1_cube.slice(l); // \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_\ell
        arma::mat dHtilde_dbeta_t_l = dHtilde_dbeta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
        arma::mat dHtilde_dbeta_t_k = dHtilde_dbeta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
        arma::mat ddRtilde_dbetadbetaprime_tminus1_kl = ddRtilde_dbetadbetaprime_tminus1_field(k, l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\beta_\ell
        arma::mat ddRtilde_dbetadbetaprime_t_kl = fc_ddRtilde_dbetadbetaprime_t_kl(l, k, m, beta_2, dRtilde_dbeta_tminus1_l, ddRtilde_dbetadbetaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
        ddRtilde_dbetadbetaprime_t_field(k, l) = ddRtilde_dbetadbetaprime_t_field(l, k) = ddRtilde_dbetadbetaprime_t_kl;
        arma::mat ddHtilde_dbetadbetaprime_t_kl = fc_ddHtilde_dbetadbetaprime_t_kl(Dtilde_t, ddRtilde_dbetadbetaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
        ddltilde_dthetadthetaprime_t(dimdelta+k,dimdelta+l) = ddltilde_dthetadthetaprime_t(dimdelta+l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_dbeta_t_l, dHtilde_dbeta_t_k, ddHtilde_dbetadbetaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\beta_\ell
      }
    }
    for(int k = 0; k < dimbeta; k++) {
      for(int l = 0; l < dimdelta; l++) {
        arma::mat dDtilde_ddelta_t_l = dDtilde_ddelta_t_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
        arma::mat dPsitilde_ddelta_tminus1_l = dPsitilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
        arma::mat dRtilde_dbeta_t_k = dRtilde_dbeta_t_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_k
        arma::mat dRtilde_ddelta_tminus1_l = dRtilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_ddelta_t_l = dHtilde_ddelta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_dbeta_t_k = dHtilde_dbeta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
        arma::mat ddRtilde_dbetaddeltaprime_tminus1_kl = ddRtilde_dbetaddeltaprime_tminus1_field(k, l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\delta_\ell
        arma::mat ddRtilde_dbetaddeltaprime_t_kl = fc_ddRtilde_dbetaddeltaprime_t_kl(l, k, m, beta_2, dPsitilde_ddelta_tminus1_l, dRtilde_ddelta_tminus1_l, ddRtilde_dbetaddeltaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
        ddRtilde_dbetaddeltaprime_t_field(k, l) = ddRtilde_dbetaddeltaprime_t_kl;
        arma::mat ddHtilde_dbetaddeltaprime_t_kl = fc_ddHtilde_dbetaddeltaprime_t_kl(Dtilde_t, dDtilde_ddelta_t_l, dRtilde_dbeta_t_k, ddRtilde_dbetaddeltaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
        ddltilde_dthetadthetaprime_t(dimdelta+k,l) = ddltilde_dthetadthetaprime_t(l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_ddelta_t_l, dHtilde_dbeta_t_k, ddHtilde_dbetaddeltaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\delta_\ell
      }
    }
    // \widehat{\Sigma}_{*, t}
    arma::mat SigmaStarhat_t = ddltilde_dthetadthetaprime_t; // \widehat{\Sigma}_{*, t}
    nSigmaStarhat = nSigmaStarhat + SigmaStarhat_t;
    // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
    varepsilontilde_tminusBkminus1TOtminus2_mat = varepsilontilde_tminusBkTOtminus1_mat;
    varepsilontilde_tminus1 = varepsilontilde_t;
    // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
    dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube;
    dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_t_mat;
    // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_k\partial\delta_\ell and \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field;
    ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube = ddlnhulinetilde_ddeltaddeltaprime_t_cube;
    // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell, \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell, and \partial^2\widetilde{R}_{t-1}(\bm\theta)/\partial\theta_\k\partial\theta_\ell
    Rtilde_tminus1 = Rtilde_t;
    dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_t_cube;
    dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_t_cube;
    ddRtilde_ddeltaddeltaprime_tminus1_field = ddRtilde_ddeltaddeltaprime_t_field;
    ddRtilde_dbetadbetaprime_tminus1_field = ddRtilde_dbetadbetaprime_t_field;
    ddRtilde_dbetaddeltaprime_tminus1_field = ddRtilde_dbetaddeltaprime_t_field;
  }
  arma::mat Sigmahat = nSigmahat / (n * 1.0); // \widehat{\Sigma} = 1/n \sum_{t=1}^{n} \widehat{\Sigma}_t
  arma::mat SigmaStarhat = nSigmaStarhat / (n * 1.0); // \widehat{\Sigma}_{*} = 1/n \sum_{t=1}^{n} \widehat{\Sigma}_{*, t}
  arma::mat CovMat = SigmaStarhat.i() * Sigmahat * SigmaStarhat.i(); // \widehat\Sigma_\star^{-1} \widehat\Sigma \widehat\Sigma_\star^{-1}
  arma::vec var_sqrtnthetahat = CovMat.diag(0); // variances of \sqrt{n} \widehat\bm\theta
  arma::vec var_thetahat = var_sqrtnthetahat / (n * 1.0); // variances of \widehat\bm\theta
  arma::vec ASD_thetahat = sqrt(var_thetahat); // ASD of \widehat\bm\theta
  return ASD_thetahat;
}

////////////////////////////////////////////////// Spillover effect //////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat fc_transfmat_theta_to_vecPhi1(int m, int r, int s, int dimtheta) {
  // transformation matrix for trasfering \theta to vec(\Phi_1) with \Phi_1 = \sum_{k=1}^{r} G_{0,k} + \sum_{k=1}^{s} G_{1,k}
  arma::mat transfmat(m*m, dimtheta); transfmat.fill(0.0);
  int dim1 = m + r + 2*s;
  for(int i = 0; i < (m*m); i++) {
    for(int k = 0; k < (r+s); k++) {
      transfmat(i, dim1+i+k*m*m) = 1.0;
    }
  }
  return transfmat;
}

// [[Rcpp::export]]
arma::vec fc_lowrank_ASD_vecPhi1hat(int n, int m, int r, int s, int Bk, arma::vec thetahat, arma::mat y_m_n, arma::mat Deltahat_theta_vartheta) {
  // ASD of vec(\Phi_1) with \Phi_1 = \sum_{k=1}^{r} G_{0,k} + \sum_{k=1}^{s} G_{1,k}
  int dimkappa = r + 2*s + r*m*m + 2*s*m*m;
  int dimbeta = 2 + m*(m-1)/2;
  int dimdelta = m + dimkappa;
  int dimtheta = dimdelta + dimbeta;
  arma::vec omegauline = thetahat.subvec(0, m-1); // \underline{\bm\omega}
  arma::vec kappa_vec = thetahat.subvec(m, m+dimkappa-1); // \bm\kappa
  arma::vec beta_vec = thetahat.subvec(m+dimkappa, m+dimkappa+dimbeta-1); // \bm\beta
  double beta_1 = beta_vec(0); double beta_2 = beta_vec(1); // \beta_1 and \beta_2
  arma::vec ruline = beta_vec.tail(dimbeta-2); // \underline\bm{r}
  arma::mat Ruline = fc_Ruline(m, ruline); // \underline{R}
  // summations using FFT algorithm
  int max_r_1 = r; if(max_r_1 == 0) {max_r_1 = 1;}
  int max_s_1 = s; if(max_s_1 == 0) {max_s_1 = 1;}
  arma::cube sum_lambdaklnyuline_0_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_0_m_nminus1_r.fill(0.0);
  arma::cube sum_lambdaklnyuline_1_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_1_m_nminus1_r.fill(0.0);
  arma::cube sum_lambdaklnyuline_2_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_2_m_nminus1_r.fill(0.0);
  arma::cube sum_gammaklnyuline_01_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_01_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_02_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_02_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_11_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_11_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_12_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_12_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_21_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_21_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_22_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_22_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_31_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_31_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_32_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_32_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_41_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_41_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_42_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_42_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_51_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_51_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_52_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_52_m_nminus1_s.fill(0.0);
  if(s == 0) {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_2_m_nminus1_r = fc_sum_lambdaklnyuline_2_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
  } else if(r == 0) {
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_31_m_nminus1_s = fc_sum_gammaklnyuline_31_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_32_m_nminus1_s = fc_sum_gammaklnyuline_32_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_41_m_nminus1_s = fc_sum_gammaklnyuline_41_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_42_m_nminus1_s = fc_sum_gammaklnyuline_42_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_51_m_nminus1_s = fc_sum_gammaklnyuline_51_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_52_m_nminus1_s = fc_sum_gammaklnyuline_52_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  } else {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_2_m_nminus1_r = fc_sum_lambdaklnyuline_2_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_31_m_nminus1_s = fc_sum_gammaklnyuline_31_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_32_m_nminus1_s = fc_sum_gammaklnyuline_32_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_41_m_nminus1_s = fc_sum_gammaklnyuline_41_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_42_m_nminus1_s = fc_sum_gammaklnyuline_42_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_51_m_nminus1_s = fc_sum_gammaklnyuline_51_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_52_m_nminus1_s = fc_sum_gammaklnyuline_52_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  }
  // n \widehat{\Sigma}
  arma::mat nSigmahat(dimtheta, dimtheta); nSigmahat.fill(0.0); // n \widehat{\Sigma} = \sum_{t=1}^{n} \widehat{\Sigma}_t = \sum_{t=1}^{n} \partial\widetilde{\ell}_t(\widehat\bm\theta)/\partial\bm\theta \partial\widetilde{\ell}_t(\widehat\bm\theta)/\partial\bm\theta'
  // n \widehat{\Sigma}_{*}
  arma::mat nSigmaStarhat(dimtheta, dimtheta); nSigmaStarhat.fill(0.0); // n \widehat{\Sigma}_{*} = \sum_{t=1}^{n} \widehat{\Sigma}_{*, t} = \sum_{t=1}^{n} \partial^2\widetilde{\ell}_t(\widehat\bm\theta) / \partial\bm\theta\partial\bm\theta'
  // t=1
  arma::vec y_1 = y_m_n.col(0); // \bm{y}_{1}
  arma::mat Dtilde_1 = fc_Dtilde_1(omegauline); // \widetilde{D}_{1}(\bm\delta)
  arma::mat inverse_Dtilde_1 = fc_inverse_Dtilde_1(omegauline); // \widetilde{D}_{1}^{-1}(\bm\delta)
  arma::mat Rtilde_1 = fc_Rtilde_1(m, beta_1, beta_2, Ruline); // \widetilde{R}_{1}(\bm\theta)
  arma::mat Htilde_1 = fc_Htilde_t(Dtilde_1, Rtilde_1); // \widetilde{H}_{1}(\bm\theta)
  arma::vec varepsilontilde_1 = inverse_Dtilde_1 * y_1; // \widetilde{\varepsilon}_{1}(\bm\delta) = \widetilde{D}_{1}^{-1}(\bm\delta) \bm{y}_{1}
  // first derivatives
  arma::mat dlnhulinetilde_ddelta_1_mat = fc_general_dlnhulinetilde_ddelta_1_mat(m, dimdelta); // \partial\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\bm\delta'_\ell, which is equal to \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'_\ell for t <= 0
  arma::cube dDtilde_ddelta_1_cube = fc_dDtilde_ddelta_t_cube(m, dimdelta, Dtilde_1, dlnhulinetilde_ddelta_1_mat); // the \ell-th slice is \partial\widetilde{D}_1(\bm\delta) / \partial\delta_\ell
  arma::cube dRtilde_ddelta_1_cube = fc_dRtilde_ddelta_1_cube(m, dimdelta); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\delta_\ell
  arma::cube dRuline_dbeta_cube = fc_dRuline_dbeta_cube(m, dimbeta); // the \ell-th slice is \partial\underline{R} / \partial\beta_\ell
  arma::cube dRtilde_dbeta_1_cube = fc_dRtilde_dbeta_1_cube(m, dimbeta, beta_1, beta_2, Ruline, dRuline_dbeta_cube); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\beta_\ell
  arma::cube dHtilde_ddelta_1_cube = fc_dHtilde_ddelta_t_cube(m, dimdelta, Dtilde_1, Rtilde_1, dDtilde_ddelta_1_cube, dRtilde_ddelta_1_cube); // the \ell-th slice is \partial\widetilde{H}_1(\bm\theta) / \partial\delta_\ell
  arma::cube dHtilde_dbeta_1_cube = fc_dHtilde_dbeta_t_cube(m, dimbeta, Dtilde_1, dRtilde_dbeta_1_cube); // the \ell-th slice is \partial\widetilde{H}_1(\bm\theta) / \partial\beta_\ell
  // \widehat{\Sigma}_1
  arma::vec dltilde_dtheta_1 = fc2_dltilde_dtheta_t(m, dimdelta, dimbeta, dimtheta, y_1, Htilde_1, dHtilde_ddelta_1_cube, dHtilde_dbeta_1_cube); // \partial\widetilde{\ell}_1(\bm\theta) / \partial\bm\theta
  arma::mat dltilde_dtheta_1_mat = fc_asmat(dltilde_dtheta_1, dimtheta, 1);
  arma::mat Sigmahat_1 = dltilde_dtheta_1 * dltilde_dtheta_1_mat.t(); // \widehat{\Sigma}_1
  nSigmahat = nSigmahat + Sigmahat_1;
  // second derivatives
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_1_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_1_cube(m, dimdelta); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\delta_k\partial\delta_\ell, which is equal to \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell for t <= 0
  arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_1_field = fc_ddRtilde_ddeltaddeltaprime_1_field(m, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddRtilde_dbetadbetaprime_1_field = fc_ddRtilde_dbetadbetaprime_1_field(m, dimbeta, beta_1, beta_2, Ruline, dRuline_dbeta_cube); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::field<arma::mat> ddRtilde_dbetaddeltaprime_1_field = fc_ddRtilde_dbetaddeltaprime_1_field(m, dimdelta, dimbeta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\beta_k\partial\delta_\ell
  // \partial^2\ln\widetilde\underline{\bm{h}}_{1-\Bbbk}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{0}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix, of which the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{-\Bbbk+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  // \partial^2\widetilde{\ell}_{1}(\bm\theta) / \partial\bm\theta \partial\bm\theta'
  arma::mat ddltilde_dthetadthetaprime_1(dimtheta, dimtheta); ddltilde_dthetadthetaprime_1.fill(0.0);
  for(int k = 0; k < dimdelta; k++) {
    for(int l = 0; l < (k+1); l++) {
      arma::vec dlnhulinetilde_ddelta_1_l = dlnhulinetilde_ddelta_1_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
      arma::vec dlnhulinetilde_ddelta_1_k = dlnhulinetilde_ddelta_1_mat.col(k); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k
      arma::mat ddlnhulinetilde_ddeltaddeltaprime_1_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_1_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\bm\delta'
      arma::vec ddlnhulinetilde_ddeltaddeltaprime_1_kl = ddlnhulinetilde_ddeltaddeltaprime_1_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat(m, Bk); ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat.fill(0.0); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{-\Bbbk+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
      ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat.each_col() = ddlnhulinetilde_ddeltaddeltaprime_1_kl;
      ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(k, l) = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(l, k) = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat;
      arma::mat dDtilde_ddelta_1_l = dDtilde_ddelta_1_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
      arma::mat dDtilde_ddelta_1_k = dDtilde_ddelta_1_cube.slice(k); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_k
      arma::mat dRtilde_ddelta_1_l = dRtilde_ddelta_1_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dRtilde_ddelta_1_k = dRtilde_ddelta_1_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_k
      arma::mat dHtilde_ddelta_1_l = dHtilde_ddelta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dHtilde_ddelta_1_k = dHtilde_ddelta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_k
      arma::mat ddDtilde_ddeltaddeltaprime_1_kl = fc_ddDtilde_ddeltaddeltaprime_t_kl(dlnhulinetilde_ddelta_1_l, dlnhulinetilde_ddelta_1_k, ddlnhulinetilde_ddeltaddeltaprime_1_kl, Dtilde_1); // \partial^2\widetilde{D}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddRtilde_ddeltaddeltaprime_1_kl = ddRtilde_ddeltaddeltaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddHtilde_ddeltaddeltaprime_1_kl = fc_ddHtilde_ddeltaddeltaprime_t_kl(Dtilde_1, Rtilde_1, dDtilde_ddelta_1_l, dDtilde_ddelta_1_k, dRtilde_ddelta_1_l, dRtilde_ddelta_1_k, ddDtilde_ddeltaddeltaprime_1_kl, ddRtilde_ddeltaddeltaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
      ddltilde_dthetadthetaprime_1(k,l) = ddltilde_dthetadthetaprime_1(l,k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_ddelta_1_l, dHtilde_ddelta_1_k, ddHtilde_ddeltaddeltaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\delta_k\partial\delta_\ell
    }
  }
  for(int k = 0; k < dimbeta; k++) {
    for(int l = 0; l < (k+1); l++) {
      arma::mat dHtilde_dbeta_1_l = dHtilde_dbeta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
      arma::mat dHtilde_dbeta_1_k = dHtilde_dbeta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
      arma::mat ddRtilde_dbetadbetaprime_1_kl = ddRtilde_dbetadbetaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
      arma::mat ddHtilde_dbetadbetaprime_1_kl = fc_ddHtilde_dbetadbetaprime_t_kl(Dtilde_1, ddRtilde_dbetadbetaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
      ddltilde_dthetadthetaprime_1(dimdelta+k,dimdelta+l) = ddltilde_dthetadthetaprime_1(dimdelta+l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_dbeta_1_l, dHtilde_dbeta_1_k, ddHtilde_dbetadbetaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\beta_\ell
    }
  }
  for(int k = 0; k < dimbeta; k++) {
    for(int l = 0; l < dimdelta; l++) {
      arma::mat dDtilde_ddelta_1_l = dDtilde_ddelta_1_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
      arma::mat dRtilde_dbeta_1_k = dRtilde_dbeta_1_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_k
      arma::mat dHtilde_ddelta_1_l = dHtilde_ddelta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dHtilde_dbeta_1_k = dHtilde_dbeta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
      arma::mat ddRtilde_dbetaddeltaprime_1_kl = ddRtilde_dbetaddeltaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
      arma::mat ddHtilde_dbetaddeltaprime_1_kl = fc_ddHtilde_dbetaddeltaprime_t_kl(Dtilde_1, dDtilde_ddelta_1_l, dRtilde_dbeta_1_k, ddRtilde_dbetaddeltaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
      ddltilde_dthetadthetaprime_1(dimdelta+k,l) = ddltilde_dthetadthetaprime_1(l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_ddelta_1_l, dHtilde_dbeta_1_k, ddHtilde_dbetaddeltaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\delta_\ell
    }
  }
  // \widehat{\Sigma}_{*, 1}
  arma::mat SigmaStarhat_1 = ddltilde_dthetadthetaprime_1; // \widehat{\Sigma}_{*, 1}
  nSigmaStarhat = nSigmaStarhat + SigmaStarhat_1;
  // t=2,...,n
  // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
  arma::mat varepsilontilde_tminusBkminus1TOtminus2_mat(m, Bk); varepsilontilde_tminusBkminus1TOtminus2_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta))
  varepsilontilde_tminusBkminus1TOtminus2_mat.each_col() = exp(-0.5*omegauline); // for s <= 0, \widetilde{\varepsilon}_s(\bm\delta) = (\exp{-1/2*\underliner{\omega}_1}, ..., \exp{-1/2*\underliner{\omega}_m})'
  arma::vec varepsilontilde_tminus1 = varepsilontilde_1; // \widetilde{\varepsilon}_{t-1}(\bm\delta)
  // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
  arma::cube dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.fill(0.0); // the i-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-2+i}(\bm\delta) / \partial\bm\delta'_\ell
  dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.each_slice() = dlnhulinetilde_ddelta_1_mat;
  arma::mat dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_1_mat;
  // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_k\partial\delta_\ell, and \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field; // the (k, \ell, i)-th element is the vector \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-2+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube = ddlnhulinetilde_ddeltaddeltaprime_1_cube;
  // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell, \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell, and \partial^2\widetilde{R}_{t-1}(\bm\theta)/\partial\theta_\k\partial\theta_\ell
  arma::mat Rtilde_tminus1 = Rtilde_1; // \widetilde{R}_{t-1}(\bm\theta)
  arma::cube dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
  arma::cube dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_\ell
  arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_tminus1_field = ddRtilde_ddeltaddeltaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddRtilde_dbetadbetaprime_tminus1_field = ddRtilde_dbetadbetaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::field<arma::mat> ddRtilde_dbetaddeltaprime_tminus1_field = ddRtilde_dbetaddeltaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\delta_\ell
  for(int t = 2; t < (n+1); t++) {
    arma::vec y_t = y_m_n.col(t-1); // \bm{y}_{t}
    arma::vec lnhulinetilde_t = fc_general_lnhulinetilde_t_general(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s); // \ln\widetilde\underline{\bm{h}}_t(\bm\delta)
    arma::mat Dtilde_t = fc_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}(\bm\delta)
    arma::mat inverse_Dtilde_t = fc_inverse_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}^{-1}(\bm\delta)
    arma::vec varepsilontilde_t = inverse_Dtilde_t * y_t; // \widetilde{\varepsilon}_{t}(\bm\delta) = \widetilde{D}_{t}^{-1}(\bm\delta) \bm{y}_{t}
    arma::mat varepsilontilde_tminusBkTOtminus1_mat(m, Bk); varepsilontilde_tminusBkTOtminus1_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk}(\bm\delta), ..., \widetilde{\varepsilon}_{t-1}(\bm\delta))
    varepsilontilde_tminusBkTOtminus1_mat.cols(0, Bk-2) = varepsilontilde_tminusBkminus1TOtminus2_mat.cols(1, Bk-1);
    varepsilontilde_tminusBkTOtminus1_mat.col(Bk-1) = varepsilontilde_tminus1;
    arma::mat Psitilde_tminus1 = fc_Psi(m, varepsilontilde_tminusBkTOtminus1_mat); // \widetilde{\Psi}_{t-1}(\bm\delta)
    arma::mat Rtilde_t = fc_Rtilde_t(m, Bk, beta_1, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1); // \widetilde{R}_{t}(\bm\theta)
    arma::mat Htilde_t = fc_Htilde_t(Dtilde_t, Rtilde_t); // \widetilde{H}_{t}(\bm\theta)
    // first derivatives
    arma::mat dlnhulinetilde_ddelta_t_mat = fc_general_dlnhulinetilde_ddelta_t_mat_general(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_lambdaklnyuline_1_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'
    arma::cube dlnhulinetilde_ddelta_tminusBkTOtminus1_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.fill(0.0); // the k-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-k)}(\bm\delta) / \partial\bm\delta'
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slices(0, Bk-2) = dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.slices(1, Bk-1);
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slice(Bk-1) = dlnhulinetilde_ddelta_tminus1_mat;
    arma::cube dDtilde_ddelta_t_cube = fc_dDtilde_ddelta_t_cube(m, dimdelta, Dtilde_t, dlnhulinetilde_ddelta_t_mat); // the \ell-th slice is \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
    arma::cube dPsitilde_ddelta_tminus1_cube = fc_dPsitilde_ddelta_tminus1_cube(m, dimdelta, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dlnhulinetilde_ddelta_tminusBkTOtminus1_cube); // the \ell-th slice is \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
    arma::cube dRtilde_ddelta_t_cube = fc2_dRtilde_ddelta_t_cube(m, dimdelta, beta_1, beta_2, dPsitilde_ddelta_tminus1_cube, dRtilde_ddelta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
    arma::cube dRtilde_dbeta_t_cube = fc_dRtilde_dbeta_t_cube(m, dimbeta, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1, dRtilde_dbeta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\beta_\ell
    arma::cube dHtilde_ddelta_t_cube = fc_dHtilde_ddelta_t_cube(m, dimdelta, Dtilde_t, Rtilde_t, dDtilde_ddelta_t_cube, dRtilde_ddelta_t_cube); // the \ell-th slice is \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
    arma::cube dHtilde_dbeta_t_cube = fc_dHtilde_dbeta_t_cube(m, dimbeta, Dtilde_t, dRtilde_dbeta_t_cube); // the \ell-th slice is \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
    // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
    arma::vec dltilde_dtheta_t = fc2_dltilde_dtheta_t(m, dimdelta, dimbeta, dimtheta, y_t, Htilde_t, dHtilde_ddelta_t_cube, dHtilde_dbeta_t_cube); // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
    // \widehat{\Sigma}_t
    arma::mat dltilde_dtheta_t_mat = fc_asmat(dltilde_dtheta_t, dimtheta, 1);
    arma::mat Sigmahat_t = dltilde_dtheta_t * dltilde_dtheta_t_mat.t(); // \widehat{\Sigma}_t
    nSigmahat = nSigmahat + Sigmahat_t;
    // second derivatives
    arma::cube ddlnhulinetilde_ddeltaddeltaprime_t_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube_general(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_1_m_nminus1_r, sum_lambdaklnyuline_2_m_nminus1_r, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s, sum_gammaklnyuline_31_m_nminus1_s, sum_gammaklnyuline_32_m_nminus1_s, sum_gammaklnyuline_41_m_nminus1_s, sum_gammaklnyuline_42_m_nminus1_s, sum_gammaklnyuline_51_m_nminus1_s, sum_gammaklnyuline_52_m_nminus1_s); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_t_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddRtilde_dbetadbetaprime_t_field(dimbeta, dimbeta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
    arma::field<arma::mat> ddRtilde_dbetaddeltaprime_t_field(dimbeta, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
    // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix, of which the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta \partial\bm\theta'
    arma::mat ddltilde_dthetadthetaprime_t(dimtheta, dimtheta); ddltilde_dthetadthetaprime_t.fill(0.0);
    for(int k = 0; k < dimdelta; k++) {
      for(int l = 0; l < (k+1); l++) {
        arma::vec dlnhulinetilde_ddelta_t_l = dlnhulinetilde_ddelta_t_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
        arma::vec dlnhulinetilde_ddelta_t_k = dlnhulinetilde_ddelta_t_mat.col(k); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k
        arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.col(l); // the i-th column is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_\ell
        arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.col(k); // the i-th column is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_k
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_t_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\bm\delta'
        arma::vec ddlnhulinetilde_ddeltaddeltaprime_t_kl = ddlnhulinetilde_ddeltaddeltaprime_t_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field(k, l); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\bm\delta'
        arma::vec ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl = ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat(m, Bk); ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.fill(0.0); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.cols(0, Bk-2) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_kl_mat.cols(1, Bk-1);
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.col(Bk-1) = ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl;
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(k, l) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(l, k) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat;
        arma::mat dDtilde_ddelta_t_l = dDtilde_ddelta_t_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
        arma::mat dDtilde_ddelta_t_k = dDtilde_ddelta_t_cube.slice(k); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_k
        arma::mat dPsitilde_ddelta_tminus1_k = dPsitilde_ddelta_tminus1_cube.slice(k); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_k
        arma::mat dRtilde_ddelta_t_l = dRtilde_ddelta_t_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dRtilde_ddelta_t_k = dRtilde_ddelta_t_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_k
        arma::mat dHtilde_ddelta_t_l = dHtilde_ddelta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_ddelta_t_k = dHtilde_ddelta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_k
        arma::mat ddDtilde_ddeltaddeltaprime_t_kl = fc_ddDtilde_ddeltaddeltaprime_t_kl(dlnhulinetilde_ddelta_t_l, dlnhulinetilde_ddelta_t_k, ddlnhulinetilde_ddeltaddeltaprime_t_kl, Dtilde_t); // \partial^2\widetilde{D}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddRtilde_ddeltaddeltaprime_tminus1_kl = ddRtilde_ddeltaddeltaprime_tminus1_field(k,l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddRtilde_ddeltaddeltaprime_t_kl = fc_ddRtilde_ddeltaddeltaprime_t_kl(m, beta_1, beta_2, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dPsitilde_ddelta_tminus1_k, dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat, dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat, ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat, ddRtilde_ddeltaddeltaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
        ddRtilde_ddeltaddeltaprime_t_field(k, l) = ddRtilde_ddeltaddeltaprime_t_field(l, k) = ddRtilde_ddeltaddeltaprime_t_kl;
        arma::mat ddHtilde_ddeltaddeltaprime_t_kl = fc_ddHtilde_ddeltaddeltaprime_t_kl(Dtilde_t, Rtilde_t, dDtilde_ddelta_t_l, dDtilde_ddelta_t_k, dRtilde_ddelta_t_l, dRtilde_ddelta_t_k, ddDtilde_ddeltaddeltaprime_t_kl, ddRtilde_ddeltaddeltaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
        ddltilde_dthetadthetaprime_t(k,l) = ddltilde_dthetadthetaprime_t(l,k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_ddelta_t_l, dHtilde_ddelta_t_k, ddHtilde_ddeltaddeltaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\delta_k\partial\delta_\ell
      }
    }
    for(int k = 0; k < dimbeta; k++) {
      for(int l = 0; l < (k+1); l++) {
        arma::mat dRtilde_dbeta_tminus1_l = dRtilde_dbeta_tminus1_cube.slice(l); // \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_\ell
        arma::mat dHtilde_dbeta_t_l = dHtilde_dbeta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
        arma::mat dHtilde_dbeta_t_k = dHtilde_dbeta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
        arma::mat ddRtilde_dbetadbetaprime_tminus1_kl = ddRtilde_dbetadbetaprime_tminus1_field(k, l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\beta_\ell
        arma::mat ddRtilde_dbetadbetaprime_t_kl = fc_ddRtilde_dbetadbetaprime_t_kl(l, k, m, beta_2, dRtilde_dbeta_tminus1_l, ddRtilde_dbetadbetaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
        ddRtilde_dbetadbetaprime_t_field(k, l) = ddRtilde_dbetadbetaprime_t_field(l, k) = ddRtilde_dbetadbetaprime_t_kl;
        arma::mat ddHtilde_dbetadbetaprime_t_kl = fc_ddHtilde_dbetadbetaprime_t_kl(Dtilde_t, ddRtilde_dbetadbetaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
        ddltilde_dthetadthetaprime_t(dimdelta+k,dimdelta+l) = ddltilde_dthetadthetaprime_t(dimdelta+l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_dbeta_t_l, dHtilde_dbeta_t_k, ddHtilde_dbetadbetaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\beta_\ell
      }
    }
    for(int k = 0; k < dimbeta; k++) {
      for(int l = 0; l < dimdelta; l++) {
        arma::mat dDtilde_ddelta_t_l = dDtilde_ddelta_t_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
        arma::mat dPsitilde_ddelta_tminus1_l = dPsitilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
        arma::mat dRtilde_dbeta_t_k = dRtilde_dbeta_t_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_k
        arma::mat dRtilde_ddelta_tminus1_l = dRtilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_ddelta_t_l = dHtilde_ddelta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_dbeta_t_k = dHtilde_dbeta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
        arma::mat ddRtilde_dbetaddeltaprime_tminus1_kl = ddRtilde_dbetaddeltaprime_tminus1_field(k, l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\delta_\ell
        arma::mat ddRtilde_dbetaddeltaprime_t_kl = fc_ddRtilde_dbetaddeltaprime_t_kl(l, k, m, beta_2, dPsitilde_ddelta_tminus1_l, dRtilde_ddelta_tminus1_l, ddRtilde_dbetaddeltaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
        ddRtilde_dbetaddeltaprime_t_field(k, l) = ddRtilde_dbetaddeltaprime_t_kl;
        arma::mat ddHtilde_dbetaddeltaprime_t_kl = fc_ddHtilde_dbetaddeltaprime_t_kl(Dtilde_t, dDtilde_ddelta_t_l, dRtilde_dbeta_t_k, ddRtilde_dbetaddeltaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
        ddltilde_dthetadthetaprime_t(dimdelta+k,l) = ddltilde_dthetadthetaprime_t(l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_ddelta_t_l, dHtilde_dbeta_t_k, ddHtilde_dbetaddeltaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\delta_\ell
      }
    }
    // \widehat{\Sigma}_{*, t}
    arma::mat SigmaStarhat_t = ddltilde_dthetadthetaprime_t; // \widehat{\Sigma}_{*, t}
    nSigmaStarhat = nSigmaStarhat + SigmaStarhat_t;
    // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
    varepsilontilde_tminusBkminus1TOtminus2_mat = varepsilontilde_tminusBkTOtminus1_mat;
    varepsilontilde_tminus1 = varepsilontilde_t;
    // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
    dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube;
    dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_t_mat;
    // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_k\partial\delta_\ell and \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field;
    ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube = ddlnhulinetilde_ddeltaddeltaprime_t_cube;
    // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell, \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell, and \partial^2\widetilde{R}_{t-1}(\bm\theta)/\partial\theta_\k\partial\theta_\ell
    Rtilde_tminus1 = Rtilde_t;
    dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_t_cube;
    dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_t_cube;
    ddRtilde_ddeltaddeltaprime_tminus1_field = ddRtilde_ddeltaddeltaprime_t_field;
    ddRtilde_dbetadbetaprime_tminus1_field = ddRtilde_dbetadbetaprime_t_field;
    ddRtilde_dbetaddeltaprime_tminus1_field = ddRtilde_dbetaddeltaprime_t_field;
  }
  arma::mat Sigmahat = nSigmahat / (n * 1.0); // \widehat{\Sigma} = 1/n \sum_{t=1}^{n} \widehat{\Sigma}_t
  arma::mat SigmaStarhat = nSigmaStarhat / (n * 1.0); // \widehat{\Sigma}_{*} = 1/n \sum_{t=1}^{n} \widehat{\Sigma}_{*, t}
  arma::mat Sigmahat_transf = Deltahat_theta_vartheta.t() * Sigmahat * Deltahat_theta_vartheta; // \Delta_\theta' \widehat{\Sigma} \Delta_\theta
  arma::mat SigmaStarhat_transf = Deltahat_theta_vartheta.t() * SigmaStarhat * Deltahat_theta_vartheta; // \Delta_\theta' \widehat{\Sigma}_{*} \Delta_\theta
  arma::mat CovMat_thetahat = Deltahat_theta_vartheta * pinv(SigmaStarhat_transf) * Sigmahat_transf * pinv(SigmaStarhat_transf) * Deltahat_theta_vartheta.t(); // \Delta_\theta (\Delta_\theta' \widehat{\Sigma}_{*} \Delta_\theta)^{-1} (\Delta_\theta' \widehat{\Sigma} \Delta_\theta) (\Delta_\theta' \widehat{\Sigma}_{*} \Delta_\theta)^{-1} \Delta_\theta'
  arma::mat transfmat = fc_transfmat_theta_to_vecPhi1(m, r, s, dimtheta); // transformation matrix for trasfering \theta to vec(\Phi_1)
  arma::mat CovMat_vecPhi1 = transfmat * CovMat_thetahat * transfmat.t(); // covariance matrix of vec(\Phi_1)
  arma::vec var_sqrtnvecPhi1 = CovMat_vecPhi1.diag(0); // variances of \sqrt{n} vec(\Phi_1)
  arma::vec var_vecPhi1 = var_sqrtnvecPhi1 / (n * 1.0); // variances of vec(\Phi_1)
  arma::vec ASD_vecPhi1 = sqrt(var_vecPhi1); // ASD of vec(\Phi_1)
  return ASD_vecPhi1;
}

// [[Rcpp::export]]
arma::vec fc_general_ASD_vecPhi1hat(int n, int m, int r, int s, int Bk, arma::vec thetahat, arma::mat y_m_n) {
  // ASD of \Phi_1 = \sum_{k=1}^{r} G_{0,k} + \sum_{k=1}^{s} G_{1,k}
  int dimkappa = r + 2*s + r*m*m + 2*s*m*m;
  int dimbeta = 2 + m*(m-1)/2;
  int dimdelta = m + dimkappa;
  int dimtheta = dimdelta + dimbeta;
  arma::vec omegauline = thetahat.subvec(0, m-1); // \underline{\bm\omega}
  arma::vec kappa_vec = thetahat.subvec(m, m+dimkappa-1); // \bm\kappa
  arma::vec beta_vec = thetahat.subvec(m+dimkappa, m+dimkappa+dimbeta-1); // \bm\beta
  double beta_1 = beta_vec(0); double beta_2 = beta_vec(1); // \beta_1 and \beta_2
  arma::vec ruline = beta_vec.tail(dimbeta-2); // \underline\bm{r}
  arma::mat Ruline = fc_Ruline(m, ruline); // \underline{R}
  // summations using FFT algorithm
  int max_r_1 = r; if(max_r_1 == 0) {max_r_1 = 1;}
  int max_s_1 = s; if(max_s_1 == 0) {max_s_1 = 1;}
  arma::cube sum_lambdaklnyuline_0_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_0_m_nminus1_r.fill(0.0);
  arma::cube sum_lambdaklnyuline_1_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_1_m_nminus1_r.fill(0.0);
  arma::cube sum_lambdaklnyuline_2_m_nminus1_r(m, n-1, max_r_1); sum_lambdaklnyuline_2_m_nminus1_r.fill(0.0);
  arma::cube sum_gammaklnyuline_01_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_01_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_02_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_02_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_11_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_11_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_12_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_12_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_21_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_21_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_22_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_22_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_31_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_31_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_32_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_32_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_41_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_41_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_42_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_42_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_51_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_51_m_nminus1_s.fill(0.0);
  arma::cube sum_gammaklnyuline_52_m_nminus1_s(m, n-1, max_s_1); sum_gammaklnyuline_52_m_nminus1_s.fill(0.0);
  if(s == 0) {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_2_m_nminus1_r = fc_sum_lambdaklnyuline_2_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
  } else if(r == 0) {
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_31_m_nminus1_s = fc_sum_gammaklnyuline_31_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_32_m_nminus1_s = fc_sum_gammaklnyuline_32_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_41_m_nminus1_s = fc_sum_gammaklnyuline_41_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_42_m_nminus1_s = fc_sum_gammaklnyuline_42_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_51_m_nminus1_s = fc_sum_gammaklnyuline_51_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_52_m_nminus1_s = fc_sum_gammaklnyuline_52_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  } else {
    sum_lambdaklnyuline_0_m_nminus1_r = fc_sum_lambdaklnyuline_0_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_1_m_nminus1_r = fc_sum_lambdaklnyuline_1_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_lambdaklnyuline_2_m_nminus1_r = fc_sum_lambdaklnyuline_2_m_nminus1_r(n, m, r, kappa_vec, y_m_n);
    sum_gammaklnyuline_01_m_nminus1_s = fc_sum_gammaklnyuline_01_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_02_m_nminus1_s = fc_sum_gammaklnyuline_02_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_11_m_nminus1_s = fc_sum_gammaklnyuline_11_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_12_m_nminus1_s = fc_sum_gammaklnyuline_12_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_21_m_nminus1_s = fc_sum_gammaklnyuline_21_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_22_m_nminus1_s = fc_sum_gammaklnyuline_22_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_31_m_nminus1_s = fc_sum_gammaklnyuline_31_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_32_m_nminus1_s = fc_sum_gammaklnyuline_32_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_41_m_nminus1_s = fc_sum_gammaklnyuline_41_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_42_m_nminus1_s = fc_sum_gammaklnyuline_42_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_51_m_nminus1_s = fc_sum_gammaklnyuline_51_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
    sum_gammaklnyuline_52_m_nminus1_s = fc_sum_gammaklnyuline_52_m_nminus1_s(n, m, r, s, kappa_vec, y_m_n);
  }
  // n \widehat{\Sigma}
  arma::mat nSigmahat(dimtheta, dimtheta); nSigmahat.fill(0.0); // n \widehat{\Sigma} = \sum_{t=1}^{n} \widehat{\Sigma}_t = \sum_{t=1}^{n} \partial\widetilde{\ell}_t(\widehat\bm\theta)/\partial\bm\theta \partial\widetilde{\ell}_t(\widehat\bm\theta)/\partial\bm\theta'
  // n \widehat{\Sigma}_{*}
  arma::mat nSigmaStarhat(dimtheta, dimtheta); nSigmaStarhat.fill(0.0); // n \widehat{\Sigma}_{*} = \sum_{t=1}^{n} \widehat{\Sigma}_{*, t} = \sum_{t=1}^{n} \partial^2\widetilde{\ell}_t(\widehat\bm\theta) / \partial\bm\theta\partial\bm\theta'
  // t=1
  arma::vec y_1 = y_m_n.col(0); // \bm{y}_{1}
  arma::mat Dtilde_1 = fc_Dtilde_1(omegauline); // \widetilde{D}_{1}(\bm\delta)
  arma::mat inverse_Dtilde_1 = fc_inverse_Dtilde_1(omegauline); // \widetilde{D}_{1}^{-1}(\bm\delta)
  arma::mat Rtilde_1 = fc_Rtilde_1(m, beta_1, beta_2, Ruline); // \widetilde{R}_{1}(\bm\theta)
  arma::mat Htilde_1 = fc_Htilde_t(Dtilde_1, Rtilde_1); // \widetilde{H}_{1}(\bm\theta)
  arma::vec varepsilontilde_1 = inverse_Dtilde_1 * y_1; // \widetilde{\varepsilon}_{1}(\bm\delta) = \widetilde{D}_{1}^{-1}(\bm\delta) \bm{y}_{1}
  // first derivatives
  arma::mat dlnhulinetilde_ddelta_1_mat = fc_general_dlnhulinetilde_ddelta_1_mat(m, dimdelta); // \partial\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\bm\delta'_\ell, which is equal to \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'_\ell for t <= 0
  arma::cube dDtilde_ddelta_1_cube = fc_dDtilde_ddelta_t_cube(m, dimdelta, Dtilde_1, dlnhulinetilde_ddelta_1_mat); // the \ell-th slice is \partial\widetilde{D}_1(\bm\delta) / \partial\delta_\ell
  arma::cube dRtilde_ddelta_1_cube = fc_dRtilde_ddelta_1_cube(m, dimdelta); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\delta_\ell
  arma::cube dRuline_dbeta_cube = fc_dRuline_dbeta_cube(m, dimbeta); // the \ell-th slice is \partial\underline{R} / \partial\beta_\ell
  arma::cube dRtilde_dbeta_1_cube = fc_dRtilde_dbeta_1_cube(m, dimbeta, beta_1, beta_2, Ruline, dRuline_dbeta_cube); // the \ell-th slice is \partial\widetilde{R}_1(\bm\theta) / \partial\beta_\ell
  arma::cube dHtilde_ddelta_1_cube = fc_dHtilde_ddelta_t_cube(m, dimdelta, Dtilde_1, Rtilde_1, dDtilde_ddelta_1_cube, dRtilde_ddelta_1_cube); // the \ell-th slice is \partial\widetilde{H}_1(\bm\theta) / \partial\delta_\ell
  arma::cube dHtilde_dbeta_1_cube = fc_dHtilde_dbeta_t_cube(m, dimbeta, Dtilde_1, dRtilde_dbeta_1_cube); // the \ell-th slice is \partial\widetilde{H}_1(\bm\theta) / \partial\beta_\ell
  // \widehat{\Sigma}_1
  arma::vec dltilde_dtheta_1 = fc2_dltilde_dtheta_t(m, dimdelta, dimbeta, dimtheta, y_1, Htilde_1, dHtilde_ddelta_1_cube, dHtilde_dbeta_1_cube); // \partial\widetilde{\ell}_1(\bm\theta) / \partial\bm\theta
  arma::mat dltilde_dtheta_1_mat = fc_asmat(dltilde_dtheta_1, dimtheta, 1);
  arma::mat Sigmahat_1 = dltilde_dtheta_1 * dltilde_dtheta_1_mat.t(); // \widehat{\Sigma}_1
  nSigmahat = nSigmahat + Sigmahat_1;
  // second derivatives
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_1_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_1_cube(m, dimdelta); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_1(\bm\delta) / \partial\delta_k\partial\delta_\ell, which is equal to \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell for t <= 0
  arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_1_field = fc_ddRtilde_ddeltaddeltaprime_1_field(m, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddRtilde_dbetadbetaprime_1_field = fc_ddRtilde_dbetadbetaprime_1_field(m, dimbeta, beta_1, beta_2, Ruline, dRuline_dbeta_cube); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::field<arma::mat> ddRtilde_dbetaddeltaprime_1_field = fc_ddRtilde_dbetaddeltaprime_1_field(m, dimdelta, dimbeta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_1(\bm\theta) / \partial\beta_k\partial\delta_\ell
  // \partial^2\ln\widetilde\underline{\bm{h}}_{1-\Bbbk}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{0}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix, of which the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{-\Bbbk+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  // \partial^2\widetilde{\ell}_{1}(\bm\theta) / \partial\bm\theta \partial\bm\theta'
  arma::mat ddltilde_dthetadthetaprime_1(dimtheta, dimtheta); ddltilde_dthetadthetaprime_1.fill(0.0);
  for(int k = 0; k < dimdelta; k++) {
    for(int l = 0; l < (k+1); l++) {
      arma::vec dlnhulinetilde_ddelta_1_l = dlnhulinetilde_ddelta_1_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
      arma::vec dlnhulinetilde_ddelta_1_k = dlnhulinetilde_ddelta_1_mat.col(k); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k
      arma::mat ddlnhulinetilde_ddeltaddeltaprime_1_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_1_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\bm\delta'
      arma::vec ddlnhulinetilde_ddeltaddeltaprime_1_kl = ddlnhulinetilde_ddeltaddeltaprime_1_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat(m, Bk); ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat.fill(0.0); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{-\Bbbk+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
      ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat.each_col() = ddlnhulinetilde_ddeltaddeltaprime_1_kl;
      ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(k, l) = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field(l, k) = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_kl_mat;
      arma::mat dDtilde_ddelta_1_l = dDtilde_ddelta_1_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
      arma::mat dDtilde_ddelta_1_k = dDtilde_ddelta_1_cube.slice(k); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_k
      arma::mat dRtilde_ddelta_1_l = dRtilde_ddelta_1_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dRtilde_ddelta_1_k = dRtilde_ddelta_1_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_k
      arma::mat dHtilde_ddelta_1_l = dHtilde_ddelta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dHtilde_ddelta_1_k = dHtilde_ddelta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_k
      arma::mat ddDtilde_ddeltaddeltaprime_1_kl = fc_ddDtilde_ddeltaddeltaprime_t_kl(dlnhulinetilde_ddelta_1_l, dlnhulinetilde_ddelta_1_k, ddlnhulinetilde_ddeltaddeltaprime_1_kl, Dtilde_1); // \partial^2\widetilde{D}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddRtilde_ddeltaddeltaprime_1_kl = ddRtilde_ddeltaddeltaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
      arma::mat ddHtilde_ddeltaddeltaprime_1_kl = fc_ddHtilde_ddeltaddeltaprime_t_kl(Dtilde_1, Rtilde_1, dDtilde_ddelta_1_l, dDtilde_ddelta_1_k, dRtilde_ddelta_1_l, dRtilde_ddelta_1_k, ddDtilde_ddeltaddeltaprime_1_kl, ddRtilde_ddeltaddeltaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
      ddltilde_dthetadthetaprime_1(k,l) = ddltilde_dthetadthetaprime_1(l,k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_ddelta_1_l, dHtilde_ddelta_1_k, ddHtilde_ddeltaddeltaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\delta_k\partial\delta_\ell
    }
  }
  for(int k = 0; k < dimbeta; k++) {
    for(int l = 0; l < (k+1); l++) {
      arma::mat dHtilde_dbeta_1_l = dHtilde_dbeta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
      arma::mat dHtilde_dbeta_1_k = dHtilde_dbeta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
      arma::mat ddRtilde_dbetadbetaprime_1_kl = ddRtilde_dbetadbetaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
      arma::mat ddHtilde_dbetadbetaprime_1_kl = fc_ddHtilde_dbetadbetaprime_t_kl(Dtilde_1, ddRtilde_dbetadbetaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
      ddltilde_dthetadthetaprime_1(dimdelta+k,dimdelta+l) = ddltilde_dthetadthetaprime_1(dimdelta+l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_dbeta_1_l, dHtilde_dbeta_1_k, ddHtilde_dbetadbetaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\beta_\ell
    }
  }
  for(int k = 0; k < dimbeta; k++) {
    for(int l = 0; l < dimdelta; l++) {
      arma::mat dDtilde_ddelta_1_l = dDtilde_ddelta_1_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
      arma::mat dRtilde_dbeta_1_k = dRtilde_dbeta_1_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_k
      arma::mat dHtilde_ddelta_1_l = dHtilde_ddelta_1_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
      arma::mat dHtilde_dbeta_1_k = dHtilde_dbeta_1_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
      arma::mat ddRtilde_dbetaddeltaprime_1_kl = ddRtilde_dbetaddeltaprime_1_field(k,l); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
      arma::mat ddHtilde_dbetaddeltaprime_1_kl = fc_ddHtilde_dbetaddeltaprime_t_kl(Dtilde_1, dDtilde_ddelta_1_l, dRtilde_dbeta_1_k, ddRtilde_dbetaddeltaprime_1_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
      ddltilde_dthetadthetaprime_1(dimdelta+k,l) = ddltilde_dthetadthetaprime_1(l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_1, Htilde_1, dHtilde_ddelta_1_l, dHtilde_dbeta_1_k, ddHtilde_dbetaddeltaprime_1_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\delta_\ell
    }
  }
  // \widehat{\Sigma}_{*, 1}
  arma::mat SigmaStarhat_1 = ddltilde_dthetadthetaprime_1; // \widehat{\Sigma}_{*, 1}
  nSigmaStarhat = nSigmaStarhat + SigmaStarhat_1;
  // t=2,...,n
  // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
  arma::mat varepsilontilde_tminusBkminus1TOtminus2_mat(m, Bk); varepsilontilde_tminusBkminus1TOtminus2_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta))
  varepsilontilde_tminusBkminus1TOtminus2_mat.each_col() = exp(-0.5*omegauline); // for s <= 0, \widetilde{\varepsilon}_s(\bm\delta) = (\exp{-1/2*\underliner{\omega}_1}, ..., \exp{-1/2*\underliner{\omega}_m})'
  arma::vec varepsilontilde_tminus1 = varepsilontilde_1; // \widetilde{\varepsilon}_{t-1}(\bm\delta)
  // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
  arma::cube dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.fill(0.0); // the i-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-2+i}(\bm\delta) / \partial\bm\delta'_\ell
  dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.each_slice() = dlnhulinetilde_ddelta_1_mat;
  arma::mat dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_1_mat;
  // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_k\partial\delta_\ell, and \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field = ddlnhulinetilde_ddeltaddeltaprime_1minusBkTO0_field; // the (k, \ell, i)-th element is the vector \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-2+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
  arma::cube ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube = ddlnhulinetilde_ddeltaddeltaprime_1_cube;
  // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell, \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell, and \partial^2\widetilde{R}_{t-1}(\bm\theta)/\partial\theta_\k\partial\theta_\ell
  arma::mat Rtilde_tminus1 = Rtilde_1; // \widetilde{R}_{t-1}(\bm\theta)
  arma::cube dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
  arma::cube dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_1_cube; // the \ell-th slice is \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_\ell
  arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_tminus1_field = ddRtilde_ddeltaddeltaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_k\partial\delta_\ell
  arma::field<arma::mat> ddRtilde_dbetadbetaprime_tminus1_field = ddRtilde_dbetadbetaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\beta_\ell
  arma::field<arma::mat> ddRtilde_dbetaddeltaprime_tminus1_field = ddRtilde_dbetaddeltaprime_1_field; // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\delta_\ell
  for(int t = 2; t < (n+1); t++) {
    arma::vec y_t = y_m_n.col(t-1); // \bm{y}_{t}
    arma::vec lnhulinetilde_t = fc_general_lnhulinetilde_t_general(t, m, r, s, omegauline, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s); // \ln\widetilde\underline{\bm{h}}_t(\bm\delta)
    arma::mat Dtilde_t = fc_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}(\bm\delta)
    arma::mat inverse_Dtilde_t = fc_inverse_Dtilde_t(lnhulinetilde_t); // \widetilde{D}_{t}^{-1}(\bm\delta)
    arma::vec varepsilontilde_t = inverse_Dtilde_t * y_t; // \widetilde{\varepsilon}_{t}(\bm\delta) = \widetilde{D}_{t}^{-1}(\bm\delta) \bm{y}_{t}
    arma::mat varepsilontilde_tminusBkTOtminus1_mat(m, Bk); varepsilontilde_tminusBkTOtminus1_mat.fill(0.0); // (\widetilde{\varepsilon}_{t-\Bbbk}(\bm\delta), ..., \widetilde{\varepsilon}_{t-1}(\bm\delta))
    varepsilontilde_tminusBkTOtminus1_mat.cols(0, Bk-2) = varepsilontilde_tminusBkminus1TOtminus2_mat.cols(1, Bk-1);
    varepsilontilde_tminusBkTOtminus1_mat.col(Bk-1) = varepsilontilde_tminus1;
    arma::mat Psitilde_tminus1 = fc_Psi(m, varepsilontilde_tminusBkTOtminus1_mat); // \widetilde{\Psi}_{t-1}(\bm\delta)
    arma::mat Rtilde_t = fc_Rtilde_t(m, Bk, beta_1, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1); // \widetilde{R}_{t}(\bm\theta)
    arma::mat Htilde_t = fc_Htilde_t(Dtilde_t, Rtilde_t); // \widetilde{H}_{t}(\bm\theta)
    // first derivatives
    arma::mat dlnhulinetilde_ddelta_t_mat = fc_general_dlnhulinetilde_ddelta_t_mat_general(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_0_m_nminus1_r, sum_lambdaklnyuline_1_m_nminus1_r, sum_gammaklnyuline_01_m_nminus1_s, sum_gammaklnyuline_02_m_nminus1_s, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\bm\delta'
    arma::cube dlnhulinetilde_ddelta_tminusBkTOtminus1_cube(m, dimdelta, Bk); dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.fill(0.0); // the k-th slice is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-k)}(\bm\delta) / \partial\bm\delta'
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slices(0, Bk-2) = dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube.slices(1, Bk-1);
    dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.slice(Bk-1) = dlnhulinetilde_ddelta_tminus1_mat;
    arma::cube dDtilde_ddelta_t_cube = fc_dDtilde_ddelta_t_cube(m, dimdelta, Dtilde_t, dlnhulinetilde_ddelta_t_mat); // the \ell-th slice is \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
    arma::cube dPsitilde_ddelta_tminus1_cube = fc_dPsitilde_ddelta_tminus1_cube(m, dimdelta, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dlnhulinetilde_ddelta_tminusBkTOtminus1_cube); // the \ell-th slice is \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
    arma::cube dRtilde_ddelta_t_cube = fc2_dRtilde_ddelta_t_cube(m, dimdelta, beta_1, beta_2, dPsitilde_ddelta_tminus1_cube, dRtilde_ddelta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
    arma::cube dRtilde_dbeta_t_cube = fc_dRtilde_dbeta_t_cube(m, dimbeta, beta_2, Ruline, Psitilde_tminus1, Rtilde_tminus1, dRtilde_dbeta_tminus1_cube); // the \ell-th slice is \partial\widetilde{R}_t(\bm\theta) / \partial\beta_\ell
    arma::cube dHtilde_ddelta_t_cube = fc_dHtilde_ddelta_t_cube(m, dimdelta, Dtilde_t, Rtilde_t, dDtilde_ddelta_t_cube, dRtilde_ddelta_t_cube); // the \ell-th slice is \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
    arma::cube dHtilde_dbeta_t_cube = fc_dHtilde_dbeta_t_cube(m, dimbeta, Dtilde_t, dRtilde_dbeta_t_cube); // the \ell-th slice is \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
    // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
    arma::vec dltilde_dtheta_t = fc2_dltilde_dtheta_t(m, dimdelta, dimbeta, dimtheta, y_t, Htilde_t, dHtilde_ddelta_t_cube, dHtilde_dbeta_t_cube); // \partial\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta
    // \widehat{\Sigma}_t
    arma::mat dltilde_dtheta_t_mat = fc_asmat(dltilde_dtheta_t, dimtheta, 1);
    arma::mat Sigmahat_t = dltilde_dtheta_t * dltilde_dtheta_t_mat.t(); // \widehat{\Sigma}_t
    nSigmahat = nSigmahat + Sigmahat_t;
    // second derivatives
    arma::cube ddlnhulinetilde_ddeltaddeltaprime_t_cube = fc_general_ddlnhulinetilde_ddeltaddeltaprime_t_cube_general(t, m, r, s, dimdelta, kappa_vec, sum_lambdaklnyuline_1_m_nminus1_r, sum_lambdaklnyuline_2_m_nminus1_r, sum_gammaklnyuline_11_m_nminus1_s, sum_gammaklnyuline_12_m_nminus1_s, sum_gammaklnyuline_21_m_nminus1_s, sum_gammaklnyuline_22_m_nminus1_s, sum_gammaklnyuline_31_m_nminus1_s, sum_gammaklnyuline_32_m_nminus1_s, sum_gammaklnyuline_41_m_nminus1_s, sum_gammaklnyuline_42_m_nminus1_s, sum_gammaklnyuline_51_m_nminus1_s, sum_gammaklnyuline_52_m_nminus1_s); // the \ell-th column in the k-th slice is \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddRtilde_ddeltaddeltaprime_t_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddRtilde_dbetadbetaprime_t_field(dimbeta, dimbeta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
    arma::field<arma::mat> ddRtilde_dbetaddeltaprime_t_field(dimbeta, dimdelta); // the (k,\ell)-th element is the matrix \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
    // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    arma::field<arma::mat> ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(dimdelta, dimdelta); // the (k,\ell)-th element is the matrix, of which the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1+i}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\bm\theta \partial\bm\theta'
    arma::mat ddltilde_dthetadthetaprime_t(dimtheta, dimtheta); ddltilde_dthetadthetaprime_t.fill(0.0);
    for(int k = 0; k < dimdelta; k++) {
      for(int l = 0; l < (k+1); l++) {
        arma::vec dlnhulinetilde_ddelta_t_l = dlnhulinetilde_ddelta_t_mat.col(l); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_\ell
        arma::vec dlnhulinetilde_ddelta_t_k = dlnhulinetilde_ddelta_t_mat.col(k); // \partial\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k
        arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.col(l); // the i-th column is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_\ell
        arma::mat dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube.col(k); // the i-th column is \partial\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_k
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_t_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_t_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\bm\delta'
        arma::vec ddlnhulinetilde_ddeltaddeltaprime_t_kl = ddlnhulinetilde_ddeltaddeltaprime_t_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field(k, l); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl_mat = ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube.slice(k); // \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\bm\delta'
        arma::vec ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl = ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl_mat.col(l); // \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat(m, Bk); ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.fill(0.0); // the i-th column is \partial^2\ln\widetilde\underline{\bm{h}}_{t-(\Bbbk+1-i)}(\bm\delta) / \partial\delta_k\partial\delta_\ell
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.cols(0, Bk-2) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_kl_mat.cols(1, Bk-1);
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat.col(Bk-1) = ddlnhulinetilde_ddeltaddeltaprime_tminus1_kl;
        ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(k, l) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field(l, k) = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat;
        arma::mat dDtilde_ddelta_t_l = dDtilde_ddelta_t_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
        arma::mat dDtilde_ddelta_t_k = dDtilde_ddelta_t_cube.slice(k); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_k
        arma::mat dPsitilde_ddelta_tminus1_k = dPsitilde_ddelta_tminus1_cube.slice(k); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_k
        arma::mat dRtilde_ddelta_t_l = dRtilde_ddelta_t_cube.slice(l); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dRtilde_ddelta_t_k = dRtilde_ddelta_t_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\delta_k
        arma::mat dHtilde_ddelta_t_l = dHtilde_ddelta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_ddelta_t_k = dHtilde_ddelta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_k
        arma::mat ddDtilde_ddeltaddeltaprime_t_kl = fc_ddDtilde_ddeltaddeltaprime_t_kl(dlnhulinetilde_ddelta_t_l, dlnhulinetilde_ddelta_t_k, ddlnhulinetilde_ddeltaddeltaprime_t_kl, Dtilde_t); // \partial^2\widetilde{D}_t(\bm\delta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddRtilde_ddeltaddeltaprime_tminus1_kl = ddRtilde_ddeltaddeltaprime_tminus1_field(k,l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_k\partial\delta_\ell
        arma::mat ddRtilde_ddeltaddeltaprime_t_kl = fc_ddRtilde_ddeltaddeltaprime_t_kl(m, beta_1, beta_2, varepsilontilde_tminusBkTOtminus1_mat, Psitilde_tminus1, dPsitilde_ddelta_tminus1_k, dlnhulinetilde_ddelta_tminusBkTOtminus1_l_mat, dlnhulinetilde_ddelta_tminusBkTOtminus1_k_mat, ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_kl_mat, ddRtilde_ddeltaddeltaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
        ddRtilde_ddeltaddeltaprime_t_field(k, l) = ddRtilde_ddeltaddeltaprime_t_field(l, k) = ddRtilde_ddeltaddeltaprime_t_kl;
        arma::mat ddHtilde_ddeltaddeltaprime_t_kl = fc_ddHtilde_ddeltaddeltaprime_t_kl(Dtilde_t, Rtilde_t, dDtilde_ddelta_t_l, dDtilde_ddelta_t_k, dRtilde_ddelta_t_l, dRtilde_ddelta_t_k, ddDtilde_ddeltaddeltaprime_t_kl, ddRtilde_ddeltaddeltaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\delta_k\partial\delta_\ell
        ddltilde_dthetadthetaprime_t(k,l) = ddltilde_dthetadthetaprime_t(l,k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_ddelta_t_l, dHtilde_ddelta_t_k, ddHtilde_ddeltaddeltaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\delta_k\partial\delta_\ell
      }
    }
    for(int k = 0; k < dimbeta; k++) {
      for(int l = 0; l < (k+1); l++) {
        arma::mat dRtilde_dbeta_tminus1_l = dRtilde_dbeta_tminus1_cube.slice(l); // \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_\ell
        arma::mat dHtilde_dbeta_t_l = dHtilde_dbeta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_\ell
        arma::mat dHtilde_dbeta_t_k = dHtilde_dbeta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
        arma::mat ddRtilde_dbetadbetaprime_tminus1_kl = ddRtilde_dbetadbetaprime_tminus1_field(k, l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\beta_\ell
        arma::mat ddRtilde_dbetadbetaprime_t_kl = fc_ddRtilde_dbetadbetaprime_t_kl(l, k, m, beta_2, dRtilde_dbeta_tminus1_l, ddRtilde_dbetadbetaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
        ddRtilde_dbetadbetaprime_t_field(k, l) = ddRtilde_dbetadbetaprime_t_field(l, k) = ddRtilde_dbetadbetaprime_t_kl;
        arma::mat ddHtilde_dbetadbetaprime_t_kl = fc_ddHtilde_dbetadbetaprime_t_kl(Dtilde_t, ddRtilde_dbetadbetaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\beta_\ell
        ddltilde_dthetadthetaprime_t(dimdelta+k,dimdelta+l) = ddltilde_dthetadthetaprime_t(dimdelta+l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_dbeta_t_l, dHtilde_dbeta_t_k, ddHtilde_dbetadbetaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\beta_\ell
      }
    }
    for(int k = 0; k < dimbeta; k++) {
      for(int l = 0; l < dimdelta; l++) {
        arma::mat dDtilde_ddelta_t_l = dDtilde_ddelta_t_cube.slice(l); // \partial\widetilde{D}_t(\bm\delta) / \partial\delta_\ell
        arma::mat dPsitilde_ddelta_tminus1_l = dPsitilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde\Psi_{t-1}(\bm\delta) / \partial\delta_\ell
        arma::mat dRtilde_dbeta_t_k = dRtilde_dbeta_t_cube.slice(k); // \partial\widetilde{R}_t(\bm\theta) / \partial\beta_k
        arma::mat dRtilde_ddelta_tminus1_l = dRtilde_ddelta_tminus1_cube.slice(l); // \partial\widetilde{R}_{t-1}(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_ddelta_t_l = dHtilde_ddelta_t_cube.slice(l); // \partial\widetilde{H}_t(\bm\theta) / \partial\delta_\ell
        arma::mat dHtilde_dbeta_t_k = dHtilde_dbeta_t_cube.slice(k); // \partial\widetilde{H}_t(\bm\theta) / \partial\beta_k
        arma::mat ddRtilde_dbetaddeltaprime_tminus1_kl = ddRtilde_dbetaddeltaprime_tminus1_field(k, l); // \partial^2\widetilde{R}_{t-1}(\bm\theta) / \partial\beta_k\partial\delta_\ell
        arma::mat ddRtilde_dbetaddeltaprime_t_kl = fc_ddRtilde_dbetaddeltaprime_t_kl(l, k, m, beta_2, dPsitilde_ddelta_tminus1_l, dRtilde_ddelta_tminus1_l, ddRtilde_dbetaddeltaprime_tminus1_kl); // \partial^2\widetilde{R}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
        ddRtilde_dbetaddeltaprime_t_field(k, l) = ddRtilde_dbetaddeltaprime_t_kl;
        arma::mat ddHtilde_dbetaddeltaprime_t_kl = fc_ddHtilde_dbetaddeltaprime_t_kl(Dtilde_t, dDtilde_ddelta_t_l, dRtilde_dbeta_t_k, ddRtilde_dbetaddeltaprime_t_kl); // \partial^2\widetilde{H}_t(\bm\theta) / \partial\beta_k\partial\delta_\ell
        ddltilde_dthetadthetaprime_t(dimdelta+k,l) = ddltilde_dthetadthetaprime_t(l,dimdelta+k) = fc_ddltilde_dthetadthetaprime_t_kl(m, y_t, Htilde_t, dHtilde_ddelta_t_l, dHtilde_dbeta_t_k, ddHtilde_dbetaddeltaprime_t_kl); // \partial^2\widetilde{\ell}_{t}(\bm\theta) / \partial\beta_k\partial\delta_\ell
      }
    }
    // \widehat{\Sigma}_{*, t}
    arma::mat SigmaStarhat_t = ddltilde_dthetadthetaprime_t; // \widehat{\Sigma}_{*, t}
    nSigmaStarhat = nSigmaStarhat + SigmaStarhat_t;
    // \widetilde{\varepsilon}_{t-\Bbbk-1}(\bm\delta), ..., \widetilde{\varepsilon}_{t-2}(\bm\delta), and \widetilde{\varepsilon}_{t-1}(\bm\delta)
    varepsilontilde_tminusBkminus1TOtminus2_mat = varepsilontilde_tminusBkTOtminus1_mat;
    varepsilontilde_tminus1 = varepsilontilde_t;
    // \partial\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_\ell, ..., \partial\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_\ell, and \partial\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_\ell
    dlnhulinetilde_ddelta_tminusBkminus1TOtminus2_cube = dlnhulinetilde_ddelta_tminusBkTOtminus1_cube;
    dlnhulinetilde_ddelta_tminus1_mat = dlnhulinetilde_ddelta_t_mat;
    // \partial^2\ln\widetilde\underline{\bm{h}}_{t-\Bbbk-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell, ..., \partial^2\ln\widetilde\underline{\bm{h}}_{t-2}(\bm\delta) / \partial\delta_k\partial\delta_\ell and \partial^2\ln\widetilde\underline{\bm{h}}_{t-1}(\bm\delta) / \partial\delta_k\partial\delta_\ell
    ddlnhulinetilde_ddeltaddeltaprime_tminusBkminus1TOtminus2_field = ddlnhulinetilde_ddeltaddeltaprime_tminusBkTOtminus1_field;
    ddlnhulinetilde_ddeltaddeltaprime_tminus1_cube = ddlnhulinetilde_ddeltaddeltaprime_t_cube;
    // \widetilde{R}_{t-1}(\bm\theta), \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\delta_\ell, \partial\widetilde{R}_{t-1}(\bm\theta)/\partial\beta_\ell, and \partial^2\widetilde{R}_{t-1}(\bm\theta)/\partial\theta_\k\partial\theta_\ell
    Rtilde_tminus1 = Rtilde_t;
    dRtilde_ddelta_tminus1_cube = dRtilde_ddelta_t_cube;
    dRtilde_dbeta_tminus1_cube = dRtilde_dbeta_t_cube;
    ddRtilde_ddeltaddeltaprime_tminus1_field = ddRtilde_ddeltaddeltaprime_t_field;
    ddRtilde_dbetadbetaprime_tminus1_field = ddRtilde_dbetadbetaprime_t_field;
    ddRtilde_dbetaddeltaprime_tminus1_field = ddRtilde_dbetaddeltaprime_t_field;
  }
  arma::mat Sigmahat = nSigmahat / (n * 1.0); // \widehat{\Sigma} = 1/n \sum_{t=1}^{n} \widehat{\Sigma}_t
  arma::mat SigmaStarhat = nSigmaStarhat / (n * 1.0); // \widehat{\Sigma}_{*} = 1/n \sum_{t=1}^{n} \widehat{\Sigma}_{*, t}
  arma::mat CovMat_thetahat = SigmaStarhat.i() * Sigmahat * SigmaStarhat.i(); // \widehat\Sigma_\star^{-1} \widehat\Sigma \widehat\Sigma_\star^{-1}
  arma::mat transfmat = fc_transfmat_theta_to_vecPhi1(m, r, s, dimtheta); // transformation matrix for trasfering \theta to vec(\Phi_1)
  arma::mat CovMat_vecPhi1 = transfmat * CovMat_thetahat * transfmat.t(); // covariance matrix of vec(\Phi_1)
  arma::vec var_sqrtnvecPhi1 = CovMat_vecPhi1.diag(0); // variances of \sqrt{n} vec(\Phi_1)
  arma::vec var_vecPhi1 = var_sqrtnvecPhi1 / (n * 1.0); // variances of vec(\Phi_1)
  arma::vec ASD_vecPhi1 = sqrt(var_vecPhi1); // ASD of vec(\Phi_1)
  return ASD_vecPhi1;
}


