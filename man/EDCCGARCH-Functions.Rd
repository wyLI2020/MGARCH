\name{EDCCGARCH-Functions}
\alias{fc_derivative_of_transformation_function}
\alias{fc_dmathcalLtilde_dtheta}
\alias{fc_dmathcalLtilde_dvartheta}
\alias{fc_general_ASD_thetahat}
\alias{fc_general_ASD_vecPhi1hat}
\alias{fc_general_linear_constr_ci}
\alias{fc_general_linear_constr_ui}
\alias{fc_general_mathcalLtilde_n}
\alias{fc_lowrank_ASD_thetahat}
\alias{fc_lowrank_ASD_vecPhi1hat}
\alias{fc_lowrank_linear_constr_ci}
\alias{fc_lowrank_linear_constr_ui}
\alias{fc_lowrank_mathcalLtilde_n}
\alias{fc_transfmat_theta_to_vecPhi1}
\alias{fc_varthetaTOtheta_general}
\alias{fr_BIC}
\alias{fr_FFT}
\alias{fr_SGARCH_varthetaest}
\alias{fr_dmathcalLtilde_dtheta}
\alias{fr_dmathcalLtilde_dvartheta}
\alias{fr_general_BIC_mequal2}
\alias{fr_general_BIC_mequal5_rplus2slessthan2}
\alias{fr_general_SGARCH_est}
\alias{fr_general_SGARCH_est_ASD}
\alias{fr_general_SGARCH_est_ASD_vec}
\alias{fr_general_SGARCH_est_largem}
\alias{fr_general_mathcalLtilde_n}
\alias{fr_lowrank_BIC_mequal2}
\alias{fr_lowrank_BIC_mequal5_rplus2slessthan2}
\alias{fr_lowrank_SGARCH_est}
\alias{fr_lowrank_SGARCH_est_ASD}
\alias{fr_lowrank_SGARCH_est_ASD_vec}
\alias{fr_lowrank_SGARCH_est_largem}
\alias{fr_lowrank_mathcalLtilde_n}
\alias{fr_general_SGARCH_est_ASD}
\alias{fr_lowrank_SGARCH_est_ASD}
\title{Estimation of the efficient multivariate volatility model}
\description{
  Estimation of the efficient multivariate volatility model
}
\usage{
  fr_general_SGARCH_est_ASD(n, m, r, s, Bk, y_m_n, theta_ini);
  fr_lowrank_SGARCH_est_ASD(n, m, r, s, Bk, y_m_n, vartheta_ini);
}
\arguments{
  \item{n}{integer, the time dimension}
  \item{m}{integer, the dimension of multivariate time series}
  \item{r}{integer, the number of nonzero real eigenvalues}
  \item{s}{integer, the number of conjugate pairs of nonzero complex eigenvalues}
  \item{Bk}{integer, the number of lagged innovations in calculating the sample correlation matrix}
  \item{y_m_n}{(m, n) matrix, response}
  \item{theta_ini}{initial value}
  \item{vartheta_ini}{initial value}
}
\value{
  \code{fr_general_SGARCH_est_ASD} and \code{fr_lowrank_SGARCH_est_ASD} return lists.
}
\details{
  With the preproccess in following "Examples", we can get the estimators of the efficient multivariate volatility model.
}
\examples{
  library(EDCCGARCH)
  data("y_5_times_2581")
  m = 5; nplus = 2581; n = nplus; y_m_n <- y_m5_2581
  
  r=2; s=0; Bk=m
  omegauline_ini <- rep(1.5, m)
  lambda_1_ini <- 0.8; lambda_2_ini <- 0.7; lambda_ini <- c(lambda_1_ini, lambda_2_ini)
  set.seed(123); g_011_ini <- runif(m, min = 0.2, max = 0.3); g_012_ini <- runif(m, min = 0.02, max = 0.03)
  set.seed(123); g_021_ini <- runif(m, min = 0.3, max = 0.6); g_022_ini <- runif(m, min = 0.03, max = 0.06001)
  g_01_ini <- c(g_011_ini, g_012_ini); g_02_ini <- c(g_021_ini, g_022_ini); g_0_ini <- c(g_01_ini, g_02_ini)
  % G_01_ini <- g_011_ini %*% t(g_012_ini); G_02_ini <- g_021_ini %*% t(g_022_ini)
  G_01_ini <- tcrossprod(g_011_ini, g_012_ini); G_02_ini <- tcrossprod(g_021_ini, g_022_ini)
  beta_1_ini <- 0.1; beta_2_ini = 0.8
  Ruline_ini <- matrix(0.5, nrow = m, ncol = m); diag(Ruline_ini) <- 1.0
  ruline_ini <- Ruline_ini[lower.tri(Ruline_ini)]
  kappa_ini <- c(lambda_ini, g_0_ini)
  vartheta_ini <- c(omegauline_ini, lambda_ini, g_0_ini, beta_1_ini, beta_2_ini, ruline_ini)
  theta_ini <- c(omegauline_ini, lambda_ini, as.vector(G_01_ini), as.vector(G_02_ini), beta_1_ini, beta_2_ini, ruline_ini)
  
  list_general_est_asd <- fr_general_SGARCH_est_ASD(n, m, r, s, Bk, y_m_n, theta_ini)
  results_thetahat_general_table <- matrix(NA, nrow = length(theta_ini), ncol = 4)
  results_thetahat_general_table[,1] <- list_general_est_asd$thetahat
  results_thetahat_general_table[,2] <- list_general_est_asd$ASD_thetahat
  results_thetahat_general_table[,3] <- results_thetahat_general_table[,1] / results_thetahat_general_table[,2]
  results_thetahat_general_table[,4] <- 2*(1-pnorm(abs(results_thetahat_general_table[,1] / results_thetahat_general_table[,2])))
  colnames(results_thetahat_general_table) <- c("est", "ASD", "z", "p-value")
  rownames(results_thetahat_general_table) <- c(paste0("omegauline", 1:m), paste0("lambda", 1:r), rep("G_01", m^2), rep("G_02", m^2), paste0("beta", 1:2), rep("ruline", m*(m-1)/2))
  results_vecPhi1hat_general_table <- matrix(NA, nrow = m^2, ncol = 4)
  results_vecPhi1hat_general_table[,1] <- list_general_est_asd$vecPhi1hat
  results_vecPhi1hat_general_table[,2] <- list_general_est_asd$ASD_vecPhi1hat
  results_vecPhi1hat_general_table[,3] <- results_vecPhi1hat_general_table[,1] / results_vecPhi1hat_general_table[,2]
  results_vecPhi1hat_general_table[,4] <- 2*(1-pnorm(abs(results_vecPhi1hat_general_table[,1] / results_vecPhi1hat_general_table[,2])))
  colnames(results_vecPhi1hat_general_table) <- c("est", "ASD", "z", "p-value")
  rownames(results_vecPhi1hat_general_table) <- rep("Phi_1", m^2)
}
% \references{
  % See the documentation for Armadillo, and RcppArmadillo, for more details.
% }
\author{Wenyu Li}

