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
  data("y_5_times_2456")
  data("theta_ini_m5_r1s0")
  m = 5; n = 2456; y_m_n <- y_m5_2456
  
  r=1; s=0; Bk=m
  list_general_est_asd <- fr_general_SGARCH_est_ASD(n, m, r, s, Bk, y_m_n, theta_ini)
  
  results_thetahat_general_table <- matrix(NA, nrow = length(theta_ini), ncol = 4)
  results_thetahat_general_table[,1] <- list_general_est_asd$thetahat
  results_thetahat_general_table[,2] <- list_general_est_asd$ASD_thetahat
  results_thetahat_general_table[,3] <- results_thetahat_general_table[,1] / results_thetahat_general_table[,2]
  results_thetahat_general_table[,4] <- 2*(1-pnorm(abs(results_thetahat_general_table[,1] / results_thetahat_general_table[,2])))
  colnames(results_thetahat_general_table) <- c("est", "ASD", "z", "p-value")
  rownames(results_thetahat_general_table) <- c(paste0("omegauline", 1:m), paste0("lambda", 1:r), rep("G_01", m^2), paste0("beta", 1:2), rep("ruline", m*(m-1)/2))
  
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

