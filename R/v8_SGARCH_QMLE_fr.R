# library(Rcpp)
# library(RcppArmadillo)
# library("nloptr")

# QMLE
##### \widehat{\bm\theta}_n = \argmin_{\bm\theta \in \Theta} \widetilde{\mathcal{L}}_n (\bm\theta)
##### In general,
##### \bm\theta = (\bm\delta', \bm\beta')'
#####           = (\underline\bm\omega'(m*1), \bm\kappa', \beta_1, \beta_2, \underline\bmr'(m(m-1)/2*1))'
#####           = (\underline\bm\omega'(m*1), \bm\lambda'(r*1), \bm\gamma'(s*1), \bm\varphi'(s*1), \bm{g}_0'(rm^2*1), \bm{g}_1'(sm^2*1), \bm{g}_2'(sm^2*1), \beta_1, \beta_2, \underline\bm{r}'(m(m-1)/2*1))'
##### \bm{g}_0 = (\bm{g}_{01}', ..., \bm{g}_{0r}')', \bm{g}_1 = (\bm{g}_{11}', ..., \bm{g}_{1s}')', \bm{g}_2 = (\bm{g}_{21}', ..., \bm{g}_{2s}')' with 
##### \bm{g}_{0k} = vec(G_{0k}),
##### \bm{g}_{1k} = vec(G_{1k}),
##### \bm{g}_{2k} = vec(G_{2k}).
##### Under low-rank constraints rank(G_{0k}) = 1 and rank(G_{1k}) = rank(G_{2k}) = 2, 
##### \bm\vartheta = (\bm\delta', \bm\beta')'
#####              = (\underline\bm\omega'(m*1), \bm\kappa', \beta_1, \beta_2, \underline\bmr'(m(m-1)/2*1))'
#####              = (\underline\bm\omega'(m*1), \bm\lambda'(r*1), \bm\gamma'(s*1), \bm\varphi'(s*1), \bm{g}_0'(2rm*1), \bm{g}_1'(4sm*1), \bm{g}_2'(4sm*1), \beta_1, \beta_2, \underline\bm{r}'(m(m-1)/2*1))'
##### \bm{g}_0 = (\bm{g}_01', ..., \bm{g}_0r')', \bm{g}_1 = (\bm{g}_11', ..., \bm{g}_1s')', \bm{g}_2 = (\bm{g}_21', ..., \bm{g}_2s')' with 
##### \bm{g}_0k = (\bm{g}_0k1', \bm{g}_0k2')' with G_0k = \bm{g}_0k1 \bm{g}_0k2' such that rank(G_0k) = 1, 
##### \bm{g}_1k = (\bm{g}_1k1', \bm{g}_1k2', \bm{g}_1k3', \bm{g}_1k4')' with G_1k = \bm{g}_1k1 \bm{g}_1k2' + \bm{g}_1k3 \bm{g}_1k4' such that rank(G_1k) <= 2, 
##### \bm{g}_2k = (\bm{g}_2k1', \bm{g}_2k2', \bm{g}_2k3', \bm{g}_2k4')' with G_2k = \bm{g}_2k1 \bm{g}_2k2' + \bm{g}_2k3 \bm{g}_2k4' such that rank(G_2k) <= 2, 
##### dimkappa <- r + 2*s + 2*r*m + 8*s*m
##### dimdelta <- m + dimkappa
##### dimbeta <- 2 + m*(m-1)/2
##### dimvartheta <- dimdelta + dimbeta
# source("sv_v8_SGARCH_QMLE_FFT_fr.R")
# Rcpp::sourceCpp('sv_v8_SGARCH_QMLE_fc.cpp')
## objective function and its gradient
##### objective function is fc_mathcalLtilde_n(n, m, r, s, Bk, vartheta_vec, y_m_n)
fr_lowrank_mathcalLtilde_n <- function(vartheta_vec, n, m, r, s, Bk, y_m_n) { 
  # \widetilde{\mathcal{L}}_{n}(\bm\theta)
  result <- fc_lowrank_mathcalLtilde_n(n, m, r, s, Bk, vartheta_vec, y_m_n)
  return(result)
}
fr_general_mathcalLtilde_n <- function(theta_vec, n, m, r, s, Bk, y_m_n) { 
  # \widetilde{\mathcal{L}}_{n}(\bm\theta)
  result <- fc_general_mathcalLtilde_n(n, m, r, s, Bk, theta_vec, y_m_n)
  return(result)
}
##### the gradient of the objective function is fc_dmathcalLtilde_dvartheta(n, m, r, s, Bk, vartheta_vec, y_m_n)
fr_dmathcalLtilde_dvartheta <- function(vartheta_vec, n, m, r, s, Bk, y_m_n) {
  # \partial\widetilde{\mathcal{L}}_{n}(\bm\theta) / \partial\bm\vartheta
  result <- fc_dmathcalLtilde_dvartheta(n, m, r, s, Bk, vartheta_vec, y_m_n)
  return(result)
}
fr_dmathcalLtilde_dtheta <- function(theta_vec, n, m, r, s, Bk, y_m_n) {
  # \partial\widetilde{\mathcal{L}}_{n}(\bm\theta) / \partial\bm\theta
  result <- fc_dmathcalLtilde_dtheta(n, m, r, s, Bk, theta_vec, y_m_n)
  return(result)
}
## Optimization and estimation
##### estimation of parameters
fr_lowrank_SGARCH_est <- function(n, m, r, s, Bk, y_m_n, vartheta_ini) {
  dimvartheta <- m + (r + 2*s + 2*r*m + 8*s*m) + 2 + m*(m-1)/2
  ## linear constraints ui %*% theta - ci >=0
  ui <- fc_lowrank_linear_constr_ui(m, r, s, dimvartheta)
  ci <- fc_lowrank_linear_constr_ci(m, r, s)
  ## objective function and its gradient
  fr2_lowrank_mathcalLtilde_n <- function(vartheta_vec) {result <- fr_lowrank_mathcalLtilde_n(vartheta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  fr2_dmathcalLtilde_dvartheta <- function(vartheta_vec) {result <- fr_dmathcalLtilde_dvartheta(vartheta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  ## estimation
  varthetahat_list <- constrOptim(theta = vartheta_ini,
                                  f = fr2_lowrank_mathcalLtilde_n,
                                  grad = fr2_dmathcalLtilde_dvartheta,
                                  ui = ui,
                                  ci = ci )
  varthetahat <- varthetahat_list$par
  thetahat <- fc_varthetaTOtheta_general(m, r, s, varthetahat)
  return(thetahat)
}
fr_general_SGARCH_est <- function(n, m, r, s, Bk, y_m_n, theta_ini) {
  dimtheta <- m + (r + 2*s + r*m*m + 2*s*m*m) + 2 + m*(m-1)/2
  ## linear constraints ui %*% theta - ci >=0
  ui <- fc_general_linear_constr_ui(m, r, s, dimtheta)
  ci <- fc_general_linear_constr_ci(m, r, s)
  ## objective function and its gradient
  fr2_general_mathcalLtilde_n <- function(theta_vec) {result <- fr_general_mathcalLtilde_n(theta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  fr2_dmathcalLtilde_dtheta <- function(theta_vec) {result <- fr_dmathcalLtilde_dtheta(theta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  ## estimation
  thetahat_list <- constrOptim(theta = theta_ini,
                               f = fr2_general_mathcalLtilde_n,
                               grad = fr2_dmathcalLtilde_dtheta,
                               ui = ui,
                               ci = ci )
  thetahat <- thetahat_list$par
  return(thetahat)
}
fr_lowrank_SGARCH_est_largem <- function(n, m, r, s, Bk, y_m_n, vartheta_ini) {
  dimvartheta <- m + (r + 2*s + 2*r*m + 8*s*m) + 2 + m*(m-1)/2
  ## linear constraints ui %*% theta - ci >=0
  ui <- fc_lowrank_linear_constr_ui(m, r, s, dimvartheta)
  ci <- fc_lowrank_linear_constr_ci(m, r, s)
  ## objective function and its gradient
  fr2_lowrank_mathcalLtilde_n <- function(vartheta_vec) {result <- fr_lowrank_mathcalLtilde_n(vartheta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  fr2_dmathcalLtilde_dvartheta <- function(vartheta_vec) {result <- fr_dmathcalLtilde_dvartheta(vartheta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  ## estimation
  varthetahat_list <- constrOptim(theta = vartheta_ini,
                                  f = fr2_lowrank_mathcalLtilde_n,
                                  grad = fr2_dmathcalLtilde_dvartheta,
                                  ui = ui,
                                  ci = ci,
                                  control = list(reltol = 1e-5) )
  varthetahat <- varthetahat_list$par
  thetahat <- fc_varthetaTOtheta_general(m, r, s, varthetahat)
  return(thetahat)
}
fr_general_SGARCH_est_largem <- function(n, m, r, s, Bk, y_m_n, theta_ini) {
  dimtheta <- m + (r + 2*s + r*m*m + 2*s*m*m) + 2 + m*(m-1)/2
  ## linear constraints ui %*% theta - ci >=0
  ui <- fc_general_linear_constr_ui(m, r, s, dimtheta)
  ci <- fc_general_linear_constr_ci(m, r, s)
  ## objective function and its gradient
  fr2_general_mathcalLtilde_n <- function(theta_vec) {result <- fr_general_mathcalLtilde_n(theta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  fr2_dmathcalLtilde_dtheta <- function(theta_vec) {result <- fr_dmathcalLtilde_dtheta(theta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  ## estimation
  thetahat_list <- constrOptim(theta = theta_ini,
                               f = fr2_general_mathcalLtilde_n,
                               grad = fr2_dmathcalLtilde_dtheta,
                               ui = ui,
                               ci = ci,
                               control = list(reltol = 1e-5) )
  thetahat <- thetahat_list$par
  return(thetahat)
}
##### estimation of parameters and ASDs
fr_lowrank_SGARCH_est_ASD <- function(n, m, r, s, Bk, y_m_n, vartheta_ini) {
  dimvartheta <- m + (r + 2*s + 2*r*m + 8*s*m) + 2 + m*(m-1)/2
  dimtheta <- m + (r + 2*s + r*m*m + 2*s*m*m) + 2 + m*(m-1)/2
  ## linear constraints ui %*% theta - ci >=0
  ui <- fc_lowrank_linear_constr_ui(m, r, s, dimvartheta)
  ci <- fc_lowrank_linear_constr_ci(m, r, s)
  ## objective function and its gradient
  fr2_lowrank_mathcalLtilde_n <- function(vartheta_vec) {result <- fr_lowrank_mathcalLtilde_n(vartheta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  fr2_dmathcalLtilde_dvartheta <- function(vartheta_vec) {result <- fr_dmathcalLtilde_dvartheta(vartheta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  ## estimation
  varthetahat_list <- constrOptim(theta = vartheta_ini,
                                  f = fr2_lowrank_mathcalLtilde_n,
                                  grad = fr2_dmathcalLtilde_dvartheta,
                                  ui = ui,
                                  ci = ci )
  varthetahat <- varthetahat_list$par
  thetahat <- fc_varthetaTOtheta_general(m, r, s, varthetahat)
  Deltahat_theta_vartheta <- fc_derivative_of_transformation_function(m, r, s, varthetahat)
  ASD_thetahat <- fc_lowrank_ASD_thetahat(n, m, r, s, Bk, thetahat, y_m_n, Deltahat_theta_vartheta)
  transfmat_theta_to_vecPhi1 <- fc_transfmat_theta_to_vecPhi1(m, r, s, dimtheta)
  vecPhi1hat <- transfmat_theta_to_vecPhi1 %*% thetahat
  ASD_vecPhi1hat <- fc_lowrank_ASD_vecPhi1hat(n, m, r, s, Bk, thetahat, y_m_n, Deltahat_theta_vartheta)
  result_list <- list(thetahat=thetahat, ASD_thetahat=ASD_thetahat, vecPhi1hat=vecPhi1hat, ASD_vecPhi1hat=ASD_vecPhi1hat)
  return(result_list)
}
fr_general_SGARCH_est_ASD <- function(n, m, r, s, Bk, y_m_n, theta_ini) {
  dimtheta <- m + (r + 2*s + r*m*m + 2*s*m*m) + 2 + m*(m-1)/2
  ## linear constraints ui %*% theta - ci >=0
  ui <- fc_general_linear_constr_ui(m, r, s, dimtheta)
  ci <- fc_general_linear_constr_ci(m, r, s)
  ## objective function and its gradient
  fr2_general_mathcalLtilde_n <- function(theta_vec) {result <- fr_general_mathcalLtilde_n(theta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  fr2_dmathcalLtilde_dtheta <- function(theta_vec) {result <- fr_dmathcalLtilde_dtheta(theta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  ## estimation
  thetahat_list <- constrOptim(theta = theta_ini,
                               f = fr2_general_mathcalLtilde_n,
                               grad = fr2_dmathcalLtilde_dtheta,
                               ui = ui,
                               ci = ci )
  thetahat <- thetahat_list$par
  ASD_thetahat <- fc_general_ASD_thetahat(n, m, r, s, Bk, thetahat, y_m_n)
  transfmat_theta_to_vecPhi1 <- fc_transfmat_theta_to_vecPhi1(m, r, s, dimtheta)
  vecPhi1hat <- transfmat_theta_to_vecPhi1 %*% thetahat
  ASD_vecPhi1hat <- fc_general_ASD_vecPhi1hat(n, m, r, s, Bk, thetahat, y_m_n)
  result_list <- list(thetahat=thetahat, ASD_thetahat=ASD_thetahat, vecPhi1hat=vecPhi1hat, ASD_vecPhi1hat=ASD_vecPhi1hat)
  return(result_list)
}
##### estimation of parameters and ASDs, in one replication of simulation
fr_lowrank_SGARCH_est_ASD_vec <- function(n, m, r, s, Bk, y_m_n, vartheta_ini) {
  dimvartheta <- m + (r + 2*s + 2*r*m + 8*s*m) + 2 + m*(m-1)/2
  ## linear constraints ui %*% theta - ci >=0
  ui <- fc_lowrank_linear_constr_ui(m, r, s, dimvartheta)
  ci <- fc_lowrank_linear_constr_ci(m, r, s)
  ## objective function and its gradient
  fr2_lowrank_mathcalLtilde_n <- function(vartheta_vec) {result <- fr_lowrank_mathcalLtilde_n(vartheta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  fr2_dmathcalLtilde_dvartheta <- function(vartheta_vec) {result <- fr_dmathcalLtilde_dvartheta(vartheta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  ## estimation
  varthetahat_list <- constrOptim(theta = vartheta_ini,
                                  f = fr2_lowrank_mathcalLtilde_n,
                                  grad = fr2_dmathcalLtilde_dvartheta,
                                  ui = ui,
                                  ci = ci )
  varthetahat <- varthetahat_list$par
  thetahat <- fc_varthetaTOtheta_general(m, r, s, varthetahat)
  Deltahat_theta_vartheta <- fc_derivative_of_transformation_function(m, r, s, varthetahat)
  ASD <- fc_lowrank_ASD_thetahat(n, m, r, s, Bk, thetahat, y_m_n, Deltahat_theta_vartheta)
  result_vec <- c(thetahat, ASD)
  return(result_vec)
}
fr_general_SGARCH_est_ASD_vec <- function(n, m, r, s, Bk, y_m_n, theta_ini) {
  dimtheta <- m + (r + 2*s + r*m*m + 2*s*m*m) + 2 + m*(m-1)/2
  ## linear constraints ui %*% theta - ci >=0
  ui <- fc_general_linear_constr_ui(m, r, s, dimtheta)
  ci <- fc_general_linear_constr_ci(m, r, s)
  ## objective function and its gradient
  fr2_general_mathcalLtilde_n <- function(theta_vec) {result <- fr_general_mathcalLtilde_n(theta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  fr2_dmathcalLtilde_dtheta <- function(theta_vec) {result <- fr_dmathcalLtilde_dtheta(theta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  ## estimation
  thetahat_list <- constrOptim(theta = theta_ini,
                               f = fr2_general_mathcalLtilde_n,
                               grad = fr2_dmathcalLtilde_dtheta,
                               ui = ui,
                               ci = ci )
  thetahat <- thetahat_list$par
  ASD <- fc_general_ASD_thetahat(n, m, r, s, Bk, thetahat, y_m_n)
  result_vec <- c(thetahat, ASD)
  return(result_vec)
}


# BIC and AIC
fr_SGARCH_varthetaest <- function(n, m, r, s, Bk, y_m_n, vartheta_ini) {
  dimvartheta <- m + (r + 2*s + 2*r*m + 8*s*m) + 2 + m*(m-1)/2
  ## linear constraints ui %*% theta - ci >=0
  ui <- fc_lowrank_linear_constr_ui(m, r, s, dimvartheta)
  ci <- fc_lowrank_linear_constr_ci(m, r, s)
  ## objective function and its gradient
  fr2_lowrank_mathcalLtilde_n <- function(vartheta_vec) {result <- fr_lowrank_mathcalLtilde_n(vartheta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  fr2_dmathcalLtilde_dvartheta <- function(vartheta_vec) {result <- fr_dmathcalLtilde_dvartheta(vartheta_vec, n, m, r, s, Bk, y_m_n); return(result)}
  ## estimation
  varthetahat_list <- constrOptim(theta = vartheta_ini,
                                  f = fr2_lowrank_mathcalLtilde_n,
                                  grad = fr2_dmathcalLtilde_dvartheta,
                                  ui = ui,
                                  ci = ci )
  varthetahat <- varthetahat_list$par
  return(varthetahat)
}

fr_BIC <- function(n, m, d, mathcalLtilde_n_thetahat) {
  BIC <- m*n*log(2*3.14) + n*mathcalLtilde_n_thetahat + d*log(n)
  return(BIC)
}

## m=2
fr_lowrank_BIC_mequal2 <- function(n, Bk, y_m_n, vartheta_ini_case1, vartheta_ini_case2, vartheta_ini_case3) {
  varthetahat_case1 <- fr_SGARCH_varthetaest(n, 2, 1, 0, Bk, y_m_n, vartheta_ini_case1)
  varthetahat_case2 <- fr_SGARCH_varthetaest(n, 2, 2, 0, Bk, y_m_n, vartheta_ini_case2)
  varthetahat_case3 <- fr_SGARCH_varthetaest(n, 2, 0, 1, Bk, y_m_n, vartheta_ini_case3)
  mathcalLtilde_n_case1 <- fr_lowrank_mathcalLtilde_n(varthetahat_case1, n, 2, 1, 0, Bk, y_m_n)
  mathcalLtilde_n_case2 <- fr_lowrank_mathcalLtilde_n(varthetahat_case2, n, 2, 2, 0, Bk, y_m_n)
  mathcalLtilde_n_case3 <- fr_lowrank_mathcalLtilde_n(varthetahat_case3, n, 2, 0, 1, Bk, y_m_n)
  BIC_vec <- rep(NA, 3)
  d_case1 <- length(vartheta_ini_case1)
  d_case2 <- length(vartheta_ini_case2)
  d_case3 <- length(vartheta_ini_case3)
  BIC_vec[1] <- fr_BIC(n, 2, d_case1, mathcalLtilde_n_case1)
  BIC_vec[2] <- fr_BIC(n, 2, d_case2, mathcalLtilde_n_case2)
  BIC_vec[3] <- fr_BIC(n, 2, d_case3, mathcalLtilde_n_case3)
  result <- which.min(BIC_vec)
  return(result)
}
fr_general_BIC_mequal2 <- function(n, Bk, y_m_n, theta_ini_case1, theta_ini_case2, theta_ini_case3) {
  thetahat_case1 <- fr_general_SGARCH_est(n, 2, 1, 0, Bk, y_m_n, theta_ini_case1)
  thetahat_case2 <- fr_general_SGARCH_est(n, 2, 2, 0, Bk, y_m_n, theta_ini_case2)
  thetahat_case3 <- fr_general_SGARCH_est(n, 2, 0, 1, Bk, y_m_n, theta_ini_case3)
  mathcalLtilde_n_case1 <- fr_general_mathcalLtilde_n(thetahat_case1, n, 2, 1, 0, Bk, y_m_n)
  mathcalLtilde_n_case2 <- fr_general_mathcalLtilde_n(thetahat_case2, n, 2, 2, 0, Bk, y_m_n)
  mathcalLtilde_n_case3 <- fr_general_mathcalLtilde_n(thetahat_case3, n, 2, 0, 1, Bk, y_m_n)
  BIC_vec <- rep(NA, 3)
  d_case1 <- length(theta_ini_case1)
  d_case2 <- length(theta_ini_case2)
  d_case3 <- length(theta_ini_case3)
  BIC_vec[1] <- fr_BIC(n, 2, d_case1, mathcalLtilde_n_case1)
  BIC_vec[2] <- fr_BIC(n, 2, d_case2, mathcalLtilde_n_case2)
  BIC_vec[3] <- fr_BIC(n, 2, d_case3, mathcalLtilde_n_case3)
  result <- which.min(BIC_vec)
  return(result)
}

## m=5
fr_lowrank_BIC_mequal5_rplus2slessthan2 <- function(n, Bk, y_m_n, vartheta_ini_case1, vartheta_ini_case2, vartheta_ini_case3) {
  varthetahat_case1 <- fr_SGARCH_varthetaest(n, 5, 1, 0, Bk, y_m_n, vartheta_ini_case1)
  varthetahat_case2 <- fr_SGARCH_varthetaest(n, 5, 2, 0, Bk, y_m_n, vartheta_ini_case2)
  varthetahat_case3 <- fr_SGARCH_varthetaest(n, 5, 0, 1, Bk, y_m_n, vartheta_ini_case3)
  mathcalLtilde_n_case1 <- fr_lowrank_mathcalLtilde_n(varthetahat_case1, n, 5, 1, 0, Bk, y_m_n)
  mathcalLtilde_n_case2 <- fr_lowrank_mathcalLtilde_n(varthetahat_case2, n, 5, 2, 0, Bk, y_m_n)
  mathcalLtilde_n_case3 <- fr_lowrank_mathcalLtilde_n(varthetahat_case3, n, 5, 0, 1, Bk, y_m_n)
  BIC_vec <- rep(NA, 3)
  d_case1 <- length(vartheta_ini_case1)
  d_case2 <- length(vartheta_ini_case2)
  d_case3 <- length(vartheta_ini_case3)
  BIC_vec[1] <- fr_BIC(n, 5, d_case1, mathcalLtilde_n_case1)
  BIC_vec[2] <- fr_BIC(n, 5, d_case2, mathcalLtilde_n_case2)
  BIC_vec[3] <- fr_BIC(n, 5, d_case3, mathcalLtilde_n_case3)
  result <- which.min(BIC_vec)
  return(result)
}
fr_general_BIC_mequal5_rplus2slessthan2 <- function(n, Bk, y_m_n, theta_ini_case1, theta_ini_case2, theta_ini_case3) {
  thetahat_case1 <- fr_general_SGARCH_est(n, 5, 1, 0, Bk, y_m_n, theta_ini_case1)
  thetahat_case2 <- fr_general_SGARCH_est(n, 5, 2, 0, Bk, y_m_n, theta_ini_case2)
  thetahat_case3 <- fr_general_SGARCH_est(n, 5, 0, 1, Bk, y_m_n, theta_ini_case3)
  mathcalLtilde_n_case1 <- fr_general_mathcalLtilde_n(thetahat_case1, n, 5, 1, 0, Bk, y_m_n)
  mathcalLtilde_n_case2 <- fr_general_mathcalLtilde_n(thetahat_case2, n, 5, 2, 0, Bk, y_m_n)
  mathcalLtilde_n_case3 <- fr_general_mathcalLtilde_n(thetahat_case3, n, 5, 0, 1, Bk, y_m_n)
  BIC_vec <- rep(NA, 3)
  d_case1 <- length(theta_ini_case1)
  d_case2 <- length(theta_ini_case2)
  d_case3 <- length(theta_ini_case3)
  BIC_vec[1] <- fr_BIC(n, 5, d_case1, mathcalLtilde_n_case1)
  BIC_vec[2] <- fr_BIC(n, 5, d_case2, mathcalLtilde_n_case2)
  BIC_vec[3] <- fr_BIC(n, 5, d_case3, mathcalLtilde_n_case3)
  result <- which.min(BIC_vec)
  return(result)
}

# ## m=20
# fr_lowrank_BIC_mequal20_rplus2slessthan2 <- function(n, Bk, y_m_n, vartheta_ini_case1, vartheta_ini_case2, vartheta_ini_case3) {
#   varthetahat_case1 <- fr_SGARCH_varthetaest(n, 20, 1, 0, Bk, y_m_n, vartheta_ini_case1)
#   varthetahat_case2 <- fr_SGARCH_varthetaest(n, 20, 2, 0, Bk, y_m_n, vartheta_ini_case2)
#   varthetahat_case3 <- fr_SGARCH_varthetaest(n, 20, 0, 1, Bk, y_m_n, vartheta_ini_case3)
#   mathcalLtilde_n_case1 <- fr_lowrank_mathcalLtilde_n(varthetahat_case1, n, 20, 1, 0, Bk, y_m_n)
#   mathcalLtilde_n_case2 <- fr_lowrank_mathcalLtilde_n(varthetahat_case2, n, 20, 2, 0, Bk, y_m_n)
#   mathcalLtilde_n_case3 <- fr_lowrank_mathcalLtilde_n(varthetahat_case3, n, 20, 0, 1, Bk, y_m_n)
#   BIC_vec <- rep(NA, 3)
#   d_case1 <- length(vartheta_ini_case1)
#   d_case2 <- length(vartheta_ini_case2)
#   d_case3 <- length(vartheta_ini_case3)
#   BIC_vec[1] <- fr_BIC(n, 20, d_case1, mathcalLtilde_n_case1)
#   BIC_vec[2] <- fr_BIC(n, 20, d_case2, mathcalLtilde_n_case2)
#   BIC_vec[3] <- fr_BIC(n, 20, d_case3, mathcalLtilde_n_case3)
#   result <- which.min(BIC_vec)
#   return(result)
# }

