# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

fc_lowrank_linear_constr_ui <- function(m, r, s, dimvartheta) {
    .Call(`_EDCCGARCH_fc_lowrank_linear_constr_ui`, m, r, s, dimvartheta)
}

fc_lowrank_linear_constr_ci <- function(m, r, s) {
    .Call(`_EDCCGARCH_fc_lowrank_linear_constr_ci`, m, r, s)
}

fc_lowrank_mathcalLtilde_n <- function(n, m, r, s, Bk, vartheta_vec, y_m_n) {
    .Call(`_EDCCGARCH_fc_lowrank_mathcalLtilde_n`, n, m, r, s, Bk, vartheta_vec, y_m_n)
}

fc_dmathcalLtilde_dvartheta <- function(n, m, r, s, Bk, vartheta_vec, y_m_n) {
    .Call(`_EDCCGARCH_fc_dmathcalLtilde_dvartheta`, n, m, r, s, Bk, vartheta_vec, y_m_n)
}

fc_varthetaTOtheta_general <- function(m, r, s, varthetahat) {
    .Call(`_EDCCGARCH_fc_varthetaTOtheta_general`, m, r, s, varthetahat)
}

fc_derivative_of_transformation_function <- function(m, r, s, vartheta_vec) {
    .Call(`_EDCCGARCH_fc_derivative_of_transformation_function`, m, r, s, vartheta_vec)
}

fc_lowrank_ASD_thetahat <- function(n, m, r, s, Bk, thetahat, y_m_n, Deltahat_theta_vartheta) {
    .Call(`_EDCCGARCH_fc_lowrank_ASD_thetahat`, n, m, r, s, Bk, thetahat, y_m_n, Deltahat_theta_vartheta)
}

fc_general_linear_constr_ui <- function(m, r, s, dimtheta) {
    .Call(`_EDCCGARCH_fc_general_linear_constr_ui`, m, r, s, dimtheta)
}

fc_general_linear_constr_ci <- function(m, r, s) {
    .Call(`_EDCCGARCH_fc_general_linear_constr_ci`, m, r, s)
}

fc_general_mathcalLtilde_n <- function(n, m, r, s, Bk, theta_vec, y_m_n) {
    .Call(`_EDCCGARCH_fc_general_mathcalLtilde_n`, n, m, r, s, Bk, theta_vec, y_m_n)
}

fc_dmathcalLtilde_dtheta <- function(n, m, r, s, Bk, theta_vec, y_m_n) {
    .Call(`_EDCCGARCH_fc_dmathcalLtilde_dtheta`, n, m, r, s, Bk, theta_vec, y_m_n)
}

fc_general_ASD_thetahat <- function(n, m, r, s, Bk, thetahat, y_m_n) {
    .Call(`_EDCCGARCH_fc_general_ASD_thetahat`, n, m, r, s, Bk, thetahat, y_m_n)
}

fc_transfmat_theta_to_vecPhi1 <- function(m, r, s, dimtheta) {
    .Call(`_EDCCGARCH_fc_transfmat_theta_to_vecPhi1`, m, r, s, dimtheta)
}

fc_lowrank_ASD_vecPhi1hat <- function(n, m, r, s, Bk, thetahat, y_m_n, Deltahat_theta_vartheta) {
    .Call(`_EDCCGARCH_fc_lowrank_ASD_vecPhi1hat`, n, m, r, s, Bk, thetahat, y_m_n, Deltahat_theta_vartheta)
}

fc_general_ASD_vecPhi1hat <- function(n, m, r, s, Bk, thetahat, y_m_n) {
    .Call(`_EDCCGARCH_fc_general_ASD_vecPhi1hat`, n, m, r, s, Bk, thetahat, y_m_n)
}

