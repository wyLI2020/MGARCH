# MGARCH

## Installation

```R
#install.packages("devtools")
library(devtools)
install_github("wyLI2020/MGARCH")
```

## Usage

```R
# QMLE with general coefficient matrices
fr_general_SGARCH_est_ASD(n, m, r, s, Bk, y_m_n, theta_ini)
# QMLE with low-rank constraints on coefficient matrices
fr_lowrank_SGARCH_est_ASD(n, m, r, s, Bk, y_m_n, vartheta_ini)
```

- **n**: integer, the time dimension
- **m**: integer, the dimension of multivariate time series
- **r**: integer, the number of nonzero real eigenvalues with $r+2s \leq m$
- **s**: integer, the number of conjugate pairs of nonzero complex eigenvalues with $r+2s \leq m$
- **Bk**: integer, the number of lagged innovations in calculating the sample correlation matrix
- **y_m_n**:  $(m, n)$ matrix, response
- **theta_ini**: initial value of $\mathbb{\theta}$
- **vartheta_ini**:  initial value of $\mathbb{\vartheta}$

## Example

```R
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
```

