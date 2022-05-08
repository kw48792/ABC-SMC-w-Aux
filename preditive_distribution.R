rm(list = ls())
library(readr)
library(sgd)
library(deSolve)
library(MASS)
library(mvtnorm)
library(tmvtnorm)
source('Problem_SMC.R')
library(doParallel)
library(foreach)
library(tibble)
library(latex2exp)



simulate <- function(path_parm, path_weight, total_particles, N){
  sample_size <- total_particles / 2 #Number of particles
  results_whole <- read.table(path_parm, head=FALSE, sep = ',')
  results_weight <- read.table(path_weight, head=FALSE, sep = ',')
  G <- dim(results_whole)[1] / total_particles
  pred_err_aux <- matrix(0, ncol=2*n_time, nrow=G)
  pred_traj_aux <- matrix(0, ncol=2*n_time, nrow=N)
  # n_para <- length(results[,1])

  results_g <- results_whole[(total_particles*(G-1)+1):(total_particles*(G-1)+sample_size),]
  results_weight_g <- results_weight[(total_particles*(G-1)+1):(total_particles*(G-1)+sample_size),]
  r_g_ps_aux <- results_g[,1]
  k_s_ps_aux <- results_g[,2]
  k_c_ps_aux <- results_g[,3]
  r_d_ps_aux <- results_g[,4]
  # r_d_ps <- matrix(0.005,200,12)
  sd_rho_ps_aux <- results_g[,5]
  sd_I_ps_aux <- results_g[,6]

  for (s in 1:N) {
    # print(results_weight_g)
    j <- sample.int(sample_size, 1, replace=TRUE, prob=results_weight_g)

    r_g_aug_ps_aux <- c(rep(r_g_ps_aux[j],10))
    k_s_aug_ps_aux <- c(rep(k_s_ps_aux[j],10))
    k_c_aug_ps_aux <- c(rep(k_c_ps_aux[j],10))
    r_d_aug_ps_aux <- c(rep(r_d_ps_aux[j],10))
    sd_rho_aug_ps_aux <- c(rep(sd_rho_ps_aux[j],10))
    sd_I_aug_ps_aux <- c(rep(sd_I_ps_aux[j],10))
    pred_traj_aux[s,]<-run_model_t(r_g_aug_ps_aux, k_s_aug_ps_aux, k_c_aug_ps_aux, r_d_aug_ps_aux, sd_rho_aug_ps_aux, sd_I_aug_ps_aux, 1, n_time)
    # Calculate distances
  }
  rho <- pred_traj_aux[,seq(1,2*n_time,2)]
  I <- pred_traj_aux[,seq(2,2*n_time,2)]
  return(list(rho=rho, I=I))
}

simulate_all <- function (base, method='_aux', total_particles=400, N=1000) {
  rho_df <- c()
  I_df <- c()
  for (m in 1:30) {
    path_parm <- paste(base, 'results_para', method, m, '.txt', sep='')
    path_weight <- paste(base, 'results_weight', method, m, '.txt', sep='')
    print("===================")
    print(path_parm)
    print(path_weight)
    result <- simulate(path_parm, path_weight, total_particles, N)
    rho_df <- cbind(rho_df, result$rho[,11])
    I_df <- cbind(I_df, result$I[,11])
  }
  colnames(rho_df) <- apply(array(1:30, dim=c(1,30)), 1:2, function(x) paste("r", x, sep=''))
  colnames(I_df) <- apply(array(1:30, dim=c(1,30)), 1:2, function(x) paste("r", x, sep=''))
  return (list(rho_df = rho_df, I_df = I_df))
}

################################################################################
###################### True predictive posterior distribution #################
get_true_data <- function (standard_deviation) {
  data_whole <- read.csv("simulation_data_noise_03_01.csv", header=F)
  data_whole <- t(data_whole[,-1])
  rownames(data_whole) <- NULL
  n_time <- length(data_whole[,1])
  tname <- data_whole[,1]

  act_time <- c(40,40,40)

  n_factor <- 2
  n_time <- 11
  n  <- 5000

  r_g <- 0.057
  k_s <- 3.4
  k_c <- 2.6
  r_d <- 0.005
  sd_rho <- standard_deviation
  sd_I <- standard_deviation

  r_g_aug <- c(rep(r_g[1],10))
  k_s_aug <- c(rep(k_s[1],10))
  k_c_aug <- c(rep(k_c[1],10))
  r_d_aug <- c(rep(r_d[1],10))
  sd_rho_aug <- c(rep(sd_rho[1],10))
  sd_I_aug <- c(rep(sd_I[1],10))


  D_star <- matrix(0, ncol=2*n_time,nrow=n)
  for (j in 1:n) {
    D_star[j,]<-run_model_t(r_g_aug, k_s_aug, k_c_aug, r_d_aug, sd_rho_aug, sd_I_aug, 1, n_time)
    # Calculate distances
  }
  true_rho <- data.frame(D_star[,21])
  true_I <- data.frame(D_star[,22])
  colnames(true_rho) <- 'value'
  true_rho['Var2'] <- 'True Model'
  colnames(true_I) <- 'value'
  true_I['Var2'] <- 'True Model'
  return (list(true_rho=true_rho, true_I=true_I))
}

################################################################################
true_rho <- c()
true_I <- c()
for (noise in c(0.1, 0.2)) {
  result <- get_true_data(noise)
  for (m in c(3,6,20)) {
    rho_samples <- result$true_rho
    I_samples <- result$true_I
    rho_samples['noise'] <- paste('v=', noise, sep='')
    I_samples['noise'] <- paste('v=', noise, sep='')
    rho_samples['m'] <- paste("m=", m, sep='')
    I_samples['m'] <- paste("m=", m, sep='')
    true_rho <- rbind(true_rho, rho_samples)
    true_I <- rbind(true_I, I_samples)
  }
}

data_rho <- c()
data_I <- c()
for (m in c(3,6,20)) {
  for (noise in c(0.1, 0.2)) {
    # m <- 3
    # noise <- 0.2
    base <- paste('data/m', m, 'n', sub("\\D+", "", noise, perl = TRUE), '/', sep='')
    result <- simulate_all(base, method='')
    rho_df <- result$rho
    I_df <- result$I

    rho_dat<- melt(rho_df)
    rho_dat['m'] <- paste("m=", m, sep='')
    rho_dat['noise'] <- paste('v=', noise, sep='')
    I_dat <- melt(I_df)
    I_dat['m'] <- paste("m=", m, sep='')
    I_dat['noise'] <- paste('v=', noise, sep='')
    data_rho <- rbind(data_rho, rho_dat)
    data_I <- rbind(data_I, I_dat)
  }
}

data_rho_selected <- data_rho[data_rho$Var2 %in% c("r1", "r7", "r3", "r4", "r8", "r6"),]
data_I_selected <- data_I[data_I$Var2 %in% c("r1", "r7", "r3", "r4", "r8", "r6"),]

g <-ggplot(data_rho_selected, aes(x=value, fill=Var2)) + geom_density(color=NA, alpha=0.4) +
    geom_density(data=true_rho, aes(x=value), lwd=1.2, linetype=2, fill=NA, color='black', alpha=0.9) +
    ylab("Density") +
    xlab('Cell Density') +
    facet_grid(noise ~ factor(m, levels=c('m=3','m=6','m=20')), scales="free")
g + scale_fill_discrete(name = "Macro Replication") + theme(legend.position = "bottomright", strip.text = element_text(size = 20,face="bold"), axis.title=element_text(size=20,face="bold"))

g <-ggplot(data_I_selected, aes(x=value, fill=Var2)) + geom_density(color=NA, alpha=0.4) +
    geom_density(data=true_I, aes(x=value), lwd=1.2, linetype=2, fill=NA, color='black', alpha=0.9) +
    ylab("Density") +
    xlab('Inhibitor Concentration') +
    facet_grid(noise ~ factor(m, levels=c('m=3','m=6','m=20')), scales="free")
g + scale_fill_discrete(name = "Macro Replication") + theme(legend.position = "bottomright", strip.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold"))
