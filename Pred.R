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

setwd(getwd())
set.seed(122) 

R <- 3
data_whole = read.csv("simulation_data_noise_10.csv", header=F)
data_whole <- t(data_whole[,-1])
rownames(data_whole) <- NULL
n_time <- length(data_whole[,1])
ttname = data_whole[,1]
data <- matrix(0,n_time,2*R)
data[,seq(1,2*R,2)] <- data_whole[,1] 
# data[,seq(2,2*R,2)] <- data_whole[,seq(mrep*R-1,mrep*R+1,1)]

# data_test = read.csv("simulation_data_noise_03_01_batch1_3.csv", header=T)

rho = data.frame()
# R <- length(data[1,])/2
for (i in 1:R) {
  rho = rbind(rho, c(t(data[,(2*i)])))
}

names(rho) = ttname
data_new = matrix(0,R,length(data[,1])*2)

for (i in 1:length(rho)) {
  data_new[,(2*i-1)] = rho[,i]
  data_new[,(2*i)] = (0*rho)[,i]
}

Pred_distance <-rep(0, 20)
# r_g_aug <- append(append(r_g, 0, after = 2),0, after = 6)
# k_s_aug <- append(append(k_s, 0, after = 2),0, after = 6)
# k_c_aug <- append(append(k_c, 0, after = 2),0, after = 6)
# r_d_aug <- append(append(r_d, 0, after = 2),0, after = 6)
# sd_rho_aug <- append(append(sd_rho, 0, after = 2),0, after = 6)
# sd_I_aug <- append(append(sd_I, 0, after = 2),0, after = 6)
# r_g_aug <- c(rep(r_g[1],6),rep(r_g[2],4))
# k_s_aug <- c(rep(k_s[1],6),rep(k_s[2],4))
# k_c_aug <- c(rep(k_c[1],6),rep(k_c[2],4))
# r_d_aug <- c(rep(r_d[1],6),rep(r_d[2],4))
# sd_rho_aug <- c(rep(sd_rho[1],6),rep(sd_rho[2],4))
# sd_I_aug <- c(rep(sd_I[1],6),rep(sd_I[2],4))
r_g <- 0.057
k_s <- 3.4
k_c <- 2.6
r_d <- 0.005
sd_rho <- 0.1
sd_I <- 0.1


r_g_aug <- c(rep(r_g[1],10))
k_s_aug <- c(rep(k_s[1],10))
k_c_aug <- c(rep(k_c[1],10))
r_d_aug <- c(rep(r_d[1],10))
sd_rho_aug <- c(rep(sd_rho[1],10))
sd_I_aug <- c(rep(sd_I[1],10))

n_factor <- 2
n_time <- 11
n  <- 60000

N = n_factor*n_time

tname = ttname
xname = paste('X', 1:(n_factor), sep='')
c_names = c(outer(xname, tname, FUN=paste, sep='_'))

data_new = data.frame(data_new)
names(data_new) <- c_names
dat = data_new
# act_time = c(6,18,22)
act_time = c(40,40,40)
tname[act_time]
I_act_5 = c(0,0,0.3,0,0,0,0,0,0,0,0)
I_act_16 = c(0,0,0,0,0,0,0,0,0,0,0)
I_act_18 = c(0,0,0,0,0,0,1,1,0.4,0.3,0)
I_act = rbind(I_act_5,I_act_16,I_act_18)

rho_act_5 = c(0,0,0,0,0,0,0,0,0,0,0)
rho_act_16 = c(0,0,0,0.65,2,4,0,0,0,0,0)
rho_act_18 = c(0,0,0,0,0,0,3.1,0.5,3.5,0.5,0.5)
rho_act = rbind(rho_act_5,rho_act_16,rho_act_18)


D_star <- matrix(0, ncol=2*n_time,nrow=n)
for (j in 1:n) {
  D_star[j,]<-run_model_t(r_g_aug, k_s_aug, k_c_aug, r_d_aug, sd_rho_aug, sd_I_aug, 1, n_time)
  # Calculate distances 
}
# write.table(D_star[,seq(1,2*n_time,2)],"rho_noise_001.txt",sep=",",append=TRUE,eol="\n",
#             row.names=FALSE,col.names=FALSE);
# write.table(D_star[,seq(2,2*n_time,2)],"I_noise_001.txt",sep=",",append=TRUE,eol="\n",
#             row.names=FALSE,col.names=FALSE);


# true_traj <- apply(D_star,2,mean)
# #####Time
# time_consume <- rep(0,30)
# time_consume_aux <- rep(0,30)
# for (i in 1:30) {
#   results_time <- read.table(paste0("results_time",i,".txt"), head=FALSE, sep = ',')
#   time_consume[i] <- sum(results_time)
#   results_time_aux <- read.table(paste0("results_time_aux",i,".txt"), head=FALSE, sep = ',')
#   time_consume_aux[i] <- sum(results_time_aux)
# }
# write.table(time_consume,"time_total.txt",sep=",",append=TRUE,eol="\n",
#             row.names=FALSE,col.names=FALSE);
# write.table(time_consume_aux,"time_total_aux.txt",sep=",",append=TRUE,eol="\n",
#             row.names=FALSE,col.names=FALSE);

# G <- 10

sample_size <- 400 #Number of particles
sample_mrep <- 400/2*10 #Number of simulation per mrep
# n_para <- length(results[,1])
for (i in 1:30) {
  # pred_err <- matrix(0, ncol=2*n_time,nrow=G)
  pred_traj <- matrix(0, ncol=2*n_time,nrow=sample_size/2*10)
  results_whole <- read.table(paste0("results_para",i,".txt"), head=FALSE, sep = ',')
  results <- results_whole[(length(results_whole[,1])-sample_size+1):(length(results_whole[,1])-sample_size/2),]
  r_g_ps <- results[,1]
  k_s_ps <- results[,2]
  k_c_ps <- results[,3]
  r_d_ps <- results[,4]
  # r_d_ps <- matrix(0.005,200,12)
  sd_rho_ps <- results[,5]
  sd_I_ps <- results[,6]
  
  for (j in 1:sample_size/2) {
    r_g_aug_ps <- c(rep(r_g_ps[j],10))
    k_s_aug_ps <- c(rep(k_s_ps[j],10))
    k_c_aug_ps <- c(rep(k_c_ps[j],10))
    r_d_aug_ps <- c(rep(r_d_ps[j],10))
    sd_rho_aug_ps <- c(rep(sd_rho_ps[j],10))
    sd_I_aug_ps <- c(rep(sd_I_ps[j],10))
    for (r in 1:10) {
      pred_traj[((j-1)*10+r),]<-run_model_t(r_g_aug_ps, k_s_aug_ps, k_c_aug_ps, r_d_aug_ps, sd_rho_aug_ps, sd_I_aug_ps, 1, n_time)
    }
    
  }
  
  # Calculate distances 
  # pred_err <- apply(abs(t(t(pred_traj) - true_traj)),2,mean)
  # rho_err <- pred_err[seq(1,2*n_time,2)]
  # rho_err_total <- sum(rho_err)
  # I_err <- pred_err[seq(2,2*n_time,2)]
  # I_err_total <- sum(I_err)
  
  #Calculate KS stats
  ks_rho <- ks.test(pred_traj[,(n_time*2-1)],D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2-1)])$statistic
  ks_I <- ks.test(pred_traj[,(n_time*2)],D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2)])$statistic
  
  #Calculate PI
  upp_rho <- quantile(pred_traj[,(n_time*2-1)],0.975)
  low_rho <- quantile(pred_traj[,(n_time*2-1)],0.025)
  width_rho <- upp_rho - low_rho
  PI_rho <- length(which((D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2-1)]<upp_rho) & (D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2-1)] > low_rho)))/sample_mrep
  
  upp_I <- quantile(pred_traj[,(n_time*2)],0.975)
  low_I <- quantile(pred_traj[,(n_time*2)],0.025)
  width_I <- upp_I - low_I
  PI_I <- length(which((D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2)]<upp_I) & (D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2)] > low_I)))/sample_mrep
  
  
  ####with auxiliary model
  results_whole_aux <- read.table(paste0("results_para_aux",i,".txt"), head=FALSE, sep = ',')
  # results_weight_aux <- read.table("results_weight_aux3.txt", head=FALSE, sep = ',')
  # pred_err_aux <- matrix(0, ncol=2*n_time,nrow=G)
  pred_traj_aux <- matrix(0, ncol=2*n_time,nrow=sample_size/2*10)
  # n_para <- length(results[,1])
  results_aux <- results_whole_aux[(length(results_whole_aux[,1])-sample_size+1):(length(results_whole_aux[,1])-sample_size/2),]
  # weight_aux <- results_weight_aux[(200*(g-1)+1):(200*(g-1)+sample_size),]
  r_g_ps_aux <- results_aux[,1]
  k_s_ps_aux <- results_aux[,2]
  k_c_ps_aux <- results_aux[,3]
  r_d_ps_aux <- results_aux[,4]
  # r_d_ps <- matrix(0.005,200,12)
  sd_rho_ps_aux <- results_aux[,5]
  sd_I_ps_aux <- results_aux[,6]
  
  for (j in 1:sample_size/2) {
    r_g_aug_ps_aux <- c(rep(r_g_ps_aux[j],10))
    k_s_aug_ps_aux <- c(rep(k_s_ps_aux[j],10))
    k_c_aug_ps_aux <- c(rep(k_c_ps_aux[j],10))
    r_d_aug_ps_aux <- c(rep(r_d_ps_aux[j],10))
    sd_rho_aug_ps_aux <- c(rep(sd_rho_ps_aux[j],10))
    sd_I_aug_ps_aux <- c(rep(sd_I_ps_aux[j],10))
    for (r in 1:10) {
      pred_traj_aux[((j-1)*10+r),]<-run_model_t(r_g_aug_ps_aux, k_s_aug_ps_aux, k_c_aug_ps_aux, r_d_aug_ps_aux, sd_rho_aug_ps_aux, sd_I_aug_ps_aux, 1, n_time)
    }   
  }
  # Calculate distances 
  # pred_err_aux <- apply(abs(t(t(pred_traj_aux) - true_traj)),2,mean)
  # rho_err_aux <- pred_err_aux[seq(1,2*n_time,2)]
  # rho_err_total_aux <- sum(rho_err_aux)
  # I_err_aux <- pred_err_aux[seq(2,2*n_time,2)]
  # I_err_total_aux <- sum(I_err_aux)
  
  #Calculate KS stats
  ks_rho_aux <- ks.test(pred_traj_aux[,(n_time*2-1)],D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2-1)])$statistic
  ks_I_aux <- ks.test(pred_traj_aux[,(n_time*2)],D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2)])$statistic
  
  
  #Calculate PI
  upp_rho_aux <- quantile(pred_traj_aux[,(n_time*2-1)],0.975)
  low_rho_aux <- quantile(pred_traj_aux[,(n_time*2-1)],0.025)
  width_rho_aux <- upp_rho_aux - low_rho_aux
  PI_rho_aux <- length(which((D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2-1)]<upp_rho_aux) & (D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2-1)] > low_rho_aux)))/sample_mrep
  
  upp_I_aux <- quantile(pred_traj_aux[,(n_time*2)],0.975)
  low_I_aux <- quantile(pred_traj_aux[,(n_time*2)],0.025)
  width_I_aux <- upp_I_aux - low_I_aux
  PI_I_aux <- length(which((D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2)]<upp_I_aux) & (D_star[((i-1)*sample_mrep+1):(i*sample_mrep),(n_time*2)] > low_I_aux)))/sample_mrep
  
  
  write.table(t(c(ks_rho, ks_rho_aux, ks_I, ks_I_aux, width_rho, width_rho_aux, 
                  width_I, width_I_aux, PI_rho, PI_rho_aux, PI_I, PI_I_aux)),"ks_total.txt",sep=",",append=TRUE,eol="\n",
              row.names=FALSE,col.names=FALSE);
}





