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

# mrep <- 1
args <- commandArgs(trailingOnly = TRUE)
mrep <- as.integer(args[1]) 
print(mrep)
R <- 3
setwd(getwd())
set.seed(mrep*1000+R*10+1)
data_whole = read.csv("simulation_data_noise_10.csv", header=F)
data_whole <- t(data_whole[,-1])
rownames(data_whole) <- NULL
n_time <- length(data_whole[,1])
ttname = data_whole[,1]
data <- matrix(0,n_time,2*R)
data[,seq(1,2*R,2)] <- data_whole[,1] 
data[,seq(2,2*R,2)] <- data_whole[,seq(mrep*R-(R-2),mrep*R+1,1)]

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


n_factor = 2
# n_time   = length(rho)
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


# celltherapy <- function(t, state, parameters) {
#   with(as.list(c(state, parameters)),{
#     drho = r_g * rho * (1 - (1 + exp(k_s * (k_c-I)))^-1 )
#     dI = drho -r_d * I
#     list(c(drho, dI)) # return the rate of change
#   })
# }
# 
# 
# process_simulation <- function(t,t_adjust,parameters,initial.state,exchange=0,cell_adjust=0){
#   loc_adjust = which(tname == t_adjust)
#   loc_end = which(tname == t)
#   out1 <- ode(y = initial.state, times = tname[1:loc_adjust], func = celltherapy, parms = parameters, method = 'vode')
#   
#   if (t != t_adjust) {
#     initial.state2 = out1[length(out1[,1]),]
#     initial.state2[3] = initial.state2[3] * (1 - exchange)
#     if (cell_adjust != 0) {
#       initial.state2[2] = cell_adjust
#     }
#     out2 <- ode(y = initial.state2[2:3], times = tname[(loc_adjust):loc_end], func = celltherapy, parms = parameters, method = 'vode')
#     out2 = out2[-1,]
#     out1 = rbind(out1,out2)
#   }
#   return(out1)
# }


# Number of posterior sample 
# B <- 50

# Number of particles at each ABC
N_eps <- 400

# Number of neighbours for covariance matrix calculations
M <- N_eps - 50

# number of time step
n_time_real <- 2


# Empty space to store results
res.old<-array(0,c(N_eps,6*(n_time_real-1)))
res.new<-array(0,c(N_eps,6*(n_time_real-1)))

# Empty space to store weights
w.old<-array(0,c(N_eps,1))
w.new<-array(0,c(N_eps,1))

# Empty space to store distance
dis.old<-array(0,c(N_eps,1))
dis.new<-array(0,c(N_eps,1))

# # Empty space to store distance threshold epsilon
# eps <- array(0,n_time)


# The prior of parameters 
para1 <- c(0,0,0,0,0,0)  
para2 <- c(0.5,5,5,0.05,0.2,0.2)  

# Number of simulations for each parameter set
n <- 60*R

# Number of generation
# G <- 7

# Epsilon values
# epsilon <- c(2,1,0.5,0.2,0.1,0.05,0.01)


#number of total batch
# R <- 6


# ############Test the Gibbs_ABC
# bio_prod_test <- matrix(0, R, n_time*n_factor)
# r_g_test <- rep(0.057,n_time)
# k_c_test <- rep(3.4,n_time)
# k_s_test <- rep(2.6,n_time)
# r_d_test <- rep(0.01,n_time)
# sd_rho_test <- rep(0.5,n_time)
# sd_I_test <- rep(0.5,n_time)
# 
# 
# for (r in 1:R) {
#   state <- as.numeric(dat[1,1:2])
#   bio_prod_test[r,1:2] <- state
#   for (t in 1:(n_time-1)) {
#     if (t %in% act_time) {
#       ####rho adjust
#       if (rho_act[which(act_time==t),r] != 0) {
#         state[1] = rho_act[which(act_time==t),r]
#       } else {
#         state[1] = state[1]
#       }
#       ### I adjust
#       state[2] = state[2] * as.numeric(1-I_act[which(act_time==t),r])
#     } else {
#       state.old = state
#       state[1] = rnorm(1, mean = state.old[1] + (tname[t+1] - tname[t]) *
#                          r_g_test[t] * state.old[1] * (1 - (1 + exp(k_s_test[t] * (k_c_test[t]-state.old[2])))^-1),
#                        sd = sd_rho_test[t])
#       state[2] = rnorm(1, mean = state.old[2] + (tname[t+1] - tname[t]) *
#                          (r_g_test[t] * state.old[1] * (1 - (1 + exp(k_s_test[t] * (k_c_test[t]-state.old[2])))^-1) - r_d_test[t] * state[2]),
#                        sd = sd_I_test[t])
#     }
#     bio_prod_test[r,(2*t+1):(2*(t+1))] = state
#   }
#   
# }
# 
# dat_Gibbs <- matrix(0, R, n_time*n_factor)
# 
# for (i in 1:n_time) {
#   dat_Gibbs[,(2*i-1)] <- bio_prod_test[,(2*i-1)]
# }

train_index <- seq(1:R)


dat_Gibbs = dat
n_time <- length(dat_Gibbs)/2



ABC_SMC <- function(Time, g, train_index, res.old, w.old, Sigma, eps,i){
  flag <- TRUE
  while (flag) {
    if (g==1){
      # Sample from prior distributions 
      # mu_g <- runif(Time-1, min=para1[1], max = para2[1])
      # std_g <- runif(Time-1, min=para1[2], max = para2[2])
      set.seed(g*1000+mrep*100+pacc_min*100+i*111)
      r_g <- runif(Time-1, min=para1[1], max = para2[1])
      k_s <- runif(Time-1, min=para1[2], max = para2[2])
      k_c <- runif(Time-1, min=para1[3], max = para2[3])
      r_d <- runif(Time-1, min=para1[4], max = para2[4])
      sd_rho <- runif(Time-1, min = para1[5], max = para2[5])
      sd_I <- runif(Time-1, min = para1[6], max = para2[6])
    } else{
      repeat {
        p<-sample(seq(1,N_eps),1,prob=w.old/sum(w.old))		
        sigma<-Sigma[,,p]
        # diag(sigma)[diag(sigma)<0] <- 0.001
        sigma <- sigma + diag(ncol(sigma))*0.0001
        par<- rK(res.old[p,],sigma)
        # print(x)
        # print(par)
        # x = x+1
        if (!(NaN %in% par)){
          break
        }
      }
      r_g<-par[1:(Time-1)]
      # mu_g <- par[1:(Time-1)]
      # std_g <- par[(1*(Time-1)+1):(2*(Time-1))]
      k_s<-par[(1*(Time-1)+1):(2*(Time-1))]
      k_c<-par[(2*(Time-1)+1):(3*(Time-1))]
      r_d<-par[(3*(Time-1)+1):(4*(Time-1))]
      sd_rho<-par[(4*(Time-1)+1):(5*(Time-1))]
      sd_I<-par[(5*(Time-1)+1):(6*(Time-1))]
    }
    if(prior.non.zero(rbind(r_g,k_s,k_c,r_d,sd_rho,sd_I))){
      # Set number of accepted simulations to zero
      m<-0
      # distance <-matrix(0, ncol=n,nrow=R)
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
      r_g_aug <- c(rep(r_g[1],10))
      k_s_aug <- c(rep(k_s[1],10))
      k_c_aug <- c(rep(k_c[1],10))
      r_d_aug <- c(rep(r_d[1],10))
      sd_rho_aug <- c(rep(sd_rho[1],10))
      sd_I_aug <- c(rep(sd_I[1],10))
      
      # calculate parameter based on observation
      coff <- rep(0,2*(n_time-1))
      for (t in 1:(n_time-1)) {
        # dat_Gibbs_cent <- t(t(dat_Gibbs) - apply(dat_Gibbs,2,mean))
        dat_Gibbs_cent <- dat_Gibbs
        x <- dat_Gibbs_cent[1:R,(2*t-1)]
        y <- dat_Gibbs_cent[1:R,(2*t+1)]
        temp_data <- data.frame(x,y)
        # dat_Gibbs_cent <- t(t(dat_Gibbs) - apply(dat_Gibbs,2,mean))
        temp_lm <- lm(y ~ 0 + x, data = temp_data)
        coff[2*t-1] <- temp_lm$coefficients[1]
        coff[2*t] <- sum((temp_lm$residuals)**2)/R
        # coff[1] <- 0
        # print(summary(temp_lm)) 
          # dat_Gibbs - rep(apply(dat_Gibbs,2,mean),3)
        # temp_lm <- lm(y ~ 0 + x, data = temp_data)
        # coff[2*t-1] <- temp_lm$coefficients[1]
        # coff[2*t] <- sum((temp_lm$residuals)**2)/R
        # print(summary(temp_lm))
      }
      
      D_star <- matrix(0, ncol=2*n_time,nrow=n)
      for (j in 1:n) {
        # r_g <- rnorm(2,mean = mu_g,sd = std_g)
        # r_g_aug <- c(rep(r_g[1],6),rep(r_g[2],4))
        D_star[j,]<-run_model_t(r_g_aug, k_s_aug, k_c_aug, r_d_aug, sd_rho_aug, sd_I_aug, 1, n_time)
        # Calculate distances 
      }
      
      # dat_Gibbs[1:R,seq(1,2*n_time,2)]
      coff_aux <- rep(0,2*(n_time-1))
      for (t in 1:(n_time-1)) {
        # D_star_cent <- t(t(D_star) - apply(D_star,2,mean))
        D_star_cent <- D_star
        x <- D_star_cent[1:R,(2*t-1)]
        y <- D_star_cent[1:R,(2*t+1)]
        temp_data_aux <- data.frame(x,y)
        temp_lm_aux <- lm(y ~ 0 + x, data = temp_data_aux)
        coff_aux[2*t-1] <- temp_lm_aux$coefficients[1]
        coff_aux[2*t] <- sum((temp_lm_aux$residuals)**2)/R
        # coff_aux[1] <- 0
        # print(summary(temp_lm))
      }
      
      distance <- Norm.Eucl.dist(coff_aux[seq(1,2*(n_time-1),2)],coff[seq(1,2*(n_time-1),2)])
      # Norm.Eucl.dist(coff_aux,coff)
      if (distance < eps) {
        m <- m+1
      }      
      # Store results
      res <- c(r_g,k_s, k_c, r_d, sd_rho, sd_I)  
      # Calculate weights
      res.temp <- t(matrix(res,(Time-1),6))
      prod_t <- rep(0,Time-1)
      for (t in 1:(Time-1)) {
        prod_t[t] <- prod(sapply(1:dim(res.temp)[1], function(b) dunif(res.temp[b,t], min = para1[b], max = para2[b])))
      }
      w1 <- prod(prod_t)
      
      # w1<-prod(sapply(1:6, function(b) dunif(res.new[i,b,t], min=para1[b], max=para2[b])))
      if(g==1){
        w2<-1
      } else {
        w2<-sum(sapply(1:(N_eps*alpha), function(a) w.old[a] * 
                         # for (o in 1:18) {
                         #   print(o)
                         #   print(dtmvnorm(res.new[i,1:o], mean=res.old[a,1:o], sigma=sigma[1:o,1:o], lower=rep(para1,each=Time-1)[1:o], upper=rep(para2,each=Time-1)[1:o]))
                         # }
                         dtmvnorm(res, mean=res.old[a,], sigma=sigma, lower=rep(para1,each=Time-1), upper=rep(para2,each=Time-1))
        ))
      }
      # res.new[i] <- res
      # w.new[i] <- (m/n)*w1/w2
      # dis.new[i] <- mean(colMeans(distance))
      # Update counter
      # print(paste0('Timestep:', t, ' Generation: ', g, ", particle: ", i))
      flag <- FALSE 
    } 
  }
  return(c(res,w1/w2,distance))
}


n_para <- 6 *(n_time_real-1)
eps <- 4000000
alpha <- 0.5
pacc_min <- 0.15
pacc <- 1


g <- 1
num_cores <- detectCores()
print(num_cores)
cl <- makeCluster(num_cores)
# cl <- makeCluster(56)
registerDoParallel(cl)
start_time <- Sys.time()
results <- foreach(i=1:N_eps, .combine="rbind") %dopar% {
  library(tmvtnorm)
  ABC_SMC(n_time_real, g, train_index, res.old, w.old, Sigma, eps,i)
}
time_cost <-  difftime(time1 = Sys.time(), time2 = start_time, units = "secs")
stopCluster(cl)
rownames(results) <- NULL
res.new <- results[,1:n_para]
w.new[,1] <- results[,n_para+1]
dis.new[,1] <- results[,n_para+2]
# time_cost <- results[,n_para+3]

temp<-data.frame(no=seq(1,N_eps), dist = dis.new)
dis.order<-temp[order(temp$dist),]$no

res.old<-res.new[dis.order,]
w.old<-w.new[dis.order,]
# w.old <- w.old/sum(w.old)
dis.old <-  dis.new[dis.order,]

Sigma <- array(NA, dim=c(n_para,n_para,N_eps))
for(p in 1:N_eps){
  Sigma[,,p]<- getSigmaNeighbours(M, res.old[p,], res.old[1:N_eps,]) 
}

eps <-  dis.old[alpha*N_eps]
print(eps)
w.old[(alpha*N_eps+1):N_eps]  <- 0
w.old <- w.old/sum(w.old)

res.new <- res.old
w.new <- w.old
dis.new <- dis.old

write.table(res.old,paste0("results_para_aux",mrep,".txt"),sep=",",append=TRUE,eol="\n",
            row.names=FALSE,col.names=FALSE);
write.table(w.old,paste0("results_weight_aux",mrep,".txt"),sep=",",append=TRUE,eol="\n",
            row.names=FALSE,col.names=FALSE);
write.table(dis.old,paste0("results_distance_aux",mrep,".txt"),sep=",",append=TRUE,eol="\n",
            row.names=FALSE,col.names=FALSE);
write.table(time_cost,paste0("results_time_aux",mrep,".txt"),sep=",",append=TRUE,eol="\n",
            row.names=FALSE,col.names=FALSE);

while (pacc > pacc_min) {
  g <- g+1
  num_cores <- detectCores()
  print(num_cores)
  cl <- makeCluster(num_cores)
  # cl <- makeCluster(56)
  registerDoParallel(cl)
  start_time <- Sys.time()
  results <- foreach(i=(1+alpha*N_eps):N_eps, .combine="rbind") %dopar% {
    library(tmvtnorm)
    ABC_SMC(n_time_real, g, train_index, res.old, w.old, Sigma, eps,i)
  }
  time_cost <-  difftime(time1 = Sys.time(), time2 = start_time, units = "secs")
  stopCluster(cl)
  rownames(results) <- NULL
  res.new[(1+alpha*N_eps):N_eps,] <- results[,1:n_para]
  w.new[(1+alpha*N_eps):N_eps] <- results[,n_para+1]
  dis.new[(1+alpha*N_eps):N_eps] <- results[,n_para+2]
  
  pacc <- (length(which(dis.new <= eps)) - alpha*N_eps)/(N_eps - alpha*N_eps)
  print(pacc)
  temp<-data.frame(no=seq(1,N_eps), dist = dis.new)
  dis.order<-temp[order(temp$dist),]$no
  
  res.old<-res.new[dis.order,]
  w.old<-w.new[dis.order]
  # w.old <- w.old/sum(w.old)
  dis.old <-  dis.new[dis.order]
  
  Sigma <- array(NA, dim=c(n_para,n_para,N_eps))
  for(p in 1:N_eps){
    Sigma[,,p]<- getSigmaNeighbours(M, res.old[p,], res.old[1:N_eps,]) 
  }
  
  eps <-  dis.old[alpha*N_eps]
  print(eps)
  w.old[(alpha*N_eps+1):N_eps]  <- 0
  w.old <- w.old/sum(w.old)
  
  res.new <- res.old
  w.new <- w.old
  dis.new <- dis.old
  
  write.table(res.old,paste0("results_para_aux",mrep,".txt"),sep=",",append=TRUE,eol="\n",
              row.names=FALSE,col.names=FALSE);
  write.table(w.old,paste0("results_weight_aux",mrep,".txt"),sep=",",append=TRUE,eol="\n",
              row.names=FALSE,col.names=FALSE);
  write.table(dis.old,paste0("results_distance_aux",mrep,".txt"),sep=",",append=TRUE,eol="\n",
              row.names=FALSE,col.names=FALSE);
  write.table(time_cost,paste0("results_time_aux",mrep,".txt"),sep=",",append=TRUE,eol="\n",
              row.names=FALSE,col.names=FALSE);
}  
    
    





    

# num_cores <- detectCores()
# cl <- makeCluster(4)
# registerDoParallel(cl)
# 
# 
# start.time <- proc.time()
# results <- foreach(i=1:100, .combine="rbind") %dopar% {
#   ABC_SMC(Time, g, epsilon, train_index, res.old, w.old, Sigma)
# }
# 
# 
# 
# stop.time <- proc.time()
# run.time <- stop.time - start.time
# print(run.time)
# stopCluster(cl)
#2 5.57
#3 4.08
#4 3.58

# start.time <- proc.time()
# for (i in 1:N_eps) {
#   results <- ABC_SMC(Time, g, epsilon, train_index, res.old, w.old, Sigma)
#   # print(i)
#   res.new[i,] <- results[[1]]
#   w.new[i] <- results[[2]]
#   dis.new[i] <- results[[3]]
# }
# stop.time <- proc.time()
# run.time <- stop.time - start.time
# print(run.time)
#10.34


# results <- ABC_SMC(n_time,N_eps,G,epsilon,train_index)



# results[1,2]




#  ABC_SMC <- function(Time, g, epsilon, train_index, res.old, w.old, Sigma){
#   flag <- TRUE
#   while (flag) {
#     if (g==1){
#       # Sample from prior distributions 
#       r_g <- runif(Time-1, min=para1[1], max = para2[1])
#       k_s <- runif(Time-1, min=para1[2], max = para2[2])
#       k_c <- runif(Time-1, min=para1[3], max = para2[3])
#       r_d <- runif(Time-1, min=para1[4], max = para2[4])
#       sd_rho <- runif(Time-1, min = para1[5], max = para2[5])
#       sd_I <- runif(Time-1, min = para1[6], max = para2[6])
#     } else{
#       p<-sample(seq(1,N_eps),1,prob=w.old/sum(w.old))		
#       sigma<-Sigma[,,p]
#       par<- rK(res.old[p,],sigma)
#       r_g<-par[1:(Time-1)]
#       k_s<-par[(1*(Time-1)+1):(2*(Time-1))]
#       k_c<-par[(2*(Time-1)+1):(3*(Time-1))]
#       r_d<-par[(3*(Time-1)+1):(4*(Time-1))]
#       sd_rho<-par[(4*(Time-1)+1):(5*(Time-1))]
#       sd_I<-par[(5*(Time-1)+1):(6*(Time-1))]
#     }
#     if(prior.non.zero(rbind(r_g,k_s,k_c,r_d,sd_rho,sd_I))){
#       # Set number of accepted simulations to zero
#       m<-0
#       distance <-matrix(0, ncol=n,nrow=R)
#       for (j in 1:n) {
#         for(r in train_index){
#           D_star<-run_model_t(r_g, k_s, k_c, r_d, sd_rho, sd_I, r, n_time)     
#           # Calculate distances 
#           distance[r,j] <- Norm.Eucl.dist(D_star[,seq(1,2*n_time,2)],dat_Gibbs[r,seq(1,2*n_time,2)])
#         }
#         if (mean(distance[,j]) < epsilon[g]) {
#           m <- m+1
#         }
#       }
#       
#       if (m>0){
#         # Store results
#         res <- c(r_g, k_s, k_c, r_d, sd_rho, sd_I)  
#         # Calculate weights
#         res.temp <- t(matrix(res,(Time-1),6))
#         prod_t <- rep(0,n_time-1)
#         for (t in 1:(n_time-1)) {
#           prod_t[t] <- prod(sapply(1:dim(res.temp)[1], function(b) dunif(res.temp[b,t], min = para1[b], max = para2[b])))
#         }
#         w1 <- prod(prod_t)
#         
#         # w1<-prod(sapply(1:6, function(b) dunif(res.new[i,b,t], min=para1[b], max=para2[b])))
#         if(g==1){
#           w2<-1
#         } else {
#           w2<-sum(sapply(1:N_eps, function(a) w.old[a] * 
#                            # for (o in 1:18) {
#                            #   print(o)
#                            #   print(dtmvnorm(res.new[i,1:o], mean=res.old[a,1:o], sigma=sigma[1:o,1:o], lower=rep(para1,each=Time-1)[1:o], upper=rep(para2,each=Time-1)[1:o]))
#                            # }
#                            dtmvnorm(res, mean=res.old[a,], sigma=sigma, lower=rep(para1,each=Time-1), upper=rep(para2,each=Time-1))
#           ))
#         }
#         # res.new[i] <- res
#         # w.new[i] <- (m/n)*w1/w2
#         # dis.new[i] <- mean(colMeans(distance))
#         # Update counter
#         # print(paste0('Timestep:', t, ' Generation: ', g, ", particle: ", i))
#         flag <- FALSE 
#       }
#     } 
#   } 
#   return(c(res,(m/n)*w1/w2,mean(colMeans(distance))))
# }


#   # Sigma <- list(NA, N_eps)
# Sigma <- array(NA, dim=c(n_para,n_para,N_eps))
#   for(p in 1:N_eps){
#     Sigma[,,p]<- getSigmaNeighbours(M, res.new[p,], res.new[1:N_eps,]) 
#   }
#   
#   temp<-data.frame(no=seq(1,N_eps), dist = dis.new)
#   dis.order<-temp[order(temp$dist),]$no
#   
#   res.old<-res.new[dis.order,]
#   w.old<-w.new[dis.order,]
#   dis.old <-  dis.new[dis.order,]
# }


# dat_Gibbs = dat[,1:8]
# n_time <- length(dat_Gibbs)/2
# 
# 
# ############Test the Gibbs_ABC
#  bio_prod_test <- matrix(0, R, n_time*n_factor)
#  r_g_test <- rep(0.057,n_time)
#  k_c_test <- rep(3.4,n_time)
#  k_s_test <- rep(2.6,n_time)
#  r_d_test <- rep(0.01,n_time)
#  sd_rho_test <- rep(0.00001,n_time)
#  sd_I_test <- rep(0.00001,n_time)
# 
# 
#  for (r in 1:R) {
#    state <- as.numeric(dat_Gibbs[r,1:2])
#    bio_prod_test[r,1:2] <- state
#    for (t in 1:(n_time-1)) {
#      if (t %in% act_time) {
#        ####rho adjust
#        if (rho_act[which(act_time==t),r] != 0) {
#          state[1] = rho_act[which(act_time==t),r]
#        } else {
#          state[1] = state[1]
#        }
#        ### I adjust
#        state[2] = state[2] * as.numeric(1-I_act[which(act_time==t),r])
#      } else {
#        state.old = state
#        state[1] = rnorm(1, mean = state.old[1] + (tname[t+1] - tname[t]) *
#                           r_g_test[t] * state.old[1] * (1 - (1 + exp(k_s_test[t] * (k_c_test[t]-state.old[2])))^-1),
#                         sd = sd_rho_test[t])
#        state[2] = rnorm(1, mean = state.old[2] + (tname[t+1] - tname[t]) *
#                           (r_g_test[t] * state.old[1] * (1 - (1 + exp(k_s_test[t] * (k_c_test[t]-state.old[2])))^-1) - r_d_test[t] * state[2]),
#                         sd = sd_I_test[t])
#      }
#      bio_prod_test[r,(2*t+1):(2*(t+1))] = state
#    }
# 
#  }
# 
#  for (i in 1:n_time) {
#    dat_Gibbs[,(2*i-1)] <- bio_prod_test[,(2*i-1)]
#  }
# 
#   
# 
# # Number of posterior sample 
# B <- 50
# 
# # Number of neighbours for covariance matrix calculations
# M <- 20
# 
# # Number of particles at each ABC
# N_eps <- 1000
# 
# # Empty space to store results
# res.old<-array(0,c(N_eps,6*(n_time-1)))
# res.new<-array(0,c(N_eps,6*(n_time-1)))
# 
# # Empty space to store weights
# w.old<-array(0,c(N_eps,1))
# w.new<-array(0,c(N_eps,1))
# 
# # Empty space to store distance
# dis.old<-array(0,c(N_eps,1))
# dis.new<-array(0,c(N_eps,1))
# 
# # # Empty space to store distance threshold epsilon
# # eps <- array(0,n_time)
# 
# 
# # The prior of parameters 
# para1 <- c(0,0,0,0,0,0)  
# para2 <- c(0.3,5,5,0.3,0.1,0.1)  
# 
# # Number of simulations for each parameter set
# n <- 10
# 
# # Number of generation
# G <- 1
# 
# # Epsilon values
# epsilon <- c(32,24,8,6,4,2,0.5,0.2,0.1)
# 
# # eps <- 0
# 
# # Parameters of Gibbs
# warmup <-40
# hsize <- 2
# 
# Total <- warmup + hsize*B
# 
# 
# r_g_s <- matrix(0, Total, n_time)
# k_s_s <- matrix(0, Total, n_time)
# k_c_s <- matrix(0, Total, n_time)
# r_d_s <- matrix(0, Total, n_time)
# sd_rho_s <- matrix(0, Total, n_time)
# sd_I_s <- matrix(0, Total, n_time)
# 
# 
# ########################################################
# # simulate data by GIBBS-ABC
# ########################################################
# 
# # Gibbs_ABC <- function(Total,train_index){}
# 
# train_index <- seq(1:11)
# 
# 
# library(doParallel)
# 
# num_cores <- detectCores()
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# 
# start.time <- proc.time()
# results <- ABC_SMC(n_time,N_eps,G,epsilon,train_index)
# stop.time <- proc.time()
# run.time <- stop.time - start.time
# print(run.time)
# 
# stopCluster(cl)
# #122.84
# 
# write.table(results[[1]],"results.txt",sep=",",append=TRUE,eol="\n",
#             row.names=FALSE,col.names=FALSE);


# res.old[,,t] <- results[[1]]
# w.old[,,t] <- results[[2]]
# dis.old[,,t] <- results[[3]]
# 
# r_g_pr <- results[[1]][,(0*(Time-1)+1):(1*(Time-1))]
# k_s_pr <- results[[1]][,(1*(Time-1)+1):(2*(Time-1))]
# k_c_pr <- results[[1]][,(2*(Time-1)+1):(3*(Time-1))]
# r_d_pr <- results[[1]][,(3*(Time-1)+1):(4*(Time-1))]
# sd_rho_pr <- results[[1]][,(4*(Time-1)+1):(5*(Time-1))]
# sd_I_pr <- results[[1]][,(5*(Time-1)+1):(6*(Time-1))]
# 
# 
# warr_g_pr <- runif(n_time, min=para1[1], max = para2[1])
# k_s_pr <- runif(n_time, min=para1[2], max = para2[2])
# k_c_pr <- runif(n_time, min=para1[3], max = para2[3])
# r_d_pr <- runif(n_time, min=para1[4], max = para2[4])
# sd_rho_pr <- runif(n_time, min = para1[5], max = para2[5])
# sd_I_pr <- runif(n_time, min = para1[6], max = para2[6])
# SIGMA <- array(0,dim=c(6,6,N_eps,n_time))
# 
# for (b in 1:Total){
#   
#   ##############################
#   # Generate sample path of Inhibitor I
#   # bio_process <- matrix(0, R, N)
#   for (r in train_index) {
#     dat_Gibbs[r,] <- run_model_trajectory(r_g_pr,k_s_pr,k_c_pr,r_d_pr,sd_rho_pr,sd_I_pr,r)
#   }
#   
#   for (t in 1:(n_time-1)){
#     ##############################
#     if (t %in% act_time) {
#       next
#     }
#     # ABC for theta
#     if (b == 1) {
#       results <- ABC_SMC(t,N_eps,G,epsilon,train_index)
#       res.old[,,t] <- results[[1]]
#       w.old[,,t] <- results[[2]]
#       dis.old[,,t] <- results[[3]]
#       # if (t == 1) {
#       #   SIGMA <- results[[4]]
#       #   A <- array(0,dim=c(6,6,10,2))
#       # } else{
#       #   SIGMA <- c(SIGMA,results[[4]])
#       # }
#       SIGMA[,,,t] <- results[[4]]
#       r_g_pr[t] <- results[[1]][1,1]
#       k_s_pr[t] <- results[[1]][1,2]
#       k_c_pr[t] <- results[[1]][1,3]
#       r_d_pr[t] <- results[[1]][1,4]
#       sd_rho_pr[t] <- results[[1]][1,5]
#       sd_I_pr[t] <- results[[1]][1,6]
#     } else {
#       results <- ABC_SMC(t,N_eps,G,epsilon,train_index,res.old,w.old,dis.old,SIGMA[,,,t],1)
#       res.old[,,t] <- results[[1]]
#       w.old[,,t] <- results[[2]]
#       dis.old[,,t] <- results[[3]]
#       SIGMA[,,,t] <- results[[4]]
#       r_g_pr[t] <- results[[1]][1,1]
#       k_s_pr[t] <- results[[1]][1,2]
#       k_c_pr[t] <- results[[1]][1,3]
#       r_d_pr[t] <- results[[1]][1,4]
#       sd_rho_pr[t] <- results[[1]][1,5]
#       sd_I_pr[t] <- results[[1]][1,6]
#     }
# 
#   }
#   
#   # record b-th sample
#   r_g_s[b,] <- r_g_pr
#   k_s_s[b,] <- k_s_pr
#   k_c_s[b,] <- k_c_pr
#   r_d_s[b,] <- r_d_pr
#   sd_rho_s[b,] <- sd_rho_pr
#   sd_I_s[b,] <- sd_I_pr
#   write.table(t(c(r_g_pr,k_s_pr,k_c_pr,r_d_pr,sd_rho_pr,sd_I_pr)),"results.txt",sep=",",append=TRUE,eol="\n",
#               row.names=FALSE,col.names=FALSE);
#   
#   if (b %in% c(warmup,warmup/2,warmup/4)) {
#     epsilon <- epsilon[-1]
#   }
#   
#   if ((b%%10)==0){
#     print(b)
#   }
# }
# 
# parameter_ps = read.csv("Parameter_Whole.csv", header=F)
# 
# # posterior index
# # id <- seq((warmup+1), Total, by=hsize)
# # 
# # r_g_ps <- r_g_s[id,]
# # k_c_ps <- k_c_s[id,]
# # k_s_ps <- k_s_s[id,]
# # r_d_ps <- r_d_s[id,]
# # sd_rho_ps <- sd_rho_s[id,]
# # sd_I_ps <- sd_I_s[id,]
# 
# id <- seq((20+1), Total, by=hsize)
# 
# r_g_ps <- parameter_ps[id,1:50]
# k_c_ps <- parameter_ps[id,51:100]
# k_s_ps <- parameter_ps[id,101:150]
# r_d_ps <- parameter_ps[id,151:200]
# sd_rho_ps <- parameter_ps[id,201:250]
# sd_I_ps <- parameter_ps[id,251:300]
# 
# ### Prediction #####
# r <- 9
# bio_prod <- matrix(0, B, n_time*n_factor)
# 
# for (b in 1:B) {
#   state <- as.numeric(dat_Gibbs[r,1:2])
#   bio_prod[b,1:2] <- state
#   for (t in 1:(n_time-1)) {
#     if (t %in% act_time) {
#       ####rho adjust
#       if (rho_act[which(act_time==t),r] != 0) {
#         state[1] = rho_act[which(act_time==t),r]
#       } else {
#         state[1] = state[1]
#       }
#       ### I adjust
#       state[2] = state[2] * as.numeric(1-I_act[which(act_time==t),r])
#     } else {
#       state.old = state
#       state[1] = rnorm(1, mean = state.old[1] + (tname[t+1] - tname[t]) *
#                          r_g_ps[b,t] * state.old[1] * (1 - (1 + exp(k_s_ps[b,t] * (k_c_ps[b,t]-state.old[2])))^-1),
#                        sd = sd_rho_ps[b,t])
#       state[2] = rnorm(1, mean = state.old[2] + (tname[t+1] - tname[t]) * 
#                          (r_g_ps[b,t] * state.old[1] * (1 - (1 + exp(k_s_ps[b,t] * (k_c_ps[b,t]-state.old[2])))^-1) - r_d_ps[b,t] * state[2]), 
#                        sd = sd_I_ps[b,t]) 
#     }
#     bio_prod[b,(2*t+1):(2*(t+1))] = state
#   }
# }
# 
# dat_prod_rho <- bio_prod[,seq(1,ncol(bio_prod),2)]
# dat_prod_I <- bio_prod[,seq(0,ncol(bio_prod),2)]
# 
# plot(tname[1:n_time], dat[r,seq(1,ncol(dat),2)][1:n_time],ylim = c(0,11), main = "D1", xlab = "Time (hrs)", ylab = "Cells/ml (*1E6)",col = 'red')
# dat_prod_rho_95_up = apply(dat_prod_rho,2,quantile, probs = 0.975)
# dat_prod_rho_95_low = apply(dat_prod_rho,2,quantile, probs = 0.025)
# lines(tname[1:n_time],dat_prod_rho_95_up,col = "blue")
# lines(tname[1:n_time],dat_prod_rho_95_low,col = "blue")
# legend(1, 11, legend=c("raw data", "95% PI Upper Bound","95% PI Lower Bound"),
#        col=c("red", "blue","blue"), pch = c(1, NA, NA), lty=c(NA,1,1), cex=1.5)
# 
# 
# dat_prod_I_95_up = apply(dat_prod_I,2,quantile, probs = 0.975)
# dat_prod_I_95_low = apply(dat_prod_I,2,quantile, probs = 0.025)
# plot(tname[1:n_time], dat_prod_I_95_up, main = "D3", xlab = "Time (hrs)", ylab = "Inhibitor",col = 'blue','l')
# lines(tname[1:n_time],dat_prod_I_95_low,col = "blue")
# 
# 
# ###Validation######
# 
# r_g_test <- 0.0541
# k_s_test <- 2.6
# k_c_test <- 3.4
# r_d_test <- 0.014
# sd_rho_test <- 0.0001
# sd_I_test <- 0.0001
# 
# 
# t <- 30
# 
# par(mfrow=c(3,2))
# plot(x = 1:length(parameter_ps[,1]), y = parameter_ps[,t], xlab = 'steps', ylab = 'r_g',cex=3, ylim = c(0,0.5))
# abline(h = r_g_test,col = 'red')
# 
# plot(x = 1:length(parameter_ps[,1]), y = parameter_ps[,(t+50)], xlab = 'steps', ylab = 'k_s',cex=3, ylim = c(0,5))
# abline(h = k_s_test,col = 'red')
# 
# plot(x = 1:length(parameter_ps[,1]), y = parameter_ps[,(t+100)], xlab = 'steps', ylab = 'k_c',cex=3, ylim = c(0,5))
# abline(h = k_c_test,col = 'red')
# 
# plot(x = 1:length(parameter_ps[,1]), y = parameter_ps[,(t+150)], xlab = 'steps', ylab = 'r_d',cex=3, ylim = c(0,0.2))
# abline(h = r_d_test,col = 'red')
# 
# plot(x = 1:length(parameter_ps[,1]), y = parameter_ps[,(t+200)], xlab = 'steps', ylab = 'sd_rho',cex=3, ylim = c(0,0.5))
# abline(h = sd_rho_test,col = 'red')
# 
# plot(x = 1:length(parameter_ps[,1]), y = parameter_ps[,(t+250)], xlab = 'steps', ylab = 'sd_I',cex=3, ylim = c(0,0.5))
# abline(h = sd_I_test,col = 'red')
# par(mfrow=c(1,1))
# 
# 
