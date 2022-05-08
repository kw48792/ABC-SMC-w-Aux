run_model<-function(r_g, k_s, k_c, r_d, sd_rho, sd_I, r, timestep){
  t = timestep
  state <- as.numeric(dat_Gibbs[r,(2*t-1):(2*t)])
  if (t %in% act_time) {
    ####rho adjust
    if (rho_act[which(act_time==t),r] != 0) {
      state[1] = rho_act[which(act_time==t),r]
      } else {
        state[1] = state[1]
      }
    #### I adjust
    state[2] = state[2] * as.numeric(1-I_act[which(act_time==t),r])
  } else {
    state.old = state
    state[1] = rnorm(1, mean = state.old[1] + (tname[t+1] - tname[t]) * 
              r_g * state.old[1] * (1 - (1 + exp(k_s * (k_c-state.old[2])))^-1), 
              sd = sd_rho) 
    state[2] = rnorm(1, mean = state.old[2] + (tname[t+1] - tname[t]) * 
              (r_g * state.old[1] * (1 - (1 + exp(k_s * (k_c-state.old[2])))^-1) - r_d * state[2]), 
              sd = sd_I) 
  }
  return(state)  
}

run_model_t<-function(r_g, k_s, k_c, r_d, sd_rho, sd_I, r, timestep){
  # state <- as.numeric(dat_Gibbs[r,1:2])
  state <- c(3,0)
  bio_prod_t <- matrix(0, 1, timestep*n_factor)
  bio_prod_t[1,1:2] <- state
  for (t in 1:(timestep-1)) {
    if (t %in% act_time) {
      ####rho adjust
      if (rho_act[which(act_time==t),r] != 0) {
        state[1] = rho_act[which(act_time==t),r]
      } else {
        state[1] = state[1]
      }
      #### I adjust
      state[2] = state[2] * as.numeric(1-I_act[which(act_time==t),r])
    } else {
      state.old = state
      state[1] = rnorm(1, mean = state.old[1] + (tname[t+1] - tname[t]) * 
                         r_g[t] * state.old[1] * (1 - (1 + exp(k_s[t] * (k_c[t]-state.old[2])))^-1), 
                       sd = sd_rho[t]) 
      state[2] = rnorm(1, mean = state.old[2] + (tname[t+1] - tname[t]) * 
                         (r_g[t] * state.old[1] * (1 - (1 + exp(k_s[t] * (k_c[t]-state.old[2])))^-1) - r_d[t] * state[2]), 
                       sd = sd_I[t]) 
    }
    bio_prod_t[1,(2*t+1):(2*(t+1))] = state
  }
  return(bio_prod_t)  
}


run_model_trajectory<-function(r_g, k_s, k_c, r_d, sd_rho, sd_I, r){
  bio_prod <- matrix(0, 1, N)
  state <- as.numeric(dat_Gibbs[r,1:2])
  bio_prod[,1:2] <- state
  for (t in 1:(n_time-1)) {
    if (t %in% act_time) {
      # ####rho adjust
      # if (rho_act[which(act_time==t),r] != 0) {
      #   state[1] = rho_act[which(act_time==t),r]
      # } else {
      #   state[1] = state[1]
      # }
      #### I adjust
      state[2] = state[2] * as.numeric(1-I_act[which(act_time==t),r])
    } else {
      state.old = state
      # state[1] = rnorm(1, mean = state.old[1] + (tname[t+1] - tname[t]) * 
      #                    r_g[t] * state.old[1] * (1 - (1 + exp(k_s[t] * (k_c[t]-state.old[2])))^-1), 
      #                  sd = sqrt(v2_rho[t])) 
      state[2] = rnorm(1, mean = state.old[2] + (tname[t+1] - tname[t]) * 
                         (r_g[t] * state.old[1] * (1 - (1 + exp(k_s[t] * (k_c[t]-state.old[2])))^-1) - r_d[t] * state[2]), 
                       sd = sd_I[t]) 
    }
    state[1] <- as.numeric(dat_Gibbs[r,(2*t+1)])
    bio_prod[,(2*t+1):(2*(t+1))] = state
  }
  return(bio_prod)  
}


# rK <- function(mean, sigma){
#   return(rmvnorm(1,mean=mean, sigma=sigma))
# }

rK <- function(mean, sigma){
  return(rtmvnorm(1,mean=mean, sigma=sigma, lower=rep(para1,each=length(mean)/6), upper=rep(para2,each=length(mean)/6),algorithm="gibbs", burn.in.samples = 200))
}

Norm.Eucl.dist<-function(p1,p2){
  sqrt(sum((p1-p2)^2)) }

#  Covariance based on M neighbours
getSigmaNeighbours<-function(M, theta, Theta){
  N <-  dim(Theta)[1]
  dist<- sapply(1:N, function(a) Norm.Eucl.dist(as.numeric(theta), as.numeric(Theta[a,])))
  temp<-data.frame(no=seq(1,N), dist)
  temp<-temp[order(temp$dist),]
  sigma<-cov(Theta[temp$no[1:(M+1)],])
  return(sigma)
}


# calc_distance <- function(pred,obs){
#   dist<- sapply(1:30, function(a) Norm.Eucl.dist(pred[a,seq(from=1, to=N, by=2)], obs[seq(from=1, to=N, by=2)]))
#   return(dist)
# }
# seq(from=1, to=N, by=2)
# Perturbation kernel 


#  Identity function: H(x)= 1 if x=T
H <- function(x) as.numeric(x>0)

#  Test if prior is non zero
prior.non.zero<-function(par){
    n_time <- dim(par)[2]
    prod_t <- rep(0,n_time)
    for (t in 1:n_time) {
      prod_t[t] <- prod(sapply(1:dim(par)[1], function(a) H(par[a,t]-para1[a])* H(para2[a]-par[a,t])))
    }
    return(prod(prod_t))
}




