
simulate <- function(path_parm, path_weight, total_particles, N){
  sample_size <- total_particles / 2 #Number of particles
  results_whole <- read.table(path_parm, head=FALSE, sep = ',')
  results_weight <- read.table(path_weight, head=FALSE, sep = ',')
  G <- dim(results_whole)[1] / total_particles
  output <- c() # data.frame(row.names=c('r_g', 'k_s', 'k_c', 'r_d', 'v_rho', 'v_I')) #matrix(0, ncol=6, nrow=N)
  results_g <- results_whole[(total_particles*(G-1)+1):(total_particles*(G-1)+sample_size),]
  results_weight_g <- results_weight[(total_particles*(G-1)+1):(total_particles*(G-1)+sample_size),]
  for (s in 1:N) {
    # print(results_weight_g)
    i <- sample.int(sample_size, 1, replace=TRUE, prob=results_weight_g)
    output<- rbind(output, results_g[i,])
  }
  return(output)
}

simulate_all <- function (base, method='_aux', total_particles=400, N=1000) {
  r_g <- c()
  k_s <- c()
  k_c <- c()
  r_d <- c()
  for (m in 1:30) {
    path_parm <- paste(base, 'results_para', method, m, '.txt', sep='')
    path_weight <- paste(base, 'results_weight', method, m, '.txt', sep='')
    print("===================")
    print(path_parm)
    print(path_weight)
    result <- simulate(path_parm, path_weight, total_particles, N)
    r_g <- cbind(r_g, result[,1])
    k_s <- cbind(k_s, result[,2])
    k_c <- cbind(k_c, result[,3])
    r_d <- cbind(r_d, result[,4])
  }
  colnames(r_g) <- apply(array(1:30, dim=c(1,30)), 1:2, function(x) paste("r", x, sep=''))
  colnames(k_s) <- apply(array(1:30, dim=c(1,30)), 1:2, function(x) paste("r", x, sep=''))
  colnames(k_c) <- apply(array(1:30, dim=c(1,30)), 1:2, function(x) paste("r", x, sep=''))
  colnames(r_d) <- apply(array(1:30, dim=c(1,30)), 1:2, function(x) paste("r", x, sep=''))
  return (list(r_g = r_g, k_s = k_s, k_c = k_c, r_d = r_d))
}

data_param <- c()
for (m in c(3,6,20)) {
  for (noise in c(0.1, 0.2)) {
    # m <- 3
    # noise <- 0.2
    base <- paste('data/m', m, 'n', sub("\\D+", "", noise, perl = TRUE), '/', sep='')
    result <- simulate_all(base, method='')
    r_g_df <- result$r_g
    k_s_df <- result$k_s
    k_c_df <- result$k_c
    r_d_df <- result$r_d

    r_g_df<- melt(r_g_df)
    r_g_df['m'] <- paste("m=", m, sep='')
    r_g_df['noise'] <- paste('v=', noise, sep='')
    r_g_df['parameter'] <- 'r_g'

    k_s_df<- melt(k_s_df)
    k_s_df['m'] <- paste("m=", m, sep='')
    k_s_df['noise'] <- paste('v=', noise, sep='')
    k_s_df['parameter'] <- 'k_s'

    k_c_df<- melt(k_c_df)
    k_c_df['m'] <- paste("m=", m, sep='')
    k_c_df['noise'] <- paste('v=', noise, sep='')
    k_c_df['parameter'] <- 'k_c'

    r_d_df<- melt(r_d_df)
    r_d_df['m'] <- paste("m=", m, sep='')
    r_d_df['noise'] <- paste('v=', noise, sep='')
    r_d_df['parameter'] <- 'r_d'

    data_param <- rbind(data_param, r_g_df, k_s_df, k_c_df, r_d_df)
  }
}
data_param_selected <- data_param[data_param$Var2 %in% c("r1", "r2", "r7", "r4", "r8", "r6"),]

mean_values <- c(0.057, 3.4, 2.6, 0.005)
params <- c('r_g', 'k_s', 'k_c', 'r_d')

para1 <- c(0,0,0,0,0,0)
para2 <- c(0.2,5,5,0.05,0.2,0.2)


k <- 4
i <- params[k]
j <- mean_values[k]
data_selected <- data_param_selected[data_param_selected['parameter'] == i, ]
g <-ggplot(data_selected, aes(x=value, fill=Var2)) + geom_density(color=NA, alpha=0.4) +
  geom_vline(xintercept=j, size=1.2, color="black", linetype=2)+
  ylab("Density") +
  # xlab(TeX("$k_s$")) +
  xlim(para1[k], para2[k]) +
  facet_grid(noise ~ factor(m, levels=c('m=3','m=6','m=20')) , scales="free")
g + scale_fill_discrete(name = "Macro Replication") +
  theme(legend.position = "bottomright",
        axis.text.x = element_text(size = 16, angle=40),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 20,face="bold"),
        axis.title.x = element_blank(),
        axis.title=element_text(size=20,face="bold"))



