##### Code for Figure 5 #####

library(foreach)
library(doParallel)

source("algorithms.R")

### Experiments for Figure 5 ###
n = 500         # Dimension of signal vector x
num_rep = 100   # Number of Monte Carlo trials
num_try = 100   # Maximum number of allowed restarts
prec = 0.01     # Precision required to declare a run successful

# Experiment 1, k=20 fixed
k = 20            # Sparsity level of x
m = 100 * (1:10)  # Number of observations

# Parallelize computation
numcl = 10
cl <- makeCluster(numcl)
registerDoParallel(cl)
clusterExport(cl,list("gen_data_sparse", "hwf", "wf_loss", "wf_grad"))
clusterExport(cl,list("n","k", "m", "num_rep", "num_try", "prec"),envir=environment())

# Run HWF num_rep times with the specified parameters
success_k20 <- foreach(rep = 1:num_rep, .combine=rbind, .multicombine=TRUE) %dopar%{
  B = rep(Inf, length(m))
  for(i in 1:length(m)){
    # Generate data
    data = gen_data_sparse(m[i], n, k)
    
    # Signal reconstruction using HWF
    ind_sup = sort(t(data$A^2) %*% data$y, index.return = TRUE)$ix[(n-num_try+1):n]
    for(trial in 1:num_try){
      # To save computation, only run HWF when we pick a coordinate on the support
      if(data$x[ind_sup[num_try+1-trial]]!=0){
        res = hwf(data$A, data$y, ini = ind_sup[num_try+1-trial], iteration = 100000)
        # Check whether this run was successful, if so save for which trial
        if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
          B[i] = trial
          break
        }
      }
    } 
  }
  B
}
stopCluster(cl)

# Save result
saveRDS(success_k20, file = "output_numtrial_k20.rds")

# Experiment 2, m=500 fixed
k = 5 + 5*(1:23)  # Sparsity level of x
m = 500           # Number of observations

# Parallelize computation
numcl = 10
cl <- makeCluster(numcl)
registerDoParallel(cl)
clusterExport(cl,list("gen_data_sparse", "hwf", "wf_loss", "wf_grad"))
clusterExport(cl,list("n","k", "m", "num_rep", "num_try", "prec"),envir=environment())

# Run HWF num_rep times with the specified parameters
success_m500 <- foreach(rep = 1:num_rep, .combine=rbind, .multicombine=TRUE) %dopar%{
  B = rep(Inf, length(k))
  for(j in 1:length(k)){
    # Generate data
    data = gen_data_sparse(m, n, k[j])
    
    # Signal reconstruction using HWF
    ind_sup = sort(t(data$A^2) %*% data$y, index.return = TRUE)$ix[(n-num_try+1):n]
    for(trial in 1:num_try){
      # To save computation, only run HWF when we pick a coordinate on the support
      if(data$x[ind_sup[num_try+1-trial]]!=0){
        res = hwf(data$A, data$y, ini = ind_sup[num_try+1-trial], iteration = 100000)
        # Check whether this run was successful, if so save for which trial
        if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
          B[j] = trial
          break
        }
      }
    } 
  }
  B
}
stopCluster(cl)

# Save result
saveRDS(success_m500, file = "output_numtrial_m500.rds")


### Generate plot ###

# Convert number of trials needed to success/failure for given number of trials
plot_numtrial1k20 = success_k20
plot_numtrial5k20 = success_k20
plot_numtrial10k20 = success_k20
plot_numtrial50k20 = success_k20
plot_numtrial100k20 = success_k20
for(i in 1:dim(success_k20)[1]){
  for(j in 1:dim(success_k20)[2]){
    plot_numtrial1k20[i,j] = as.numeric(success_k20[i,j]<=1)
    plot_numtrial5k20[i,j] = as.numeric(success_k20[i,j]<=5)
    plot_numtrial10k20[i,j] = as.numeric(success_k20[i,j]<=10)
    plot_numtrial50k20[i,j] = as.numeric(success_k20[i,j]<=50)
    plot_numtrial100k20[i,j] = as.numeric(success_k20[i,j]<=100)
  }
}

plot_numtrial1m500 = success_m500
plot_numtrial5m500 = success_m500
plot_numtrial10m500 = success_m500
plot_numtrial50m500 = success_m500
plot_numtrial100m500 = success_m500
for(i in 1:dim(success_m500)[1]){
  for(j in 1:dim(success_m500)[2]){
    plot_numtrial1m500[i,j] = as.numeric(success_m500[i,j]<=1)
    plot_numtrial5m500[i,j] = as.numeric(success_m500[i,j]<=5)
    plot_numtrial10m500[i,j] = as.numeric(success_m500[i,j]<=10)
    plot_numtrial50m500[i,j] = as.numeric(success_m500[i,j]<=50)
    plot_numtrial100m500[i,j] = as.numeric(success_m500[i,j]<=100)
  }
}

pdf(file = "plot_numtrial.pdf", width = 8, height = 3)

par(mfrow=c(1,2), mai = c(0.7, 0.7, 0.1, 0.1), bg = "transparent")

# Left plot: k=20 fixed
plot(apply(plot_numtrial1k20, 2, mean), type = "o", lwd = 2, pch = 16, col = "red", ylim = c(0,1), ylab = "", xlab = "", xaxt = "n")
axis(1, at = (1:10), label = 100*(1:10), cex.axis = 1.2)
title(xlab = "Number of measurements m", ylab = "Success rate", line = 2.2, cex.lab = 1.4)
lines(apply(plot_numtrial5k20, 2, mean), type = "o", lwd = 2, pch = 18, col = "magenta")
lines(apply(plot_numtrial10k20, 2, mean), type = "o", lwd = 2, pch = 15, col = "blue")
lines(apply(plot_numtrial50k20, 2, mean), type = "o", lwd = 2, pch = 17, col = "black")
lines(apply(plot_numtrial100k20, 2, mean), type = "o", lwd = 2, pch = 4, col = "brown")
legend('bottomright',legend=c("B=1", "B=5", "B=10", "B=50", "B=100"), col=c("red", "magenta", "blue", "black", "brown"), pch = c(16,18,15,17,4), lwd = 2, cex = 0.9, inset = c(0.01, 0), bty = "n")

# Right plot: m=500 fixed
plot(apply(plot_numtrial1m500, 2, mean), type = "o", lwd = 2, pch = 16, col = "red", ylim = c(0,1), ylab = "", xlab = "", xaxt = "n")
axis(1, at = (1 + 2*(0:11)), label = 10*(1:12), cex.axis = 1.2)
title(xlab = "Sparsity level k", ylab = "Success rate", line = 2.2, cex.lab = 1.4)
lines(apply(plot_numtrial5m500, 2, mean), type = "o", lwd = 2, pch = 18, col = "magenta")
lines(apply(plot_numtrial10m500, 2, mean), type = "o", lwd = 2, pch = 15, col = "blue")
lines(apply(plot_numtrial50m500, 2, mean), type = "o", lwd = 2, pch = 17, col = "black")
lines(apply(plot_numtrial100m500, 2, mean), type = "o", lwd = 2, pch = 4, col = "brown")
legend('topright',legend=c("B=1", "B=5", "B=10", "B=50", "B=100"), col=c("red", "magenta", "blue", "black", "brown"), pch = c(16,18,15,17,4), lwd = 2, cex = 0.7, inset = c(0.01, 0), bty = "n")

dev.off()