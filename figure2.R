##### Code for Figure 2 ##### 

library(foreach)
library(doParallel)
library(abind)

source("algorithms.R")

acomb <- function(...) abind(..., along=2)

### Experiments for Figure 2 ###
n = 500         # Dimension of signal vector x
num_rep = 100   # Number of Monte Carlo trials  
num_try = 50    # Maximum number of allowed restarts
prec = 0.01     # Precision required to declare a run successful

# Experiment 1, k=20 fixed
k = 20            # Sparsity level of x
m = 100 * (1:10)  # Number of observations

# Parallelize computation
numcl = 10
cl <- makeCluster(numcl)
registerDoParallel(cl)
clusterExport(cl,list("gen_data_sparse", "hwf", "alg4", "sparta", "swf", "wf_loss", "wf_grad", "rwf_grad"))
clusterExport(cl,list("n","k", "m", "num_rep", "num_try", "prec"),envir=environment())

# Run different reconstruction algorithms num_rep times with the specified parameters
success_k20 <- foreach(rep = 1:num_rep, .combine='acomb', .multicombine=TRUE) %dopar%{
  B = array(0, dim = c(length(m), 1, 4))
  for(i in 1:length(m)){
    # Generate data
    data = gen_data_sparse(m[i], n, k)
    
    # Signal reconstruction using HWF
    ind_sup = sort(t(data$A^2) %*% data$y, index.return = TRUE)$ix[(n-num_try+1):n]
    for(trial in 1:num_try){
      # To save computation, only run HWF when we pick a coordinate on the support
      if(data$x[ind_sup[num_try+1-trial]]!=0){
        res = hwf(data$A, data$y, ini = ind_sup[num_try+1-trial], iteration = 100000)
        # Check whether this run was successful
        if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
          B[i, 1, 1] = 1
          break
        }
      }
    } 
    
    # Signal reconstruction using Algorithm 4
    for(trial in 1:num_try){
      res = alg4(data$A, data$y, k, trial = trial, iteration = 10000)
      # To save computation on unsuccessful runs, we run a maximum of 10000 iterations;
      # there is no benefit in running the algorithm for more iterations
      
      # Check whether this run was successful
      if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
        B[i, 1, 2] = 1
        break
      }
    } 
    
    # Signal reconstruction using SPARTA
    res = sparta(data$A, data$y, k = k, iteration = 10000)
    # Check whether this run was successful
    if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
      B[i, 1, 3] = 1
    }
    
    # Signal reconstruction using SWF
    res = swf(data$A, data$y, k = k, iteration = 10000)
    # Check whether this run was successful
    if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
      B[i, 1, 4] = 1
    }
  }
  B
}
stopCluster(cl)

# Save result
saveRDS(success_k20, file = "output_k20.rds")

# Experiment 2, m=500 fixed
k = 5 + 5*(1:23)  # Sparsity level of x
m = 500           # Number of observations

# Parallelize computation
numcl = 10
cl <- makeCluster(numcl)
registerDoParallel(cl)
clusterExport(cl,list("gen_data_sparse", "hwf", "alg4", "sparta", "swf", "wf_loss", "wf_grad", "rwf_grad"))
clusterExport(cl,list("n","k", "m", "num_rep", "num_try"),envir=environment())

# Run different reconstruction algorithms num_rep times with the specified parameters
success_m500 <- foreach(rep = 1:num_rep, .combine='acomb', .multicombine=TRUE) %dopar%{
  B = array(0, dim = c(length(k), 1, 4))
  for(j in 1:length(k)){
    # Generate data
    data = gen_data_sparse(m, n, k[j])
    
    # Signal reconstruction using HWF
    ind_sup = sort(t(data$A^2) %*% data$y, index.return = TRUE)$ix[(n-num_try+1):n]
    for(trial in 1:num_try){
      # To save computation, only run HWF when we pick a coordinate on the support
      if(data$x[ind_sup[num_try+1-trial]]!=0){
        res = hwf(data$A, data$y, ini = ind_sup[num_try+1-trial], iteration = 100000)
        # Check whether this run was successful
        if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
          B[j, 1, 1] = 1
          break
        }
      }
    } 
    
    # Signal reconstruction using Algorithm 4
    for(trial in 1:num_try){
      res = alg4(data$A, data$y, k[j], trial = trial, iteration = 10000)
      # To save computation on unsuccessful runs, we run 10000 iterations; there is 
      # no benefit in running the algorithm for more iterations
      
      # Check whether this run was successful
      if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
        B[j, 1, 2] = 1
        break
      }
    } 
    
    # Signal reconstruction using SPARTA
    res = sparta(data$A, data$y, k = k[j], iteration = 10000)
    # Check whether this run was successful
    if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
      B[j, 1, 3] = 1
    }
    
    # Signal reconstruction using SWF
    res = swf(data$A, data$y, k = k[j], iteration = 10000)
    # Check whether this run was successful
    if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
      B[j, 1, 4] = 1
    }
  }
  B
}
stopCluster(cl)

saveRDS(success_m500, file = "output_m500.rds")

# Results for PR-GAMP obtained using the code available from https://sourceforge.net/projects/gampmatlab/
prgamp_k20 = c(0, 0.28, 0.79, 0.92, 1, 1, 1, 1, 1, 1)
prgamp_m500 = c(1, 1, 1, 0.96, 0.88, 0.85, 0.75, 0.71, 0.63, 0.53, 0.45, 0.31, 0.28, 0.14, 0.12)

### Generate plot ###

pdf(file = "plot_k20m500.pdf", width = 8, height = 3)

par(mfrow=c(1,2), mai = c(0.7, 0.7, 0.1, 0.1), bg = "transparent")

# Left plot: k=20 fixed
plot(prgamp_k20, type = "o", lwd = 2, pch = 4, col = "brown", ylim = c(0,1), ylab = "", xlab = "", xaxt = "n")
axis(1, at = (1:10), label = 100*(1:10))
title(xlab = "Number of measurements m", ylab = "Success rate", line = 2.5)
lines(apply(success_k20, c(1,3), mean)[,2], type = "o", lwd = 2, pch = 18, col = "magenta")
lines(apply(success_k20, c(1,3), mean)[,3], type = "o", lwd = 2, pch = 15, col = "blue")
lines(apply(success_k20, c(1,3), mean)[,4], type = "o", lwd = 2, pch = 17, col = "black")
lines(apply(success_k20, c(1,3), mean)[,1], type = "o", lwd = 2, pch = 16, col = "red")
legend('bottomright',legend=c("HWF", "SWF", "SPARTA", "PR-GAMP", "SPARTA-support"), col=c("red", "black", "blue", "brown", "magenta"), pch = c(16,17,15,4,18), lwd = 2, cex = 0.8, inset = c(0, -0.03), bty = "n")

# Right plot: m=500 fixed
plot(prgamp_m500, type = "o", lwd = 2, pch = 4, col = "brown", ylim = c(0,1), ylab = "", xlab = "", xaxt = "n")
axis(1, at = (1 + 2*(0:11)), label = 10*(1:12))
title(xlab = "Sparsity level k", ylab = "Success rate", line = 2.5)
lines(apply(success_m500, c(1,3), mean)[,2], type = "o", lwd = 2, pch = 18, col = "magenta")
lines(apply(success_m500, c(1,3), mean)[,3], type = "o", lwd = 2, pch = 15, col = "blue")
lines(apply(success_m500, c(1,3), mean)[,4], type = "o", lwd = 2, pch = 17, col = "black")
lines(apply(success_m500, c(1,3), mean)[,1], type = "o", lwd = 2, pch = 16, col = "red")
legend('bottomleft',legend=c("HWF", "SWF", "SPARTA", "PR-GAMP", "SPARTA-support"), col=c("red", "black", "blue", "brown", "magenta"), pch = c(16,17,15,4,18), lwd = 2, cex = 0.8, inset = c(0, -0.03), bty = "n")

dev.off()