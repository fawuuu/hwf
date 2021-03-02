##### Code for Figure 3 #####

library(foreach)
library(doParallel)
library(abind)
library(latex2exp)

source("algorithms.R")

acomb <- function(...) abind(..., along=3)

### Experiments for Figure 3 ###
n = 500           # Dimension of signal vector x
k = 5 + 5*(1:23)  # Sparsity level of x
m = 500           # Number of observations
num_rep = 100     # Number of Monte Carlo trials
num_try = 50      # Maximum number of allowed restarts
prec = 0.01       # Precision required to declare a run successful
  
# Parallelize computation
numcl = 10
cl <- makeCluster(numcl)
registerDoParallel(cl)
clusterExport(cl,list("gen_data_sparse", "hwf", "alg4", "sparta", "swf", "wf_loss", "wf_grad", "rwf_grad"))
clusterExport(cl,list("n","k", "m", "num_rep", "num_try", "prec"),envir=environment())

# Run different reconstruction algorithms num_rep times with the specified parameters
success_xmax <- foreach(rep = 1:num_rep, .combine='acomb', .multicombine=TRUE) %dopar%{
  B = array(0, dim = c(length(k), 3, 1, 4))
  for(j in 1:length(k)){
    for(i in 1:3){
      # Generate data
      data = gen_data_sparse(m, n, k[j], xmax = (i+1))
      
      # Signal reconstruction using HWF
      ind_sup = sort(t(data$A^2) %*% data$y, index.return = TRUE)$ix[(n-num_try+1):n]
      for(trial in 1:num_try){
        # To save computation, only run HWF when we pick a coordinate on the support
        if(data$x[ind_sup[num_try+1-trial]]!=0){
          res = hwf(data$A, data$y, ini = ind_sup[num_try+1-trial], iteration = 100000)
          # Check whether this run was successful
          if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
            B[j, i, 1, 1] = 1
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
          B[j, i, 1, 2] = 1
          break
        }
      } 
      
      # Signal reconstruction using SPARTA
      res = sparta(data$A, data$y, k = k[j], iteration = 10000)
      # Check whether this run was successful
      if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
        B[j, i, 1, 3] = 1
      }
      
      # Signal reconstruction using SWF
      res = swf(data$A, data$y, k = k[j], iteration = 10000)
      # Check whether this run was successful
      if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
        B[j, i, 1, 4] = 1
      }
    }
  }
  B
}
stopCluster(cl)

# Save result
saveRDS(success_xmax, file = "output_xmax.rds")

# Results for PR-GAMP obtained using the code available from https://sourceforge.net/projects/gampmatlab/
pr_gamp1 = c(0.95, 0.69, 0.34, 0.16, 0.04, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0) #, 0, 0, 0, 0)
pr_gamp2 = c(1, 1, 0.99, 1, 0.95, 0.91, 0.89, 0.65, 0.63, 0.56, 0.49, 0.37, 0.38, 0.22, 0.24) #, 0.2, 0.17, 0.08, 0.07)
pr_gamp3 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99)
pr_gamp4 = c(0, 0.01, 0.12, 0.33, 0.54, 0.78, 0.9, 0.98, 1, 1, 1, 1)

### Generate plot ###

pdf(file = "plot_xmax.pdf", width = 13, height = 3)

par(mfrow=c(1,3), mai = c(0.5, 0.7, 0.1, 0.1), bg = "transparent")

# Left plot: x_max = 1/sqrt(k)
plot(pr_gamp1, type = "o", lwd = 2, pch = 4, col = "brown", ylim = c(0,1), ylab = "", xlab = "", xaxt = "n", cex.axis = 1.4)
axis(1, at = (1 + 2*(0:11)), label = 10*(1:12), cex.axis = 1.4)
title(xlab = "Sparsity level k", ylab = "Success rate", line = 2.5, cex.lab = 1.7)
lines(apply(success_xmax, c(1,2,4), mean)[,1,2], type = "o", lwd = 2, pch = 18, col = "magenta")
lines(apply(success_xmax, c(1,2,4), mean)[,1,3], type = "o", lwd = 2, pch = 15, col = "blue")
lines(apply(success_xmax, c(1,2,4), mean)[,1,4], type = "o", lwd = 2, pch = 17, col = "black")
lines(apply(success_xmax, c(1,2,4), mean)[,1,1], type = "o", lwd = 2, pch = 16, col = "red")
legend('topright',legend=c("HWF", "SWF", "SPARTA", "PR-GAMP", "SPARTA-support"), col=c("red", "black", "blue", "brown", "magenta"), pch = c(16,17,15,4,18), lwd = 2, cex = 1.3, inset = c(0.01, -0.03), bty = "n")

# Middle plot: x_max = k^(-0.25)
plot(pr_gamp2, type = "o", lwd = 2, pch = 4, col = "brown", ylim = c(0,1), ylab = "", xlab = "", xaxt = "n", cex.axis = 1.4)
axis(1, at = (1 + 2*(0:11)), label = 10*(1:12), cex.axis = 1.4)
title(xlab = "Sparsity level k", ylab = "Success rate", line = 2.5, cex.lab = 1.7)
lines(apply(success_xmax, c(1,2,4), mean)[,2,2], type = "o", lwd = 2, pch = 18, col = "magenta")
lines(apply(success_xmax, c(1,2,4), mean)[,2,3], type = "o", lwd = 2, pch = 15, col = "blue")
lines(apply(success_xmax, c(1,2,4), mean)[,2,4], type = "o", lwd = 2, pch = 17, col = "black")
lines(apply(success_xmax, c(1,2,4), mean)[,2,1], type = "o", lwd = 2, pch = 16, col = "red")
legend('bottomleft',legend=c("HWF", "SWF", "SPARTA", "PR-GAMP", "SPARTA-support"), col=c("red", "black", "blue", "brown", "magenta"), pch = c(16,17,15,4,18), lwd = 2, cex = 1.3, inset = c(0.01, -0.03), bty = "n")

# Right plot: x_max = 0.7
plot(pr_gamp3, type = "o", lwd = 2, pch = 4, col = "brown", ylim = c(0,1), ylab = "", xlab = "", xaxt = "n", cex.axis = 1.4)
axis(1, at = (1 + 2*(0:11)), label = 10*(1:12), cex.axis = 1.4)
title(xlab = "Sparsity level k", ylab = "Success rate", line = 2.5, cex.lab = 1.7)
lines(apply(success_xmax, c(1,2,4), mean)[,3,2], type = "o", lwd = 2, pch = 18, col = "magenta")
lines(apply(success_xmax, c(1,2,4), mean)[,3,3], type = "o", lwd = 2, pch = 15, col = "blue")
lines(apply(success_xmax, c(1,2,4), mean)[,3,4], type = "o", lwd = 2, pch = 17, col = "black")
lines(apply(success_xmax, c(1,2,4), mean)[,3,1], type = "o", lwd = 2, pch = 16, col = "red")
legend('bottomleft',legend=c("HWF", "SWF", "SPARTA", "PR-GAMP", "SPARTA-support"), col=c("red", "black", "blue", "brown", "magenta"), pch = c(16,17,15,4,18), lwd = 2, cex = 1.3, inset = c(0.01, -0.03), bty = "n")

dev.off()