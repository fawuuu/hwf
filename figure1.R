##### Code for Figure 1 ##### 

library(foreach)
library(doParallel)
library(abind)

source("algorithms.R")

acomb <- function(...) abind(..., along=3)

### Experiments for Figure 1 ###
n = 10000         # Dimension of signal vector x
k = 5 + 5*(1:23)  # Sparsity level of x
m = 5000          # Number of observations
num_rep = 100     # Number of Monte Carlo trials

# Parallelize computation
numcl = 10
cl <- makeCluster(numcl)
registerDoParallel(cl)
clusterExport(cl,list("gen_data_sparse", "hwf", "wf_loss", "wf_grad"))
clusterExport(cl,list("n","k", "m", "num_rep"),envir=environment())

# Run different support recovery algorithms num_rep times with the specified parameters
success <- foreach(rep = 1:num_rep, .combine='acomb', .multicombine=TRUE) %dopar%{
  B = array(0, dim = c(length(k), 4, 1, 3))
  for(j in 1:length(k)){
    for(i in 1:4){
      # Generate data
      data = gen_data_sparse(m, n, k[j], xmax = i)
      
      # Support recovery using HWF
      ind_sup = sort(t(data$A^2) %*% data$y, index.return = TRUE)$ix[n]
      res = hwf(data$A, data$y, ini = ind_sup, iteration = 2)[2,]
      # Proportion of correctly identified support variables using HWF
      B[j,i,1,1] = sum(sort(abs(res), index.return = TRUE)$ix[(n-k[j]+1):n] %in% which(data$x!=0)) / k[j]
     
      # Support recovery of SPARTA
      B[j,i,1,2] = sum(sort(t(data$A^2) %*% data$y, index.return = TRUE)$ix[(n-k[j]+1):n] %in% which(data$x!=0)) / k[j]
      
      # Support recovery of SparseAltMinPhase
      B[j,i,1,3] = sum(sort(t(abs(data$A)) %*% sqrt(data$y), index.return = TRUE)$ix[(n-k[j]+1):n] %in% which(data$x!=0)) / k[j]
    }
  }
  B
}
stopCluster(cl)

# Save results
saveRDS(success, file = "output_support.rds")


### Generate plot ###

# Compute mean and standard deviation of proportion of support variables recovered
support_mean = apply(success, c(1,2,4), mean)
support_sd = apply(success, c(1,2,4), sd)

pdf(file = "plot_support.pdf", width = 14, height = 2.5)

par(mfrow=c(1,4), mai = c(0.5, 0.7, 0.1, 0.1), bg = "transparent")

# Left plot: x_max = 1/sqrt(k)
plot(support_mean[,2,1], type = "o", lwd = 2, pch = 16, col = "red", ylim = c(0,1), ylab = "", xlab = "", xaxt = "n", cex.axis = 1.45)
axis(1, at = (1 + 2*(0:11)), label = 10*(1:12), cex.axis = 1.45)
title(xlab = "Sparsity level k", ylab = "Prop. of support recovered", line = 2.5, cex.lab = 1.7)
lines(support_mean[,2,2], type = "o", lwd = 2, pch = 15, col = "blue")
lines(support_mean[,2,3], type = "o", lwd = 2, pch = 17, col = "black")
legend('topright',legend=c("HWF", "SPARTA", "SAMP"), col=c("red", "blue", "black"), pch = c(16,15,17), lwd = 2, cex = 1.5, inset = c(0, -0.025), bty = "n")
# Add error bars
conf_low = support_mean[,2,1] - support_sd[,2,1]
conf_low[conf_low<0] = 0
conf_up = support_mean[,2,1] + support_sd[,2,1]
conf_up[conf_up>1] = 1
arrows(x0=1:23, y0=conf_low, x1=1:23, y1=conf_up, length=0, col = "red", lwd = 2)
arrows(x0=1:23, y0=support_mean[,2,2] - support_sd[,2,2],
       x1=1:23, y1=support_mean[,2,2] + support_sd[,2,2], length=0, col = "blue", lwd = 2)
arrows(x0=1:23, y0=support_mean[,2,3] - support_sd[,2,3],
       x1=1:23, y1=support_mean[,2,3] + support_sd[,2,3], length=0, col = "black", lwd = 2)

# Second from left plot: x_max = k^(-0.25)
plot(support_mean[,3,1], type = "o", lwd = 2, pch = 16, col = "red", ylim = c(0,1), ylab = "", xlab = "", xaxt = "n", cex.axis = 1.45)
axis(1, at = (1 + 2*(0:11)), label = 10*(1:12), cex.axis = 1.45)
title(xlab = "Sparsity level k", ylab = "Prop. of support recovered", line = 2.5, cex.lab = 1.7)
lines(support_mean[,3,2], type = "o", lwd = 2, pch = 15, col = "blue")
lines(support_mean[,3,3], type = "o", lwd = 2, pch = 17, col = "black")
legend('topright',legend=c("HWF", "SPARTA", "SAMP"), col=c("red", "blue", "black"), pch = c(16,15,17), lwd = 2, cex = 1.5, inset = c(0, -0.025), bty = "n")
# Add error bars
arrows(x0=1:23, y0=support_mean[,3,1] - support_sd[,3,1],
       x1=1:23, y1=support_mean[,3,1] + support_sd[,3,1], length=0, col = "red", lwd = 2)
arrows(x0=1:23, y0=support_mean[,3,2] - support_sd[,3,2],
       x1=1:23, y1=support_mean[,3,2] + support_sd[,3,2], length=0, col = "blue", lwd = 2)
arrows(x0=1:23, y0=support_mean[,3,3] - support_sd[,3,3],
       x1=1:23, y1=support_mean[,3,3] + support_sd[,3,3], length=0, col = "black", lwd = 2)

# Second from right plot: x_max = 0.7
plot(support_mean[,4,1], type = "o", lwd = 2, pch = 16, col = "red", ylim = c(0,1), ylab = "", xlab = "",xaxt = "n", cex.axis = 1.45)
axis(1, at = (1 + 2*(0:11)), label = 10*(1:12), cex.axis = 1.45)
title(xlab = "Sparsity level k", ylab = "Prop. of support recovered", line = 2.5, cex.lab = 1.7)
lines(support_mean[,4,2], type = "o", lwd = 2, pch = 15, col = "blue")
lines(support_mean[,4,3], type = "o", lwd = 2, pch = 17, col = "black")
legend('topright',legend=c("HWF", "SPARTA", "SAMP"), col=c("red", "blue", "black"), pch = c(16,15,17), lwd = 2, cex = 1.5, inset = c(0, -0.025), bty = "n")
# Add error bars
arrows(x0=1:23, y0=support_mean[,4,1] - support_sd[,4,1],
       x1=1:23, y1=support_mean[,4,1] + support_sd[,4,1], length=0, col = "red", lwd = 2)
arrows(x0=1:23, y0=support_mean[,4,2] - support_sd[,4,2],
       x1=1:23, y1=support_mean[,4,2] + support_sd[,4,2], length=0, col = "blue", lwd = 2)
arrows(x0=1:23, y0=support_mean[,4,3] - support_sd[,4,3],
       x1=1:23, y1=support_mean[,4,3] + support_sd[,4,3], length=0, col = "black", lwd = 2)

# Right plot: no restrictions on x_max
plot(support_mean[,1,1], type = "o", lwd = 2, pch = 16, col = "red", ylim = c(0,1), ylab = "", xlab = "",xaxt = "n", cex.axis = 1.45)
axis(1, at = (1 + 2*(0:11)), label = 10*(1:12), cex.axis = 1.45)
title(xlab = "Sparsity level k", ylab = "Prop. of support recovered", line = 2.5, cex.lab = 1.7)
lines(support_mean[,1,2], type = "o", lwd = 2, pch = 15, col = "blue")
lines(support_mean[,1,3], type = "o", lwd = 2, pch = 17, col = "black")
legend('topright',legend=c("HWF", "SPARTA", "SAMP"), col=c("red", "blue", "black"), pch = c(16,15,17), lwd = 2, cex = 1.5, inset = c(0, -0.025), bty = "n")
# Add error bars
arrows(x0=1:23, y0=support_mean[,1,1] - support_sd[,1,1],
       x1=1:23, y1=support_mean[,1,1] + support_sd[,1,1], length=0, col = "red", lwd = 2)
arrows(x0=1:23, y0=support_mean[,1,2] - support_sd[,1,2],
       x1=1:23, y1=support_mean[,1,2] + support_sd[,1,2], length=0, col = "blue", lwd = 2)
arrows(x0=1:23, y0=support_mean[,1,3] - support_sd[,1,3],
       x1=1:23, y1=support_mean[,1,3] + support_sd[,1,3], length=0, col = "black", lwd = 2)

dev.off()