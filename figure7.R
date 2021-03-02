##### Code for Figure 7 #####

source("algorithms.R")

# Define HWF with random initialization
# A: measurement matrix
# y: vector of observations
# step: step size
# alpha: initialization size
# iteration: maximum number of iterations
hwf_rand = function(A, y, step = 0.1, alpha = 1e-3, iteration = 100){
  m = nrow(A)
  n = ncol(A)
  
  X = matrix(0, nrow = iteration, ncol = n)
  
  # Initialization 
  # Estimate of the signal size
  theta = sqrt(mean(y))
  
  u_cur = rnorm(n, sd = alpha)
  v_cur = rnorm(n, sd = alpha)
  
  x_cur = u_cur^2 - v_cur^2
 
  X[1,] = x_cur
  
  for(t in 2:iteration){
    #print(t)
    # Gradient updates
    r = 2 * step * wf_grad(x_cur, A, y)
    
    u_cur = u_cur * (1 - r)
    v_cur = v_cur * (1 + r)
    x_cur = u_cur^2 - v_cur^2
    
    X[t,] = x_cur
  }
  
  return(X)
}

### Experiment for Figure 6 ###
set.seed(1)

n = 1000        # Dimension of signal vector x
k = 10          # Sparsity level of x
m = 700         # Number of observations

# Run HWF for 3000 iterations
iter = 3000

# Generate data
data = gen_data_sparse(m, n, k)

# HWF with our proposed initialization
ini1 = sort(t(data$A^2) %*% data$y, index.return = TRUE)$ix[n]
res = hwf(data$A, data$y, ini = ini1, alpha = 1e-3, eps = 0, iteration = iter)

# HWF with random initialization
res_rand = hwf_rand(data$A, data$y, alpha = 0.01, iteration = iter)

# Compute the distance from the true signal
error = rep(0, iter)
error_rand = rep(0, iter)
for(i in 1:iter){
  error[i] = min(sqrt(sum((res[i,]-data$x)^2)), sqrt(sum((res[i,]+data$x)^2))) 
  error_rand[i] = min(sqrt(sum((res_rand[i,]-data$x)^2)), sqrt(sum((res_rand[i,]+data$x)^2)))
}

### Generate plot ###

pdf(file = "plot_randini.pdf", width = 6, height = 4)

par(mfrow=c(1,1), mai = c(0.7, 0.7, 0.1, 0.1), bg = "transparent")
plot(log(error), type = "l", lwd = 2, col = "red", ylim = c(-6.5, 0), ylab = "", xlab = "", xaxt = "n")
axis(1, at = 500*(0:6), label = 500*(0:6), cex.axis = 1.3)
title(xlab = "Iteration t", ylab = "Relative error (log)", line = 2.2, cex.lab = 1.4)
lines(log(error_rand), col = "blue", lwd = 2)
legend('topright',legend=c("Our initialization", "Random initialization"), col=c("red", "blue"), lwd = 2, cex = 1.2, inset = c(0.025, 0.025), bty = "n")

dev.off()