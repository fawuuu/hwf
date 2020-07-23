##### Code for Figure 6 #####

source("algorithms.R")

# Define the function to generate complex data
# m: number of observations
# n: dimension of the signal x^*
# k: sparsity level
gen_data_sparse_complex = function(m, n, k){
  # Generate complex k-sparse Gaussian vector x
  x = rnorm(n) + rnorm(n) * sqrt(rep(as.complex(-1), n))
  ind = sample(1:n, k)
  x[-ind] = 0
  # Normalize to ||x||_2=1
  norm = sqrt(as.numeric(x %*% Conj(x)))
  x = x / norm
  
  # Generate complex-Gaussian measurement matríx A
  A = matrix(rnorm(m*n, sd = sqrt(0.5)) + 
               rnorm(m*n, sd = sqrt(0.5)) * sqrt(rep(as.complex(-1),m*n)),
             nrow = m, ncol = n)
  
  # Generate (real) vector of observations y
  y = Re((A %*% Conj(x)) * Conj(A %*% Conj(x)))
  
  return(list(A = A, x = x, y = y))
}

# Define complex HWF
# A: measurement matrix
# y: vector of observations
# ini: index for initialization (i_max)
# step: step size
# alpha: initialization size
# iteration: maximum number of iterations
hwf_complex = function(A, y, ini = 1, step = 0.1, alpha = 1e-3, iteration = 100){
  m = nrow(A)
  n = ncol(A)
  
  X = matrix(as.complex(0), nrow = iteration, ncol = n)
  U = matrix(as.complex(0), nrow = iteration, ncol = n)
  V = matrix(as.complex(0), nrow = iteration, ncol = n)
  
  # Initialize u_j, v_j = alpha + alpha*i
  u_cur = rep(alpha + alpha * sqrt(as.complex(-1)), n)
  v_cur = rep(alpha + alpha * sqrt(as.complex(-1)), n)
  
  # For i_max, set u_{i_max} = 3^(-0.25) + alpha * i
  u_cur[ini] = 3^(-0.25) + alpha * sqrt(as.complex(-1))
  
  # Parametrization: x_j = a_j + i*b_j, where 
  # a_j = Re(u_j)^2 - Re(v_j)^2 and b_j = Im(u_j)^2 - Im(v_j)^2
  x_cur = Re(u_cur)^2 - Re(v_cur)^2 + (Im(u_cur)^2 - Im(v_cur)^2) * sqrt(as.complex(-1))
  
  U[1,] = u_cur
  V[1,] = v_cur
  X[1,] = x_cur
  
  for(t in 2:iteration){
    # Current estimate of y using x_cur
    est = Re((A %*% Conj(x_cur)) * Conj(A %*% Conj(x_cur)))
    # Compute the gradient
    r = 2 * step * t(A) %*% ((est - y) * (Conj(A) %*% x_cur)) / m
    
    # Perform the gradient updates on u and v
    u_cur = Re(u_cur) * (1 - Re(r)) + (Im(u_cur) * (1 - Im(r))) * sqrt(as.complex(-1))
    v_cur = Re(v_cur) * (1 + Re(r)) + (Im(v_cur) * (1 + Im(r))) * sqrt(as.complex(-1))
    x_cur = Re(u_cur)^2 - Re(v_cur)^2 + (Im(u_cur)^2 - Im(v_cur)^2)* sqrt(as.complex(-1))
    
    # Save new estimate
    U[t,] = u_cur
    V[t,] = v_cur
    X[t,] = x_cur
  }
  
  return(X)
}

# Define distance between two complex numbers for a given phase phi
dist = function(phi, x, z){
  dif = exp(sqrt(as.complex(-1))*phi) * z - x
  return(sqrt(Re(dif %*% Conj(dif))))
}

### Experiment for Figure 6 ###
set.seed(1)
n = 500       # Dimension of signal vector x
k = 10        # Sparsity level of x
m = 500       # Number of observations
iter = 10000  # Number of iterations

# HWF in the real model
data = gen_data_sparse(n, m, k)
ini1 = sort(t(data$A^2) %*% data$y, index.return = TRUE)$ix[n]
res = hwf(data$A, data$y, ini = ini1, alpha = 1e-4, eps = 0, iteration = iter)

# HWF in the complex model
datac = gen_data_sparse_complex(n, m, k)
ini2 = sort(t(Re(datac$A * Conj(datac$A))) %*% datac$y, index.return = TRUE)$ix[n]
resc = hwf_complex(datac$A, datac$y, ini = ini2, alpha = 1e-4, iteration = iter)

# Compute the distance from the true signal
error = rep(0, iter)
errorc = rep(0, iter)
for(i in 1:iter){
  error[i] = min(sqrt(sum((res[i,]-data$x)^2)), sqrt(sum((res[i,]+data$x)^2))) 
  errorc[i] = optimize(dist, x = datac$x, z = resc[i,], interval = c(0, 2*pi))$objective
}


### Generate plot ###

pdf(file = "plot_conv_complex.pdf", width = 8, height = 3)

# Left plot: relative error for 10000 iterations
par(mfrow=c(1,2), mai = c(0.7, 0.7, 0.1, 0.1), bg = "transparent")
plot(log(error), type = "l", lwd = 2, col = "red", ylab = "", xlab = "", xaxt = "n")
axis(1, at = 1000*(0:10), label = 1000*(0:10))
title(xlab = "Iteration t", ylab = "Relative error (log)", line = 2.2, cex.lab = 1.2)
lines(log(errorc), col = "blue", lwd = 2)
legend('topright',legend=c("HWF-real", "HWF-complex"), col=c("red", "blue"), lwd = 2, cex = 0.9, inset = c(0.025, 0.025), bty = "n")

# Right plot: relative error for 500 iterations
plot(log(error)[1:500], type = "l", lwd = 2, col = "red", ylab = "", xlab = "", xaxt = "n", ylim = c(-3.8,0))
axis(1, at = 50*(0:10), label = c(0,"",100,"",200,"",300,"",400,"",500))
title(xlab = "Iteration t", ylab = "Relative error (log)", line = 2.2, cex.lab = 1.2)
lines(log(errorc)[1:500], col = "blue", lwd = 2)
legend('bottomleft',legend=c("HWF-real", "HWF-complex"), col=c("red", "blue"), lwd = 2, cex = 0.9, inset = c(0.025, 0.025), bty = "n")

dev.off()