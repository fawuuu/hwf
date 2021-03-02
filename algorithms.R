##### Code for "Hadamard Wirtinger Flow for Sparse Phase Retrieval" #####

# Define the function to generate data
# m: number of observations
# n: dimension of the signal x^*
# k: sparsity level
# xmax: Which restriction on x_max to enforce
# maxval: If xmax = 5, set value for x_max directly
gen_data_sparse = function(m, n, k, xmax = 1, maxval = 1){
  # Generate a k-sparse Gaussian vector x
  x = rnorm(n)
  ind = sample(1:n, k)
  x[-ind] = 0
  
  if(xmax == 1){  # No restrictions on x_max
    x = x / sqrt(sum(x^2))
  }
  
  if(xmax == 2){  # x_max = 1/sqrt(k)
    x[ind] = rep(1/sqrt(k), k) * sample(c(-1,1), k, replace = T)
  }
  
  if(xmax == 3){  # x_max = k^(-0.25)
    ind_max = sample(ind, 1)
    x[ind_max] = k^(-0.25)
    # Normalize to ||x||_2=1 
    x[-ind_max] = x[-ind_max] * sqrt((1 - x[ind_max]^2) / sum(x[-ind_max]^2))
    # Ensure that x_max is indeed the maximum component
    while(any(abs(x[-ind_max])>x[ind_max])){
      ind_decr = which(abs(x)>x[ind_max])
      x[ind_decr] = runif(length(ind_decr), 1/sqrt(k), x[ind_max]) * sign(x[ind_decr])
      x[-ind_max] = x[-ind_max] * sqrt((1 - x[ind_max]^2) / sum(x[-ind_max]^2))
    }
  }
  
  if(xmax == 4){  # x_max = 0.7
    ind_max = sample(ind, 1)
    x[ind_max] = 0.7
    # Normalize to ||x||_2=1 
    x[-ind_max] = x[-ind_max] * sqrt((1 - x[ind_max]^2) / sum(x[-ind_max]^2))
    # Ensure that x_max is indeed the maximum component
    while(any(abs(x[-ind_max])>x[ind_max])){
      ind_decr = which(abs(x)>x[ind_max])
      x[ind_decr] = runif(length(ind_decr), 1/sqrt(k), x[ind_max]) * sign(x[ind_decr])
      x[-ind_max] = x[-ind_max] * sqrt((1 - x[ind_max]^2) / sum(x[-ind_max]^2))
    }
  }
  
  if(xmax == 5){  # x_max = maxval (to specify)
    ind_max = sample(ind, 1)
    x[ind_max] = maxval
    # Normalize to ||x||_2=1 
    x[-ind_max] = x[-ind_max] * sqrt((1 - x[ind_max]^2) / sum(x[-ind_max]^2))
    # Ensure that x_max is indeed the maximum component
    if(maxval > 1/sqrt(k)){   # x_max = maxval only possible if maxval > 1/sqrt(k)
      while(any(abs(x[-ind_max])>x[ind_max])){
        ind_decr = which(abs(x)>x[ind_max])
        x[ind_decr] = runif(length(ind_decr), 1/sqrt(k), x[ind_max]) * sign(x[ind_decr])
        x[-ind_max] = x[-ind_max] * sqrt((1 - x[ind_max]^2) / sum(x[-ind_max]^2))
      }
    }
  }
  
  A = matrix(rnorm(m*n), nrow = m, ncol = n)
  
  y = as.vector((A%*%x)^2)
  
  return(list(A = A, x = x, y = y))
}

# Define the empirical risk, squared-magnitude-based loss
wf_loss = function(x, A, y){
  m = length(y)
  return(sum((y - (A %*% x)^2)^2) / (4*m))
}

# Define the gradient, squared-magnitude-based loss 
wf_grad = function(x, A, y){
  m = nrow(A)
  return(t(A) %*% (((A%*%x)^2 - y) * (A%*%x)) / m)
}

# Define the gradient, magnitude-based loss
rwf_grad = function(x, A, y){
  m = nrow(A)
  est = A %*% x
  return(t(A) %*% (est - y * sign(est)) / m)
}

# Define HWF
# A: measurement matrix
# y: vector of observations
# ini: index for initialization (i_max)
# step: step size
# alpha: initialization size
# eps: stopping criterion 
# iteration: maximum number of iterations
hwf = function(A, y, ini, step = 0.1, alpha = 1e-3, eps = 1e-7, iteration = 100){
  m = nrow(A)
  n = ncol(A)
  X = matrix(0, nrow = iteration, ncol = n)
  
  # Initialization 
  u_cur = rep(alpha, n)
  v_cur = rep(alpha, n)
  
  u_cur[ini] = sqrt(sqrt(1/3) - alpha^2) 
  
  x_cur = u_cur^2 - v_cur^2
  X[1,] = x_cur
  
  for(t in 1:iteration){
    # Gradient updates
    r = 2 * step * wf_grad(x_cur, A, y)
    
    u_cur = u_cur * (1 - r)
    v_cur = v_cur * (1 + r)
    x_cur = u_cur^2 - v_cur^2
    
    X[t,] = x_cur
    
    # Stop algorithm if estimates diverge
    if(sum(x_cur^2)>1e10){
      X[t,] = rep(0, n)
      X = X[1:t,]
      break
    } 
    
    # Stop algorithm if objective function is smaller than eps 
    # Check every 1000 iterations to save computation
    if(t/1000 == round(t/1000)){
      if(wf_loss(x_cur, A, y) < eps){
        X = X[1:t,]
        break
      } 
    }
  }

  return(X)
}

# Define SPARTA-support
# trial: number of allowed restarts
# trunc_thresh, gamma: parameters of SPARTA
sparta_sup = function(A, y, k, step = 0.1, trial = 1, trunc_thresh = 1/6, gamma = 1, eps = 1e-7, iteration = 100){
  m = nrow(A)
  n = ncol(A)
  X = matrix(0, nrow = iteration, ncol = n)
  
  # Initialization - support recovery using HWF
  ind_ini = sort(t(A^2) %*% y, index.return = TRUE)$ix[n-trial+1]
  res = hwf(A, y, ini = ind_ini, iteration = 2)[2,]
  ind_sup = sort(abs(res), index.return = TRUE)$ix[(n-k+1):n] 
  
  # Orthogonality promoting initialization from SPARTA
  trunc = ceiling(trunc_thresh * m)
  ind_ini = sort(sqrt(y) / sqrt(apply(A[,ind_sup]^2, 1, sum)), index.return = TRUE)$ix[(m-trunc+1):m]
  
  Y = t(A[ind_ini,ind_sup]) %*% (A[ind_ini,ind_sup] / apply(A[ind_ini,ind_sup]^2, 1, sum)) / trunc
  
  v <- rep(1/sqrt(k),k)
  for (s in 1:100) {
    v <- Y %*% v
    v <- v / sqrt(sum(v^2))
  }
  
  x_cur = rep(0, n)
  x_cur[ind_sup] = sqrt(mean(y)) * v
  X[1,] = x_cur
  
  for(t in 2:iteration){
    # Thresholded gradient updates from SPARTA
    ind_update = which(abs(A %*% x_cur) > sqrt(y) / (1 + gamma))
    if(any(abs(A %*% x_cur) > sqrt(y) / (1 + gamma))){
      x_cur = x_cur - step * rwf_grad(x_cur, A[ind_update,], sqrt(y[ind_update])) * length(ind_update) / m
      x_cur[sort(abs(x_cur), index.return = TRUE)$ix[1:(n-k)]] = 0 
    }
    X[t,] = x_cur
    
    # Stop algorithm if estimates diverge
    if(sum(x_cur^2)>1e10){
      X[t,] = rep(0, n)
      X = X[1:t,]
      break
    } 
    
    # Stop algorithm if objective function is smaller than eps 
    # Check every 1000 iterations to save computation
    if(t/1000 == round(t/1000)){
      if(wf_loss(x_cur, A, y) < eps){
        X = X[1:t,]
        break
      } 
    }
  }
  
  return(X)
}

# Define SPARTA
sparta = function(A, y, k, step = 1, trunc_thresh = 1/6, gamma = 1, eps = 1e-7, iteration = 100){
  m = nrow(A)
  n = ncol(A)
  X = matrix(0, nrow = iteration, ncol = n)
    
  #Initialization - support recovery
  ind_sup = sort(t(A^2) %*% y, index.return = TRUE)$ix[(n-k+1):n]
  
  # Orthogonality promoting initialization
  trunc = ceiling(trunc_thresh * m)
  ind_ini = sort(sqrt(y) / sqrt(apply(A[,ind_sup]^2, 1, sum)), index.return = TRUE)$ix[(m-trunc+1):m]
  
  Y = t(A[ind_ini,ind_sup]) %*% (A[ind_ini,ind_sup] / apply(A[ind_ini,ind_sup]^2, 1, sum)) / trunc
  
  v <- rep(1/sqrt(k),k)
  for (s in 1:100) {
    v <- Y %*% v
    v <- v / sqrt(sum(v^2))
  }
  
  x_cur = rep(0, n)
  x_cur[ind_sup] = sqrt(mean(y)) * v
  X[1,] = x_cur
  
  for(t in 2:iteration){
    # Thresholded gradient updates
    ind_update = which(abs(A %*% x_cur) > sqrt(y) / (1 + gamma))
    if(any(abs(A %*% x_cur) > sqrt(y) / (1 + gamma))){
      x_cur = x_cur - step * rwf_grad(x_cur, A[ind_update,], sqrt(y[ind_update])) * length(ind_update) / m
      x_cur[sort(abs(x_cur), index.return = TRUE)$ix[1:(n-k)]] = 0 
    }
    X[t,] = x_cur
    
    # Stop algorithm if estimates diverge
    if(sum(x_cur^2)>1e10){
      X[t,] = rep(0, n)
      X = X[1:t,]
      break
    } 
    
    # Stop algorithm if objective function is smaller than eps 
    # Check every 1000 iterations to save computation
    if(t/1000 == round(t/1000)){
      if(wf_loss(x_cur, A, y) < eps){
        X = X[1:t,]
        break
      } 
    }
  }
  
  return(X)
}

# Define SWF
swf = function(A, y, k, alpha_y = 3, eps = 1e-7, iteration = 100){
  m = nrow(A)
  n = ncol(A)
  X = matrix(0, nrow = iteration, ncol = n)
  
  #Initialization - support recovery
  ind_sup = sort(t(A^2) %*% y, index.return = TRUE)$ix[(n-k+1):n]
  
  # Spectral initialization
  phi = mean(y)
  ind_ini = which(abs(y) <= alpha_y^2 * phi)
  
  Y = t(A[ind_ini,ind_sup]) %*% (A[ind_ini,ind_sup] * y[ind_ini]) / m
  
  v <- rep(1/sqrt(k),k)
  for (s in 1:100) {
    v <- Y %*% v
    v <- v / sqrt(sum(v^2))
  }
  
  x_cur = rep(0, n)
  x_cur[ind_sup] = sqrt(phi) * v
  X[1,] = x_cur
  
  for(t in 2:iteration){
    # Computation of the (increasing) sequence of step sizes
    step = min(0.5 * (1-exp(-t/330)), 0.1)
    # Thresholded gradient updates
    x_cur = x_cur - step/phi * wf_grad(x_cur, A, y)
    x_cur[sort(abs(x_cur), index.return = TRUE)$ix[1:(n-k)]] = 0 
    X[t,] = x_cur
    
    # Stop algorithm if estimates diverge
    if(sum(x_cur^2)>1e10){
      X[t,] = rep(0, n)
      X = X[1:t,]
      break
    } 
    
    # Stop algorithm if objective function is smaller than eps 
    # Check every 1000 iterations to save computation
    if(t/1000 == round(t/1000)){
      if(wf_loss(x_cur, A, y) < eps){
        X = X[1:t,]
        break
      } 
    }
  }
  
  return(X)
}