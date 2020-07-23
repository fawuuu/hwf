##### Code for Figure 4 #####

library(foreach)
library(doParallel)
library(abind)
library(gplots)

source("algorithms.R")

acomb <- function(...) abind(..., along=3)

### Experiments for Figure 4 ###
n = 500           # Dimension of signal vector x
k = 5 + 5*(1:23)  # Sparsity level of x
m = 100 * (1:10)  # Number of observations 
num_rep = 100     # Number of Monte Carlo trials
num_try = 50      # Maximum number of allowed restarts
prec = 0.01       # Precision required to declare a run successful

# Parallelize computation
numcl = 10
cl <- makeCluster(numcl)
registerDoParallel(cl)
clusterExport(cl,list("gen_data_sparse","hwf", "wf_loss", "wf_grad"))
clusterExport(cl,list("n","k", "m", "num_rep", "num_try", "prec"),envir=environment())

# Run HWF num_rep times with the specified parameters
success <- foreach(rep = 1:num_rep, .combine='acomb', .multicombine=TRUE) %dopar%{
  B = array(0, dim = c(length(m), length(k), 1))
  for(i in 1:length(m)){
    for(j in 1:length(k)){
      # Generate data
      data = gen_data_sparse(m[i], n, k[j])
      
      # Signal reconstruction using HWF
      ind_sup = sort(t(data$A^2) %*% data$y, index.return = TRUE)$ix[(n-num_try+1):n]
      for(trial in 1:num_try){
        # To save computation, only run HWF when we pick a coordinate on the support
        if(data$x[ind_sup[num_try+1-trial]]!=0){
          res = hwf(data$A, data$y, ini = ind_sup[num_try+1-trial], iteration = 100000)
          if(min(sqrt(sum((res[nrow(res),]-data$x)^2)), sqrt(sum((res[nrow(res),]+data$x)^2))) < prec){
            B[i, j, 1] = 1
            break
          }
        }
      }
    }
  }
  B
}
stopCluster(cl)

# Save result
saveRDS(success, file = "output_hwf_heatmap.rds")

### Generate plot ###

pdf(file = "plot_heatmap.pdf", width = 7, height = 3)

result = data.frame(apply(success, c(1,2), mean))

row.names(result) <- 100 * (1:10)
colnames(result) <- 5 + 5 * (1:23)

# Estimate the mean of x_max of our signal
kxmax = rep(0, 23)
for(i in 1:23){
  for(j in 1:100000){
    temp = rnorm(5 + 5*i)
    temp = temp / sqrt(sum(temp^2))
    kxmax[i] = kxmax[i] + max(abs(temp)) / 100000 
  }
}

# Define the curve m = 1/3 * k(x^*_max)^(-2)log(n/k)
pointdata <- data.frame(
  xname = 1:23, 
  ypos = (5 + 5*(1:23)) * log(500/(5+5*(1:23))) / (3 * kxmax^2)  + 100
) 

result2 <- result %>% rownames_to_column() %>% gather(sparsity,  p_success, -rowname)

level_sparsity <- c('10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70',
                    '75', '80', '85', '90', '95', '100', '105', '110', '115', '120')
level_obs <- c('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000')

ggplot() + geom_tile(result2, mapping = aes(x = factor(sparsity, levels = level_sparsity), y = as.numeric(rowname, levels = level_obs),
                                            fill = p_success)) + geom_tile() + scale_fill_continuous(low="blue",high="red") + 
  xlab("Sparsity level k") + ylab("Number of observations m") +
  scale_x_discrete(limits=level_sparsity, breaks=level_sparsity[seq(1,length(level_sparsity),by=2)]) +
  geom_line(data = pointdata, mapping = aes(xname, ypos), size = 0.5, col = "black") +
  scale_y_continuous(breaks=c(100,200,300,400,500,600,700,800,900,1000),expand=c(0,0)) + 
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA), 
        axis.text=element_text(size=11), axis.title=element_text(size=13))

dev.off()