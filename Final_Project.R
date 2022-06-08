# Import package
# install.packages(c('usethis', 'grid', 'ggplot2', 'plot3D', 'boot', 'BB', 'devtools'))
library(grid)
library(usethis)
library(ggplot2)
library(plot3D)
library(boot)
library(BB)
library(devtools)

# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

# Plot the Real distribution
logistic <- function (beta0, beta1){
  x <- seq(-10, 10, 0.1)
  y <- 1/(1+exp(-beta0 - beta1 * x))
  
  return (cbind(x, y))
}

ggplot(data = NULL, aes(x = logistic(beta0 = 0, beta1 = 1)[ , 1], y = logistic(beta0 = 0, beta1 = 1)[ , 2])) + geom_smooth(color = 'red') + xlim(-3, 3) + ylim(0, 1)
ggplot(data = NULL, aes(x = logistic(beta0 = -0.5, beta1 = 0.5)[ , 1], y = logistic(beta0 = -0.5, beta1 = 0.5)[ , 2])) + geom_smooth(color = 'red') + xlim(-3, 3) + ylim(0, 1)




# 模拟出分布 beta0 = -0.5, beta_1 = 0.5
gdata <- function (n, beta0, beta1){
  x <- rnorm(n, mean = 0, sd = 1)
  pi <- 1 / (1 + exp(-beta0 - x * beta1))
  u <- runif(n)
  y <- rep(0, n)
  y[u < pi] = 1
  return (cbind(x, y))
}

ggplot(data = NULL, aes(x = gdata(1000, -0.5, 0.5)[ , 1], y = gdata(1000, -0.5, 0.5)[ , 2]))+ geom_point(col = "red")

gdata_mean <- function(n, m, beta0, beta1){
  x <- rnorm(n, mean = 0, sd = 1)
  pi <- 1 / (1 + exp(-beta0 - x * beta1))

  y <- matrix(data = rep(0, n * m), nrow = n, ncol = m)
  for ( i in 1:n){
    u <- runif(m)
    y[i, ][u < pi[i]] = 1
  }
  return (cbind(x, y))
}

get_F <- function(n, m, beta0, beta1){
  data <- gdata_mean(n, m, beta0, beta1)
  m <- dim(data)[2]
  data <- data[order(data[ , 1]), ]
  x <- data[ , 1]
  F <- as.matrix(apply(data[ , 1:m], 1, mean))
  return (cbind(x, F))
}

xF <- get_F(1000, 500, -0.5, 0.5)
ggplot(data = NULL) +
  geom_line(mapping = aes(x = xF[ , 1], y = xF[ , 2]), color = "blue") +
  geom_smooth(mapping = aes(x = logistic(-0.5, 0.5)[ , 1], y = logistic(-0.5, 0.5)[ , 2], color = 'red')) +
  xlim(-3, 3) + ylim(0, 1) + theme(legend.title=element_blank())





# Sensitivity analysis (via MSE)
n = 1000
beta0_array <- seq(-1.5,0.5,0.1)
beta1_array <- seq(-0.5,1.5,0.1)
list <- list()

Real_num <- function (x, beta0, beta1){
  return (1 / (1 + exp(-beta0 - beta1 * x)))
}

MSE <- function (n, x, F, beta0, beta1){
  return (sum((Real_num(x, beta0, beta1) - F)^2) / n)
}

for (beta0 in beta0_array){
  for (beta1 in beta1_array){
    xF <- get_F(1000, 500, beta0, beta1)
    x <- xF[ , 1]
    F <- xF[ , 2]
    mse <- c(MSE(n, x, F, beta0, beta1))
    list <- t(c(list, mse))
  }
}

mymatrix <- matrix(data = as.numeric(list),
                   nrow = length(beta0_array), ncol = length(beta1_array), byrow = TRUE)

rownames(mymatrix) <- round(seq(-1.5,0.5,0.1), 1)
colnames(mymatrix) <- round(seq(-0.5,1.5,0.1), 1)

write.csv(mymatrix, 'sensitivity_MSE.csv')

## plot the heatmap img
Heatmap(mymatrix, cluster_columns=FALSE, cluster_rows = FALSE, row_title = 'beta0', column_title = 'beta1')




# Nonpara - bootstrap
drive.bootstrap<-function(mdata,mlefunc, bootfunc, B){
      mle<-mlefunc(mdata)       # compute mle
      print(c("mle",mle))
      boots<-bootfunc(mdata,B)  # generate bootstrap samples
      print("quantiles")
      print(quantile(boots, c(0.005,0.025,0.5,0.975,0.995)))
      meanb<-mean(boots)
      print(c("mean",meanb))
      stderr<-sqrt(var(boots))
      print(c("stderr",stderr))
      biasc<-2*mle-meanb
      print(c("bias-corrected point estimate", biasc))


      cof <- length(boots[boots <= mle]) / B
      print(c("cof",cof))
      z0 <- qnorm(cof)      #corresponding normal quantile
      z.alpha <- qnorm(0.975)
      p.low.end <- pnorm(2*z0-z.alpha)
      p.high.end <- pnorm(2*z0+z.alpha)
      print(c(p.low.end,p.high.end))
      boots <- sort(boots)
      print("bias corrected C.I.")
      print(c(boots[B * p.low.end],boots[B * p.high.end]))
      
      
}

############## Function defining computation of \hat{\theta} ####

meanratio<-function(mdata){
   mean(mdata[ , 2])/var(mdata[, 1])
}

########### Function carrying out desired type of bootstrap ######

nonparboot.ratio<-function(mydat,B){
 # nonparametric bootstrap for ratio of means
 # data object must contain 2 columns of data
 # returns B bootstrap estimates of ratio of means

    if(ncol(mydat)!=2){
        print("input matrix must have 2 columns of numeric data")
    }
    else{
        bootratio<-numeric()
        n<-nrow(mydat)      # number of observation
        for (i in 1:B){
             index1<-sample(n,replace=T)
    # sample of size n from integers 1:n with replacement

             boot1<-mydat[index1,]
    # bootstrap sample; rows of data corresponding to index1

             bootratio<-c(bootratio, mean(boot1[ , 2])/var(boot1[ , 1]))
         }
         return(bootratio)
    }
}

drive.bootstrap(gdata(1000, -0.5, 0.5),meanratio,nonparboot.ratio,1000)




# Para - bootstrap
fun <- function(beta) {
  data <- gdata(1000, -0.5, 0.5)
  x <- data[ , 1]
  y <- data[ , 2]
  f <- numeric(length(beta))
  f[1] <- sum(y - exp(beta[1] + beta[2] * x) / (1 + exp(beta[1] + beta[2] * x)))
  f[2] <- sum(x * (y - exp(beta[1] + beta[2] * x) / (1 + exp(beta[1] + beta[2] * x))))
  f
}

beta <- dfsane(c(0, 0), fun, control=list(maxit=2000,trace = FALSE))$par

bootfunc <- function(B){
  beta0_mle <- c()
  beta1_mle <- c()
  for (i in 1:B){
    set.seed(i)
    beta_mle <- dfsane(c(0, 0),fun,control=list(maxit=2000,trace = FALSE))$par
    beta0_mle <- c(beta0_mle, beta_mle[1])
    beta1_mle <- c(beta1_mle, beta_mle[2])
  }

  return (list(beta0 = beta0_mle, beta1 = beta1_mle))
}

para_bootstrap <- function (data, mlefunc, bootfunc, B){
  mle <- mlefunc
  print(c("mle", mle))

  boots <- bootfunc(B)
  list0 <<- boots$beta0
  list1 <<- boots$beta1

  print("quantiles")
  print(quantile(boots$beta0, c(0.005,0.025,0.5,0.975,0.995)))
  print(quantile(boots$beta1, c(0.005,0.025,0.5,0.975,0.995)))


  print(c("mean0", mean(boots$beta0)))
  print(c("mean1", mean(boots$beta1)))
  mean0 <<- mean(boots$beta0)
  mean1 <<- mean(boots$beta1)


  print(c("stderr0", sqrt(var(boots$beta0))))
  print(c("stderr1", sqrt(var(boots$beta1))))

  biasc <- 2 * mle - c(mean(boots$beta0), mean(boots$beta1))
  print("bias-corrected pointed estimate")
  print(biasc)

  boots0 <- boots$beta0
  boots1 <- boots$beta1
  mle0 <- mle[1]
  mle1 <- mle[2]

  ### beta0
  cof <- length(boots0[boots0 <= mle0]) / B
  print(c("cof",cof))
  z0 <- qnorm(cof)      #corresponding normal quantile
  z.alpha <- qnorm(0.975)
  p.low.end <- pnorm(2*z0-z.alpha)
  p.high.end <- pnorm(2*z0+z.alpha)
  print(c(p.low.end,p.high.end))
  boots0 <- sort(boots0)
  print("bias corrected C.I.")
  print(c(boots0[B * p.low.end],boots0[B * p.high.end]))

  ### beta1
  cof <- length(boots1[boots1 <= mle1]) / B
  print(c("cof",cof))
  z0 <- qnorm(cof)      #corresponding normal quantile
  z.alpha <- qnorm(0.975)
  p.low.end <- pnorm(2*z0-z.alpha)
  p.high.end <- pnorm(2*z0+z.alpha)
  print(c(p.low.end,p.high.end))
  boots1 <- sort(boots1)
  print("bias corrected C.I.")
  print(c(boots1[B * p.low.end],boots1[B * p.high.end]))
  
}

para_bootstrap(data, beta, bootfunc, 100)

ggplot(data = NULL, aes(x = seq(1, 100), y = list0)) + geom_point(color = 'red') +
      geom_smooth(mapping = aes(x = seq(1, 100), y = mean0))
ggplot(data = NULL, aes(x = seq(1, 100), y = list1)) + geom_point(color = 'red') +
      geom_smooth(mapping = aes(x = seq(1, 100), y = mean1))

mgdata <- function (n, beta){
  m <- length(beta) - 1
  x <- matrix(rnorm(m * n, 0, 1), nrow = n)
  X <- cbind(matrix(rep(1, n), nrow = n, ncol = 1), x)

  pi <- 1 / (1 + exp(-X %*% beta))
  u <- runif(n)
  y <- rep(0, n)
  y[u < pi] = 1

  return (cbind(X, y))
}
mdata <- mgdata(100, c(1, 2, 3))
x <- mdata[ , 2]
y <- mdata[ , 3]
z <- mdata[ , 4]
df <- data.frame(x, y, z)
scatter3D(x, y, z, ticktype = "detailed",bty = "b2")
ggplot(data = df, aes(x = x, y = y)) + geom_point(aes(colour = z)) + scale_x_discrete("x1") + scale_y_discrete("x2")

fun <- function(beta) {
  mdata <- mgdata(1000, c(-0.5, 0.5, 1))
  m <- dim(mdata)[2] - 2
  X <- mdata[ , c(1:(m + 1))]
  y <- mdata[ , dim(mdata)[2]]
  f <- numeric(length(beta))
  for (i in 1:(m + 1))
    f[i] <- (t(X) %*% (y - exp(X %*% beta) / (1 + exp(X %*% beta))))[i, ]
  f
}
beta <- dfsane(c(0, 0, 0), fun, control=list(maxit=2500,trace = FALSE))$par
beta