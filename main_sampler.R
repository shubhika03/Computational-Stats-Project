
# for different use case different implementation of this function
# should be done, but the function definition i.e parameters and return
# value shouldn't change

pdf <- function(x){
  return(0.5*dnorm(x, mean=-1, sd=1) + 0.5*dnorm(x, mean=0, sd=0.5))
}
slice_sampler <- function(start_val, samples){
  x0 <- start_val  # Average of the modes of the individual distribution
  w <- 1
  x <- numeric(samples)
  for(i in 1:samples){ 
    u <- runif(1, 0, pdf(x0))
    u_prime <- runif(1, 0, 1 )
    xmin <- x0 - u_prime*w
    xmax <- xmin + w
    while(u < pdf(xmin)){
      xmin <- xmin - w
    }
    while(u < pdf(xmax)){
      xmax <- xmax + w
    }
    while(TRUE){
      u_prime <- runif(1,0,1)
      x_new <- xmin + u_prime*(xmax - xmin)
      if(u < pdf(x_new)){
        x0 <- x_new
        x[i] <- x0
        break
      }
      if(x_new < x0){
        xmin <- x_new
      }
      else{
        xmax <- x_new
      }
    }
  }
  return(x)
}

reservoirMCMC <- function(reservoir_size=50, L=10, iterates=2000){
  n <- reservoir_size # reservoir size
  l <- L
  m <- iterates # number of iteration
  # matrix to store reservoir at each iteration
  R.m <- matrix(nrow = m, ncol = n)
  # samples <- numerics 
  # matrix to store counting at each iteration
  N.x.m <- matrix(nrow = m, ncol = n)
  
  # step 1
  R.m[1, 1] = rchisq(1, df=1) # start value preferably mean of distribution
  # print(R.m)
  
  
  R.m[1,] <- slice_sampler(R.m[1, 1], n)
  N.x.m[1,] <- rep(1, times = n)
  
  # from iteration 2 to m
  for(i in 2:m){
    
    # step 2
    # generating L more samples from the chain
    y <- slice_sampler(R.m[i-1, n], l+1)[2:(l+1)]
    
    # defining extended extended reservoir and the counting vector
    R.e <- numeric(n+l-1)
    N.e <- numeric(n+l-1)
    
    temp <- n+l-1
    
    R.e[1:n] = R.m[i-1,]
    if(l > 1){
      R.e[(n+1):(n+l-1)] = y[1:l-1]
    }
    N.e[1:n] = N.x.m[i-1,]
    if(l > 1){
      N.e[(n+1):(n+l-1)] = 1
    }
    
    # storing last value of the sample
    y_l = y[l]
    
    # step 3
    # defining initial index vector
    A <- 1:temp
    for(j in 1:l){
      probability <- N.e[A]
      probability <- probability / sum(probability)
      remove_idx <- sample(A, size=1, replace=TRUE, prob=probability)
      # shouldn't we sample without replacement? 
      # removing chosen index
      A <- A[A != remove_idx]
    }
    
    # reservoir and counting vector with n-1 elements
    R.e = R.e[A]
    N.e = N.e[A]
    # adding 1 to count for retained values from previous
    # iteration
    for(j in 1:(n-1)){
      if(A[j] <= n){
        N.e[j] <- N.e[j] + 1
      }
    }
    
    # updating reservoir
    R.m[i,1:n-1] <- R.e
    N.x.m[i,1:n-1] <- N.e
    
    # adding last value to reservoir
    R.m[i,n] <- y_l
    N.x.m[i,n] <- 1
  }
  # print(R.m)
  # print(N.x.m)
  
  x <- c(R.m)
}


##
## Monte Carlo Mean of Autocorrelation
##
B <- 2000
acfarray <- numeric(length = B*10)
acfarray <- matrix(acfarray, nrow=B, ncol=10)
MMcorrel <- numeric(length=3*10)
MMcorrel <- matrix(MMcorrel, nrow=3, ncol=10)
mu_b <- numeric(length = B)
j <- 1

for(l in c(5,10,20)){
  pb <- progress_bar$new(
    format = "  downloading [:bar] :percent eta: :eta",
    total = B, clear = FALSE, width= 60)
  for(i in 1:B){
    x <- reservoirMCMC(L=l)
    samples <- length(x)
    
    temp <- c(acf(x[samples-50:samples], plot=FALSE)[1:10])
    acfarray[i, ] <- temp$acf[,,1]
    pb$tick()
  }
  print(colMeans(acfarray))
  MMcorrel[j,] <- colMeans(acfarray)
  j <- j + 1
}

# Plotting To be done
numlags <- 3
xrange <- range(1:6)
yrange <- range(c(MMcorrel))
plot(xrange,yrange, type='n', ylab = "Autocorrelation", xlab = "Lag")
colors <- rainbow(numlags)
linetype <- c(1: numlags)
plotchar <- seq(18, 18 + numlags, 1)
for(i in 1:numlags){
  lines(1:6, MMcorrel[i,1:6], type='b', lwd=1.5, lty=linetype, col=colors[i], pch = plotchar)
}
title("Monte Carlo mean of Autocorrelation")
legend("topright", legend=c(5,10,20), cex=0.8, col=colors, pch=plotchar,
        lty = linetype, title = "ACF plot")
# plot(MMcorrel[2,], type = 'b', add = TRUE)
# plot(MMcorrel[3,], type = 'b', add = TRUE)


