
# for different use case different implementation of this function
# should be done, but the function definition i.e parameters and return
# value shouldn't change
mh.sampler <- function(start_val, num_of_samples){
  # start_val - first sample of the chain
  # num_of_samples - number of samples to generate
  
  m <- num_of_samples
  sigma <- 4
  x <- numeric(m)
  
  # chain will start from the input start value
  x[1] <- start_val

  u <- runif(m)
  for (i in 2:m) {
    xt <- x[i-1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df = y)
    den <- f(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= num/den) x[i] <- y else {
      x[i] <- xt
    }
  }
  return(x)
}

n <- 10 # sample size
l <- 5
m <- 100 # number of iteration

# matrix to store reservoir at each iteration
R.m <- matrix(nrow = m, ncol = n)

# matrix to store counting at each iteration
N.x.m <- matrix(nrow = m, ncol = n)

# step 1
R.m[1, 1] = rchisq(1, df=1) # start value preferably mean of distribution
# print(R.m)


R.m[1,] <- mh.sampler(R.m[1, 1], n)
N.x.m[1,] <- rep(1, times = n)

# from iteration 2 to m
for(i in 2:m){
  
  # step 2
  # generating L more samples from the chain
  y <- mh.sampler(R.m[i-1, n], l)

  # defining extended extended reservoir and the counting vector
  R.e <- numeric(n+l-1)
  N.e <- numeric(n+l-1)
  
  temp <- n+l-1
  
  R.e[1:n] = R.m[i-1,]
  R.e[(n+1):(n+l-1)] = y[1:l-1]
  N.e[1:n] = N.x.m[i-1,]
  N.e[(n+1):(n+l-1)] = 1
  
  # storing last value of the sample
  y_l = y[l]

  # step 3
  # defining initial index vector
  A <- 1:temp
  for(j in 1:l){
    probability <- N.e[A]
    probability <- probability / sum(probability)
    remove_idx <- sample(A, size=1, replace=TRUE, prob=probability)
    
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

print("Final sample:")
print(R.m[m,])