#Aim to obtain representive sample of complete network.
data <- read.csv("~/Desktop/Comp Stats Codes/Computational-Stats-Project/HU_edges.csv", header = TRUE)
el <- as.matrix(data) 
el[,1] <- as.character(el[,1])
el[,2] <- as.character(el[,2])
library(igraph)
g <- graph.edgelist(el, directed = FALSE)
# g <- graph.ring(11)

# change the following for graph 
MHRW <- function(start_node = 3,numsamples=50000){
   #Take samples c(50, 100)*2000
  x = rep(0,numsamples)
  # sample <- numeric(numsamples)
  
  x[1] <- start_node 
  # sample[1] <- degree(g, x[1])
  # print(sample[1])
  #initialize; I've set arbitrarily set this to 3
  pb <- progress_bar$new(
    format = "  completed [:bar] :percent eta: :eta",
    total = numsamples, clear = TRUE, width= 60)
  # Metropolis Hastings Random Walk
  rejected <- 0
  for(i in 2:numsamples){
    # progress(i, progress.bar = TRUE)
    current_x <- x[i-1]
    neighbors <- neighbors(g, current_x)
    p <- runif(1)
    proposed <- sample(neighbors, size=1)
  
    k_w <- degree(g, proposed)
    k_v <- degree(g, current_x)
    if((k_v/k_w) >= p){
      x[i] = proposed  
      # sample[i] <- degree(g, x[i])
    } else {
      rejected <- rejected + 1
      x[i] = current_x  
      # sample[i] <- sample[i-1]
    }
    pb$tick()
  }
  cat("% Rejected: ", rejected/numsamples)
  return(x)
}

reservoirMCMC <- function(reservoir_size=50, L=10, iterates=2000){
  reservoir_size=50
  L=10
  iterates=100
   n <- reservoir_size # reservoir size
  l <- L
  m <- iterates # number of iteration
  # matrix to store reservoir at each iteration
  R.m <- matrix(nrow = m, ncol = n)
  # samples <- numerics 
  # matrix to store counting at each iteration
  N.x.m <- matrix(nrow = m, ncol = n)
  
  # step 1
  R.m[1, 1] = 2 # start value preferably mean of distribution
  # print(R.m)
  
  
  R.m[1,] <- MHRW(R.m[1, 1], n)
  N.x.m[1,] <- rep(1, times = n)
  
  # from iteration 2 to m
  for(i in 2:m){
    
    # step 2
    # generating L more samples from the chain
    y <- MHRW(R.m[i-1, n], l+1)[2:(l+1)]

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

numsamples < 5000
sample <- MHRW(numsamples = 5000)
res_samples <- reservoirMCMC(reservoir_size = 50, L = 10, iterates = 100)
q1 <- 7
q2 <- 15
q3 <- 25
p1 <- which(degree(g, sample) < q1)

prop<- length(p1)/numsamples # Calculating the proportions given in the paper
hist(degree(g,sample),xlim=range(degree(g, sample)),probability = TRUE, main="Histogram of values of x visited by MH algorithm")
# lines(xx,target(xx),col="red")

hist(degree(g, res_samples))
p1_ <- which(degree(g, res_samples) < q1)
prop_ <- length(p1_)/numsamples
print(prop)
print(prop_)
