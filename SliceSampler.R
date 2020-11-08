samples <- 100000
pdf <- function(x){
  return(0.5*dnorm(x, mean=-1, sd=1) + 0.5*dnorm(x, mean=0, sd=0.5))
}
# Reference 1: https://projecteuclid.org/download/pdf_1/euclid.aos/1056562461 (Neal(2003))
# Reference 2: https://www.ias.informatik.tu-darmstadt.de/uploads/Teaching/RobotLearningSeminar/Dittmar_RLS_2013.pdf
# Code : http://www.cs.toronto.edu/~radford/ftp/slice-R-prog
slice_sampler <- function(samples){
  x0 <- -0.5  # Average of the modes of the individual distributions
  w <- 1 # Estimate of step size (Static)
  x <- numeric(samples) # Return these samples
  for(i in 1:samples){ 
    u <- runif(1, 0, pdf(x0))
    u_prime <- runif(1, 0, 1 )
    xmin <- x0 - u_prime*w
    xmax <- xmin + w
    # Expand the interval until its ends are outside the slice, or until
    # the limit on steps is reached.
    while(u < pdf(xmin)){
      xmin <- xmin - w
    }
    while(u < pdf(xmax)){
      xmax <- xmax + w
    } 
    while(TRUE){
      u_prime <- runif(1,0,1)
      x_new <- xmin + u_prime*(xmax - xmin) # Draw a new sample from the slice
      if(u < pdf(x_new)){ # Check the threshold
        x0 <- x_new
        x[i] <- x0 # Take this sample
        break
      }
      if(x_new < x0){ 
        xmin <- x_new # Adjust the lower bound
      }
      else{
        xmax <- x_new # Adjust the upper bound
      }
    }
  }
  return(x)
}
# Create Burn-in samples
burnin <- 1000
x <- slice_sampler(samples)
chain <- x[burnin + 1:samples]
hist(chain, breaks =100, freq=FALSE )
lines(density(chain), col='red')