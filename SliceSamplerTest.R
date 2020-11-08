samples <- 100000
trudist <- numeric(samples)
for(i in 1:samples){
  p <- runif(1, 0, 1)
  if(p<0.5){
    trudist[i] <- rnorm(1, -1, 1)
  }
  else{
    trudist[i] <- rnorm(1, 0, 0.5)
  }
}
traditional <- hist(trudist, breaks=100, plot=FALSE)

pdf <- function(x){
  return(0.5*dnorm(x, mean=-1, sd=1) + 0.5*dnorm(x, mean=0, sd=0.5))
}
slice_sampler <- function(samples){
  x0 <- -0.5  # Average of the modes of the individual distribution
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

x <- slice_sampler(samples)
chain <- x[1001:samples]
MCMC <- hist(chain, breaks =100, plot=FALSE )
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
plot(traditional, freq = FALSE, col = c1)
plot(MCMC, freq = FALSE, col = c2, add = TRUE)
lines(density(trudist), col='blue')
lines(density(chain), col='red')
