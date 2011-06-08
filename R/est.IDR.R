est.IDR <-
function(x, mu, sigma, rho, p, eps, max.ite=30){
  
  x1 <- x[,1]
  x2 <- x[,2]
  
  x1.cdf.func <- ecdf(x1)
  x2.cdf.func <- ecdf(x2)
  afactor <- length(x1)/(length(x1)+1)
  x1.cdf <- x1.cdf.func(x1)*afactor
  x2.cdf <- x2.cdf.func(x2)*afactor
  
  # initialization
  para <- list()
  para$mu <- mu
  para$sigma <- sigma
  para$rho <- rho
  para$p <- p  

  j <- 1
  to.run <- TRUE
  loglik.trace <- c()
  loglik.inner.trace <- c()
 
  z.1 <- get.pseudo.mix(x1.cdf, para$mu, para$sigma, para$rho, para$p)
  z.2 <- get.pseudo.mix(x2.cdf, para$mu, para$sigma, para$rho, para$p)

  while(to.run){
    
    # get pseudo value in each cycle

    i <- 1
    while(to.run){
      
      # EM for latent structure
      e.z <- e.step.2normal(z.1, z.2, para$mu, para$sigma, para$rho, para$p)
      para <- m.step.2normal(z.1, z.2, e.z)
 
      if(i > 1)
        l.old <- l.new
    
      # this is just the mixture likelihood of two-component Gaussian
      l.new <- loglik.2binormal(z.1, z.2, para$mu, para$sigma, para$rho, para$p)

      loglik.inner.trace[i] <- l.new 

      if(i > 1){
        to.run <- loglik.inner.trace[i]-loglik.inner.trace[i-1]>eps         
      }
        
#      cat("loglik.inner.trace[", i, "]=", loglik.inner.trace[i], "\n")
#    cat("mu=", para$mu, "sigma=", para$sigma, "p=", para$p, "rho=", para$rho, "\n\n")
      
      i <- i+1
    }
    

    # get pseudo value in each cycle
    z.1 <- get.pseudo.mix(x1.cdf, para$mu, para$sigma, para$rho, para$p)
    z.2 <- get.pseudo.mix(x2.cdf, para$mu, para$sigma, para$rho, para$p)

    if(j > 1)
      l.old.outer <- l.new.outer

    l.new.outer <- loglik.2binormal(z.1, z.2, para$mu, para$sigma, para$rho, para$p)

    loglik.trace[j] <- l.new.outer
    
    if(j == 1)
      to.run <- TRUE
    else{ # stop when iteration>max.ite
      if(j > max.ite)
        to.run <- FALSE
      else
        to.run <- l.new.outer - l.old.outer > eps
    }

#    if(j %% 10==0)
#      cat("loglik.trace[", j, "]=", loglik.trace[j], "\n")
#    cat("mu=", para$mu, "sigma=", para$sigma, "p=", para$p, "rho=", para$rho, "\n")
    
    j <- j+1
  }

#  bic <- -2*l.new + 4*log(length(z.1))

  idr <- 1-e.z
  o <- order(idr)
  idr.o <- idr[o]
  idr.rank <- rank(idr.o, ties.method="max")

  top.mean <- function(index, x){mean(x[1:index])}

  # IDR of selected peaks
  IDR.o <- sapply(idr.rank, top.mean, idr.o)
  IDR <- idr # spaceholder
  IDR[o] <- IDR.o
    
  return(list(para=list(p=para$p, rho=para$rho, mu=para$mu, sigma=para$sigma),
              loglik=l.new, loglik.trace=loglik.trace, idr=1-e.z, IDR=IDR))
}

