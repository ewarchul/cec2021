cma_es_csa <- function(par, fn, ..., lower, upper, if_CMA = TRUE, control=list()) {
  norm <- function(x)
    drop(sqrt(crossprod(x)))
  
  controlParam <- function(name, default) {
    v <- control[[name]]
    if (is.null(v))
      return (default)
    else
      return (v)
  }
  
  ## Inital solution:
  xmean <- par
  N <- length(xmean)
  ## Box constraints:
  if (missing(lower))
    lower <- rep(-Inf, N)
  else if (length(lower) == 1)  
    lower <- rep(lower, N)
  
  if (missing(upper))
    upper <- rep(Inf, N)
  else if (length(upper) == 1)  
    upper <- rep(upper, N)
  
  ## Parameters:
  trace       <- controlParam("trace", FALSE)
  fnscale     <- controlParam("fnscale", 1)
  stopfitness <- controlParam("stopfitness", -Inf)
  budget      <- controlParam("budget", 200000)                     ## The maximum number of fitness function calls
  sigma       <- controlParam("sigma", 1)
  sc_tolx     <- controlParam("stop.tolx", 1e-12 * sigma) ## Undocumented stop criterion
  keep.best   <- controlParam("keep.best", TRUE)
  vectorized  <- controlParam("vectorized", FALSE)
  
  ## Logging options:
  log.all    <- controlParam("diag", TRUE)
  log.sigma  <- controlParam("diag.sigma", log.all)
  log.budget  <- controlParam("diag.budget", log.all)
  log.eigen  <- controlParam("diag.eigen", log.all)
  log.value  <- controlParam("diag.value", log.all)
  log.pop    <- controlParam("diag.pop", log.all)
  log.bestVal<- controlParam("diag.bestVal", log.all)
  
  ## Strategy parameter setting (defaults as recommended by Nicolas Hansen):
  lambda      <- controlParam("lambda", 4 * N)
  maxiter     <- controlParam("maxit", round(budget/lambda))
  mu          <- controlParam("mu", floor(lambda/2))
  weights     <- controlParam("weights", log(mu+1) - log(1:mu))
  weights     <- weights/sum(weights)
  mueff       <- controlParam("mueff", sum(weights)^2/sum(weights^2))
  cc          <- controlParam("ccum", 4/(N+4))
  cs          <- controlParam("cs", (mueff+2)/(N+mueff+3))
  mucov       <- controlParam("ccov.mu", mueff)
  ccov        <- controlParam("ccov.1",
                              (1/mucov) * 2/(N+1.4)^2
                              + (1-1/mucov) * ((2*mucov-1)/((N+2)^2+2*mucov)))
  damps       <- controlParam("damps",
                              1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs)
  
  ## Safety checks:
  stopifnot(length(upper) == N)  
  stopifnot(length(lower) == N)
  stopifnot(all(lower < upper))
  stopifnot(length(sigma) == 1)
  
  ## Bookkeeping variables for the best solution found so far:
  best.fit <- Inf
  best.par <- NULL
  
  ## Preallocate logging structures:
  if (log.sigma)
    sigma.log <- numeric(maxiter)
  if (log.budget)
    budget.log <- numeric(maxiter)
  if (log.eigen)
    eigen.log <- matrix(0, nrow=maxiter, ncol=N)
  if (log.value)
    value.log <- matrix(0, nrow=maxiter, ncol=mu)
  if (log.pop)
    pop.log <- array(0, c(N, mu, maxiter))
  if(log.bestVal)
    bestVal.log <-  matrix(0, nrow=0, ncol=1)
  
  ## Initialize dynamic (internal) strategy parameters and constants
  pc <- rep(0.0, N)
  ps <- rep(0.0, N)
  B <- diag(N)
  D <- diag(N)
  BD <- B %*% D
  C <- BD %*% t(BD)
  
  chiN <- sqrt(N) * (1-1/(4*N)+1/(21*N^2))
  
  iter <- 0L      ## Number of iterations
  counteval <- 0L ## Number of function evaluations
  cviol <- 0L     ## Number of constraint violations
  msg <- NULL     ## Reason for terminating
  nm <- names(par) ## Names of parameters

  
  ## Preallocate work arrays:
  arx <- matrix(0.0, nrow=N, ncol=lambda)
  #arx <-  replicate(lambda, runif(N,lower,upper))
  arfitness <- rep(Inf, lambda)
  #arfitness <- apply(arx, 2, function(x) fn(x, ...) * fnscale)
  while (counteval < budget) {
    iter <- iter + 1L
    
    if (!keep.best) {
      best.fit <- Inf
      best.par <- NULL
    }
    if (log.sigma)
      sigma.log[iter] <- sigma
    
    ## Generate new population:
    arz <- matrix(rnorm(N*lambda), ncol=lambda)
    arx <- xmean + sigma * (BD %*% arz)

    vx <- ifelse(arx > lower, ifelse(arx < upper, arx, upper), lower)

    if (!is.null(nm))
      rownames(vx) <- nm
    pen <- 1 + colSums((arx - vx)^2)
    pen[!is.finite(pen)] <- .Machine$double.xmax / 2
    cviol <- cviol + sum(pen > 1)

    if (vectorized) {
      y <- fn(vx, ...) * fnscale
    } else {
      y <- apply(vx, 2, function(x) fn(x, ...) * fnscale)
    }
    counteval <- counteval + lambda

    arfitness <- y * pen

    valid <- pen <= 1
    if (any(valid)) {
      wb <- which.min(y[valid])
      if (y[valid][wb] < best.fit) {
        best.fit <- y[valid][wb]
        best.par <- arx[,valid,drop=FALSE][,wb]
      }
    }
    
    ## Order fitness:
    arindex <- order(arfitness)
    arfitness <- arfitness[arindex]

    if (log.bestVal) 
      bestVal.log <- rbind(bestVal.log,min(suppressWarnings(min(bestVal.log)), min(arfitness)))

    
    aripop <- arindex[1:mu]
    selx <- arx[,aripop]
    xmean <- drop(selx %*% weights)
    selz <- arz[,aripop]
    zmean <- drop(selz %*% weights)
    
    ## Save selected x value:
    if (log.pop) pop.log[,,iter] <- selx
    if (log.value) value.log[iter,] <- arfitness[aripop]
    
    ## Cumulation: Update evolutionary paths
    ps <- (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * (B %*% zmean)
    hsig <- drop((norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN) < (1.4 + 2/(N+1)))
    pc <- (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * drop(BD %*% zmean)

    ## Adapt Covariance Matrix:
    BDz <- BD %*% selz

    if (if_CMA) 
      C <- (1-ccov) * C + ccov * (1/mucov) *
        (pc %o% pc + (1-hsig) * cc*(2-cc) * C) +
        ccov * (1-1/mucov) * BDz %*% diag(weights) %*% t(BDz)
    else
      C = C
    
    ## Adapt step size sigma:
    sigma <- sigma * exp((norm(ps)/chiN - 1)*cs/damps)

    
    e <- eigen(C, symmetric=TRUE)
    eE <- eigen(cov(t(arx)))
    if (log.eigen)
      eigen.log[iter,] <- rev(sort(eE$values))
    
    if (!all(e$values >= sqrt(.Machine$double.eps) * abs(e$values[1]))) {      
      msg <- "Covariance matrix 'C' is numerically not positive definite."
      #break
    }

    if (log.budget)
      budget.log[iter] <- counteval
    
    B <- e$vectors
    D <- diag(sqrt(e$values), length(e$values))
    BD <- B %*% D
    
    ## break if fit:
    if (arfitness[1] <= stopfitness * fnscale) {
      msg <- "Stop fitness reached."
      B <- diag(N)
      D <- diag(N)
      BD <- B %*% D
      C <- BD %*% t(BD)
   #   break
    }
    
    ## Check stop conditions:
    
    ## Condition 1 (sd < tolx in all directions):
#    if (all(D < sc_tolx) && all(sigma * pc < sc_tolx)) {
#      msg <- "All standard deviations smaller than tolerance."
#      break
#    }
    
    ## Escape from flat-land:
    if (arfitness[1] == arfitness[min(1+floor(lambda/2), 2+ceiling(lambda/4))]) { 
      sigma <- sigma * exp(0.2+cs/damps);
    if (trace)
      message("Flat fitness function. Increasing sigma.")
    }

    if (trace)
      message(sprintf("Iteration %i of %i: current fitness %f",
                      iter, maxiter, arfitness[1] * fnscale))
  }
  cnt <- c(`function`=as.integer(counteval), gradient=NA)
  
  log <- list()
  ## Subset lognostic data to only include those iterations which
  ## where actually performed.
  if (log.value) log$value <- value.log[1:iter,]
  if (log.sigma) log$sigma <- sigma.log[1:iter]
  if (log.budget) log$budget <- budget.log[1:iter]
  if (log.eigen) log$eigen <- eigen.log[1:iter,]
  if (log.pop)   log$pop   <- pop.log[,,1:iter]
  if (log.bestVal) log$bestVal <- bestVal.log
  
  ## Drop names from value object
  names(best.fit) <- NULL
  res <- list(par=best.par,
              value=best.fit / fnscale,
              counts=cnt,
              iter = iter,
              convergence=ifelse(iter >= maxiter, 1L, 0L),
              message=msg,
              label = "cma-es-csa",
              constr.violations=cviol,
              diagnostic=log
  )
  class(res) <- "cma_es_csa.result"
  return(res)
}
