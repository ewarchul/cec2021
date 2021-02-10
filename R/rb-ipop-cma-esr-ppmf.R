library(BBmisc)
library(checkmate)
library(magrittr)
source("cmaesr-utils.R")

#' RB-IPOP-CMAES-PPMF without sigma supress mechanism
#'
#' @description Implementation of PPMF with IPOP rule based on https://ieeexplore.ieee.org/document/7969479
#' and {cmaesr} package  https://github.com/jakobbossek/cmaesr
#' @param last_its_type Type of LastIts restart trigger: mean | best | ave :: String

rb_ipop_cma_esr_ppmf_no_ss = function(
  par,
  fn,
  ...,
  lower = -100, 
  upper = 100,
  control = list()
  ) {
	# extract relevant data
  lb = lower 
  ub = upper
  n = length(par)
	
  # get stopping conditions
  budget = if (n == 10) {
    200000
  } else {
    1000000
  }
  stop_ons_list = 
    c(
      list(stopOnMaxEvals(budget)),
      list(stopOnTolX()),
      list(stopOnCondCov()),
      list(stopOnIndefCovMat()),
      list(stopOnNoEffectAxis()),
      list(stopOnNoEffectCoord()),
      list(stopOnLastIts())
    )
  stop.ons = getCMAESParameter(control, "stop.ons", stop_ons_list)
  if (is.null(stop.ons)) {
    stopf("There must be at least one stopping condition!")
  }
  assertList(stop.ons, min.len = 1L, types = "cma_stopping_condition")
  # alwas check for indefinit covariance matrix first
  stop.ons = c(list(stopOnIndefCovMat()), stop.ons)

  # restart mechanism (IPOP-CMA-ES)
  # RB-IPOP: new restart triggers ->  lastIts and sigSupress 
  restart.triggers = list("conditionCov", "noEffectCoord", "noEffectAxis", "tolX", "indefCovMat", "lastIts")

  stop.ons.names = sapply(stop.ons, function(stop.on) stop.on$code)
  if (!isSubset(restart.triggers, stop.ons.names)) {
    stopf("Only codes / short names of active stopping conditions allowed as restart trigger, but '%s' are no stopping conditions.", collapse(setdiff(restart.triggers, stop.ons.names), sep = ", "))
  }
  restart.multiplier = getCMAESParameter(control, "restart.multiplier", 2)
  assertNumber(restart.multiplier, lower = 1, na.ok = FALSE, finite = TRUE)
  max.restarts = getCMAESParameter(control, "max.restarts", 100)
  assertInt(max.restarts)

  #FIXME: default value should be derived from bounds
  sigma = getCMAESParameter(control, "sigma", 7)
  assertNumber(sigma, lower = 0L, finite = TRUE)
  d_param = getCMAESParameter(control, "d_param", 2)
  p_target = getCMAESParameter(control, "p_target", 0.1)
  last_its_type = getCMAESParameter(control, "last_its_type", "mean")

  # Precompute E||N(0,I)||
	chi.n = sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n^2))

	# bookkeep best individual
	best.param = rep(NA, n)
	best.fitness = Inf

  # set initial distribution mean

  # logs
  log.population = getCMAESParameter(control, "log.population", FALSE)
  assertFlag(log.population, na.ok = FALSE)
  population.trace = list()

  # init some termination criteria stuff
	iter = 0L
  n.evals = 0L
	start.time = Sys.time()

  bestVal.log = matrix(0, nrow=0, ncol=1)

  # somehow dirty trick to "really quit" if stopping condition is met and
  # now more restart should be triggered.
  do.terminate = FALSE

  max_dx = (ub - lb) / 5
  last_its = 10 + ceiling(30 * n / (4 * n))
  print(sigma)
  print(d_param)
  print(p_target)

  for (run in 0:max.restarts) {
    # population and offspring size
    if (run == 0) {
      lambda = getCMAESParameter(control, "lambda", 4 * n)
      assertInt(lambda, lower = 4)
      mu = getCMAESParameter(control, "mu", floor(lambda / 2))
      assertInt(mu)
	    m = par
    } else {
      lambda = getCMAESParameter(control, "lambda", 4 * n)
      # increase population size (IPOP-CMA-ES)
      lambda = ceiling(restart.multiplier^run * lambda)
      mu = floor(lambda / 2)
      m = runif(n, lb, ub)
    }

    # path for covariance matrix C and stepsize sigma
    pc = rep(0, n)
    ps = rep(0, n)

    # initialize recombination weights
    weights = getCMAESParameter(control, "weights", log(mu + 1) - log(1:mu))
    if (any(weights < 0)) {
      stopf("All weights need to be positive, but there are %i negative ones.", sum(which(weights < 0)))
    }
    if (length(weights) != mu) {
      stopf("You need to pass %i 'weights', but you passed %i.", mu, length(weights))
    }

    # normalize weights
    weights = weights / sum(weights)

    # variance-effectiveness / variance effective selection mass of sum w_i x_i
    mu.eff = sum(weights)^2 / sum(weights^2) # chosen such that mu.eff ~ lambda/4

    # step-size control
    cs = (mu.eff + 2) / (n + mu.eff + 3)
    ds = 1 + 2 * max(0, sqrt((mu.eff - 1) / (n + 1)) - 1) + cs # damping factor
    # covariance matrix Adaptation parameters
    cc = 4 / (n + 4)
    alpha.mu = 2L
    cmu = mu.eff 
    ccov = (1/cmu) * 2/(n+1.4)^2 + (1-1/cmu) * ((2*cmu-1)/((n+2)^2+2*cmu))

    # covariance matrix
    sigma = getCMAESParameter(control, "sigma", 7)
    B = diag(n)
    D = diag(n)
    BD = B %*% D
    C = BD %*% t(BD) # C = B D^2 B^T = B B^T, since D equals I_n
    Cinvsqrt = B %*% diag(1 / sqrt(diag(D))) %*% t(B)

    # no restart trigger fired until now
    restarting = FALSE

    # RB-IPOP variables 
    refPointFit = Inf
    sigma.repair.cnt = 0

    # PPMF variables 
    eval_mean = Inf
    eval_meanOld = Inf

    # break inner loop if terminating stopping condition active or
    # restart triggered
  	while (!restarting) {
      iter = iter + 1L

      # create new population of search points
  		arz = matrix(rnorm(n * lambda), ncol = lambda) # ~ N(0, I)
      ary = BD %*% arz # ~ N(0, C)
      arx = m + sigma * ary # ~ N(m, sigma^2 C)


      # Here we apply a penalization of violated bounds
      arx.repaired = ifelse(arx < lb, lb, ifelse(arx > ub, ub, arx))

      # Prepare penalization based on distance to repaired points (see Eq. 51)
      penalty.alpha = 1L
      penalties = penalty.alpha + colSums((arx - arx.repaired)^2)
      penalties[is.infinite(penalties)] = .Machine$double.max / 2

      # compute fitness values of repaired points
      fitn.repaired = 
        apply(arx.repaired, 2, function(x) fn(x))

      # apply penalization (see Eq. 51)
      fitn = fitn.repaired * penalties
    
      # update evaluation
  		n.evals = n.evals + lambda

      # order fitness values
      fitn.ordered.idx = order(fitn, decreasing = FALSE)
      fitn.ordered = fitn[fitn.ordered.idx]

      bestVal.log = rbind(bestVal.log, min(suppressWarnings(min(bestVal.log)), min(fitn.ordered)))

      # update best solution so far
      valid = (penalties <= 1)
      if (any(valid)) {
        #stopifnot(all(fitn[valid] == fitn.repaired[valid]))
        #stopifnot(all(arx[, valid, drop = FALSE] == arx.repaired[, valid, drop = FALSE]))
        min.valid.idx = which.min(fitn.repaired[valid])
        if (fitn.repaired[valid][min.valid.idx] < best.fitness) {
          best.fitness = fitn.repaired[valid][min.valid.idx]
          best.param = arx[, valid, drop = FALSE][, min.valid.idx]
        }
      }

      # update mean value / center of mass
      new.pop.idx = fitn.ordered.idx[1:mu]
      x.best = arx[, new.pop.idx, drop = FALSE]
      m.old = m
      m = drop(x.best %*% weights)

      
      y.best = ary[, new.pop.idx, drop = FALSE]
      y.w = drop(y.best %*% weights)
      z.best = arz[, new.pop.idx, drop = FALSE]
      z.w = drop(z.best %*% weights)

      # log population
      if (log.population) {
        population.trace[[iter]] = x.best
      }

  		# Update evolution path with cumulative step-size Adaptation (CSA) / path length control
      # For an explanation of the last factor see appendix A in https://www.lri.fr/~hansen/cmatutorial.pdf
      ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mu.eff) * (B %*% z.w)
      h.sigma <- drop((norm2(ps)/sqrt(1-(1-cs)^(2*n.evals/lambda))/chi.n) < (1.4 + 2/(n+1)))
  		# Update covariance matrix
      pc = (1 - cc) * pc + h.sigma * sqrt(cc * (2 - cc) * mu.eff) * y.w
      y = BD %*% z.best

      C = (1 - ccov) * C +
        ccov * (1/cmu) * (pc %o% pc + (1-h.sigma) * cc*(2-cc) * C) +
        ccov * (1-1/cmu) * y %*% diag(weights) %*% t(y)

      eval_meanOld = eval_mean
      mean_point = apply(arx.repaired, 1, mean) %>% t() %>% t()
      eval_mean = apply(mean_point, 2, function(x) fn(x))
      n.evals = n.evals + 1

      if (iter %% last_its == 0) {
        oldRefPointFit = refPointFit
        if (last_its_type == "mean") {
          refPointFit = fn(m)
          n.evals = n.evals + 1
        } else if (last_its_type == "ave") {
          refPointFit = eval_mean
        } else if (last_its_type == "best") {
          refPointFit = best.fitness
        } else {
          stop("last_its_type = mean | ave | best")
        }
      }

      # Update step-size sigma
      p_succ = 
        length(which(fitn.ordered < eval_meanOld))/lambda
      sigma = 
        sigma * exp(d_param * (p_succ - p_target) / (1 - p_target))
      
      # Finally do decomposition C = B D^2 B^T
      e = eigen(C, symmetric = TRUE)
      B = e$vectors
      D = diag(sqrt(e$values), length(e$values))
      BD = B %*% D
      Cinvsqrt = B %*% diag(1 / diag(D)) %*% t(B) # update C^-1/2

      # escape flat fitness values
      # RB-IPOP: Increase sigma multiplier by 2
      if (fitn.ordered[1] == fitn.ordered[min(1+floor(lambda/2), 2+ceiling(lambda/4))]) {
        sigma = 20 * sigma * exp(0.2 + cs / ds)
      }

      # CHECK STOPPING CONDITIONS
      # =========================
      stop.obj = checkStoppingConditions(stop.ons)

      n.stop.codes = length(stop.obj$codes)
      if (max.restarts > 0L && any(stop.obj$codes %in% restart.triggers)) {
        n.stop.codes = sum(!(stop.obj$codes %in% restart.triggers))
        restarting = TRUE
      }

      # check if CMA-ES should really quit, i.e., is there a stopping condition,
      # that is active and does not trigger a restart?
      if (!restarting && (n.stop.codes > 0L)) {
        do.terminate = TRUE
        break
      }
  	}

    # really quit without more restarts
    if (do.terminate) {
      break
    }
  }
  
  log = list()
  log$bestVal = bestVal.log

	return(list(
      best.param = best.param,
      best.fitness = best.fitness,
      n.evals = n.evals,
      past.time = as.integer(difftime(Sys.time(), start.time, units = "secs")),
      n.iters = iter - 1L,
      n.restarts = run,
      label = "rb_ipop_cma_esr_ppmf_no_ss",
      population.trace = population.trace,
      diagnostic = log,
      message = stop.obj$stop.msgs,
      classes = "cma_result"
	  )
  )
}

