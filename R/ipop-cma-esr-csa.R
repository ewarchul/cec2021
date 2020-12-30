library(BBmisc)
library(checkmate)

ipop_cma_esr_csa = function(
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
      list(stopOnTolX()),
      list(stopOnCondCov()),
      list(stopOnIndefCovMat()),
      list(stopOnNoEffectAxis()),
      list(stopOnNoEffectCoord()),
      list(stopOnMaxEvals(budget))
    )
  stop.ons = getCMAESParameter(control, "stop.ons", stop_ons_list)
  if (is.null(stop.ons)) {
    stopf("There must be at least one stopping condition!")
  }
  assertList(stop.ons, min.len = 1L, types = "cma_stopping_condition")
  # alwas check for indefinit covariance matrix first
  stop.ons = c(list(stopOnIndefCovMat()), stop.ons)
  # restart mechanism (IPOP-CMA-ES)
  restart.triggers = list("conditionCov", "noEffectCoord", "noEffectAxis", "tolX", "indefCovMat")
#  restart.triggers = getCMAESParameter(control, "restart.triggers", character(0L))
  stop.ons.names = sapply(stop.ons, function(stop.on) stop.on$code)
  if (!isSubset(restart.triggers, stop.ons.names)) {
    stopf("Only codes / short names of active stopping conditions allowed as restart trigger, but '%s' are no stopping conditions.", collapse(setdiff(restart.triggers, stop.ons.names), sep = ", "))
  }
  restart.multiplier = getCMAESParameter(control, "restart.multiplier", 2)
  assertNumber(restart.multiplier, lower = 1, na.ok = FALSE, finite = TRUE)
  max.restarts = getCMAESParameter(control, "max.restarts", 100)
  assertInt(max.restarts)

  #FIXME: default value should be derived from bounds
  sigma = getCMAESParameter(control, "sigma", 1)
  assertNumber(sigma, lower = 0L, finite = TRUE)

  # Precompute E||N(0,I)||
	chi.n = sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n^2))

	# bookkeep best individual
	best.param = rep(NA, n)
	best.fitness = Inf

  # set initial distribution mean
	m = par

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

  for (run in 0:max.restarts) {
    # population and offspring size
    if (run == 0) {
      lambda = getCMAESParameter(control, "lambda", 4 * n)
      assertInt(lambda, lower = 4)
      mu = getCMAESParameter(control, "mu", floor(lambda / 2))
      assertInt(mu)
    } else {
      lambda = getCMAESParameter(control, "lambda", 4 * n)
      # increase population size (IPOP-CMA-ES)
      lambda = ceiling(restart.multiplier^run * lambda)
      mu = floor(lambda / 2)
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
    sigma = getCMAESParameter(control, "sigma", 1)
    B = diag(n)
    D = diag(n)
    BD = B %*% D
    C = BD %*% t(BD) # C = B D^2 B^T = B B^T, since D equals I_n
    Cinvsqrt = B %*% diag(1 / sqrt(diag(D))) %*% t(B)

    # no restart trigger fired until now
    restarting = FALSE

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

      # Update step-size sigma
      sigma <- sigma * exp((norm2(ps)/chi.n - 1)*cs/ds)

      
      # Finally do decomposition C = B D^2 B^T
      e = eigen(C, symmetric = TRUE)
      B = e$vectors
      D = diag(sqrt(e$values), length(e$values))
      BD = B %*% D
      Cinvsqrt = B %*% diag(1 / diag(D)) %*% t(B) # update C^-1/2

      # escape flat fitness values
      if (fitn.ordered[1] == fitn.ordered[min(1+floor(lambda/2), 2+ceiling(lambda/4))]) {
        sigma = sigma * exp(0.2 + cs / ds)
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
      label = "ipop_cma_esr_csa",
      population.trace = population.trace,
      diagnostic = log,
      message = stop.obj$stop.msgs,
      classes = "cma_result"
	  )
  )
}

makeStoppingCondition = function(name, message, stop.fun, code = name, control = list()) {
  assertString(name, na.ok = FALSE)
  assertString(message, na.ok = FALSE)
  assertFunction(stop.fun, args = "envir")
  assertString(code, na.ok = FALSE)
  assertList(control)
  makeS3Obj(
    name = name,
    message = message,
    stop.fun = stop.fun,
    code = code,
    control = control,
    classes = "cma_stopping_condition"
  )
}

shouldStop = function(x, envir) {
  UseMethod("shouldStop")
}

shouldStop.cma_stopping_condition = function(x, envir) {
  return(x$stop.fun(envir))
}

checkStoppingConditions = function(stop.ons, envir = parent.frame()) {
  assertList(stop.ons, min.len = 1L, types = "cma_stopping_condition")
  stop.msgs = character(0L)
  codes = character(0L)
  for (stop.on in stop.ons) {
    if (shouldStop(stop.on, envir = envir)) {
      stop.msgs = c(stop.msgs, stop.on$message)
      codes = c(codes, stop.on$code)
      # since some stopping conditions need a "correct" covariance matrix
      # we stop here if the first condition is met (infefCovMat is first)
      break
    }
  }
  return(list(stop.msgs = stop.msgs, codes = codes))
}

getDefaultStoppingConditions = function() {
  return(
    list(
      stopOnTolX(),
      stopOnNoEffectAxis(),
      stopOnNoEffectCoord(),
      stopOnCondCov()
    )
  )
}

stopOnMaxIters = function(max.iter = 100L) {
  assertInt(max.iter, na.ok = FALSE)
  force(max.iter)
  return(makeStoppingCondition(
    name = "maxIter",
    message = sprintf("MaxIter: reached maximal number of iterations/generations %i.", max.iter),
    stop.fun = function(envir = parent.frame()) {
      return(envir$iter > max.iter)
    }
  ))
}

# @title Stopping condition: indefinite covariance matrix.
#
# @description Stop if covariance matrix is not positive definite anymore.
#
# @return [\code{cma_stopping_condition}]
# @family stopping conditions
#NOTE: this one is not exported. However, it is always prepended to the list of
# stopping conditions, since
stopOnIndefCovMat = function() {
  return(makeStoppingCondition(
    name = "indefCovMat",
    message = "Covariance matrix is not numerically positive definite.",
    stop.fun = function(envir = parent.frame()) {
      e.values = envir$e$values
      return(any(is.na(e.values)) || any(e.values <= sqrt(.Machine$double.eps) * abs(e.values[1L])))
    }
  ))
}

#' @title Stopping condition: optimal params.
#'
#' @description Stop if euclidean distance of parameter is below
#' some tolerance value.
#'
#' @param opt.param [\code{numeric}]\cr
#'   Known optimal parameter settings.
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \eqn{1e^{-8}}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnOptParam = function(opt.param, tol = 1e-8) {
  assertNumeric(opt.param, any.missing = FALSE, all.missing = FALSE)
  assertNumber(tol, lower = 0, na.ok = FALSE, finite = TRUE)
  force(opt.param)
  force(tol)
  return(makeStoppingCondition(
    name = "optParamTol",
    message = sprintf("Optimal parameters approximated nicely (gap < %.2f).", tol),
    stop.fun = function(envir = parent.frame()) {
      return(sqrt(sum(envir$best.param - opt.param)^2) < tol)
    }
  ))
}

#' @title Stopping condition: optimal objective value.
#'
#' @description Stop if best solution is close to optimal objective value.
#'
#' @param opt.value [\code{numeric(1)}]\cr
#'   Known optimal objective function value.
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \eqn{1e^{-8}}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnOptValue = function(opt.value, tol = 1e-8) {
  assertNumber(opt.value, na.ok = FALSE)
  assertNumber(tol, lower = 0, na.ok = FALSE, finite = TRUE)
  force(opt.value)
  force(tol)
  return(makeStoppingCondition(
    name = "optValTol",
    message = sprintf("Optimal function value approximated nicely (gap < %.10f).", tol),
    stop.fun = function(envir = parent.frame()) {
      return(abs(envir$best.fitness - opt.value) < tol)
    }
  ))
}


stopOnTimeBudget = function(budget) {
  assertInt(budget, na.ok = FALSE, lower = 1L)
  force(budget)
  return(makeStoppingCondition(
    name = "timeBudget",
    message = sprintf("Time budget of %i [secs] reached.", budget),
    stop.fun = function(envir = parent.frame()) {
      return(difftime(Sys.time(), envir$start.time, units = "secs") > budget)
    }
  ))
}

stopOnMaxEvals = function(max.evals) {
  assertInt(max.evals, na.ok = FALSE, lower = 1L)
  force(max.evals)
  return(makeStoppingCondition(
    name = "maxEvals",
    message = sprintf("Maximal number of %i function evaluations reached.", max.evals),
    stop.fun = function(envir = parent.frame()) {
      return(envir$n.evals >= max.evals)
    }
  ))
}


stopOnTolX = function(tol = 1e-12) {
  assertInt(tol, na.ok = FALSE)
  force(tol)
  return(makeStoppingCondition(
    name = "tolX",
    message = sprintf("Standard deviation below tolerance in all coordinates."),
    stop.fun = function(envir = parent.frame()) {
      return(all(envir$D < tol) && all((envir$sigma * envir$p.c) < tol))
    }
  ))
}


stopOnNoEffectAxis = function() {
  return(makeStoppingCondition(
    name = "noEffectAxis",
    message = "Addition of 0.1 times sigma does not change mean value.",
    stop.fun = function(envir = parent.frame()) {
      ii = (envir$iter %% envir$n) + 1L
      ui = envir$e$vectors[, ii]
      lambdai = sqrt(envir$e$values[ii])
      m = envir$m
      return(sum((m - (m + 0.1 * envir$sigma * lambdai * ui))^2) < .Machine$double.eps)
    }
  ))
}


stopOnNoEffectCoord = function() {
  return(makeStoppingCondition(
    name = "noEffectCoord",
    message = "Addition of 0.2 times sigma in any coordinate does not change mean value.",
    stop.fun = function(envir = parent.frame()) {
      m = envir$m
      return(sum((m - (m + 0.2 * envir$sigma))^2) < .Machine$double.eps)
    }
  ))
}


stopOnCondCov = function(tol = 1e14) {
  assertNumber(tol, na.ok = FALSE, lower = 0, finite = TRUE)
  force(tol)
  return(makeStoppingCondition(
    name = "conditionCov",
    message = sprintf("Condition number of covariance matrix exceeds %f", tol),
    stop.fun = function(envir = parent.frame()) {
      return(kappa(envir$C) > tol)
    }
  ))
}

getCMAESParameter = function(control, what, default) {
  return(BBmisc::coalesce(control[[what]], default))
}

norm2 = function(x) {
  return(drop(sqrt(crossprod(x))))
}
