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

#' RB-IPOP restart mechanism (sigma repair)

stopOnSigSupress = function(max_repairs = 5) {
  assertInt(max_repairs, na.ok = FALSE)
  force(max_repairs)
  return(makeStoppingCondition(
    name = "sigSupress",
    message = sprintf("Sigma was repaired more than %f times", max_repairs),
    stop.fun = function(envir = parent.frame()) {
      return(envir$sigma.repair.cnt > max_repairs)
    } 
  ))
}

#' RB-IPOP restart mechanism

stopOnLastIts = function(treshold = 1e-8) {
  assertInt(treshold, na.ok = FALSE)
  force(treshold)
  return(makeStoppingCondition(
    name = "lastIts",
    message = sprintf("The range of the reference point objective function values of lastIts iteration is zero."),
    stop.fun = function(envir = parent.frame()) {
      ref_delta = Inf
      if (envir$iter %% envir$last_its == 0) {
        ref_delta = abs(envir$oldRefPointFit - envir$refPointFit)
      } 
      return(ref_delta < treshold)
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
      return(sum((m - (m + 0.3 * envir$sigma * lambdai * ui))^2) < .Machine$double.eps)
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
