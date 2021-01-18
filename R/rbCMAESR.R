library(BBmisc)
library(magrittr)
library(checkmate)
source("bounds.R")



myCMAESR <- function(objFun, lower, upper, initSol, repairF, MaxFES, stateDumpAtEv, resurrState = NULL, genAddSolsFun = NULL) {
  nDim <- length(lower)
  SCREEN_WIDTH <- 78
  mod2putProgr <- floor(MaxFES / SCREEN_WIDTH)


  if (length(stateDumpAtEv) > 0) {
    bestSoFarSolHist <- matrix(0, nDim, length(stateDumpAtEv)) # najlepszy osobnik
    bestSoFarFitHist <- vector("numeric", length(stateDumpAtEv)) # najlepsze fitness
    whenFoundBestSolHist <- vector("numeric", length(stateDumpAtEv))
    if (!is.null(genAddSolsFun)) {
      bestSoFarAdditSolsHist <- array(0, c(nDim, ADD_SOLS_SIZE, length(stateDumpAtEv)))
      bestSoFarAdditSolsFitHist <- matrix(0, ADD_SOLS_SIZE, length(stateDumpAtEv))
      whenFoundAdditBestSolsHist <- matrix(0, ADD_SOLS_SIZE, length(stateDumpAtEv))
    }
    cecDumpIndx <- 1
  }


  objective.fun <- objFun
  start.point <- initSol
  stop.ons <- c(list(stopOnMaxEvals(MaxFES)), list(stopOnOptValue(cec17optimumValues[funNmbr], STOP_FIT)), list(stopOnIndefCovMat()), getDefaultStoppingConditions())

  # extract relevant data
  lb <- lower
  ub <- upper
  n <- nDim

  control <- list(
    restart.triggers = list("conditionCov", "noEffectCoord", "noEffectAxis", "tolX", "indefCovMat"),
    max.restarts = 100, sigma = 3, mu = floor(4L + floor(3 * log(n)) / 2), lambda = 4 * (4L + floor(3 * log(n)))
  )

  bestSoFarFit <- Inf
  bestSoFarAdditSolsFit <- rep(Inf, ADD_SOLS_SIZE)
  bestSoFarAdditSols <- matrix(NA, nrow = n, ncol = ADD_SOLS_SIZE)
  whenFoundAdditBestSols <- rep(0, ADD_SOLS_SIZE)



  # restart mechanism (IPOP-CMA-ES)
  restart.triggers <- getCMAESParameter(control, "restart.triggers", character(0L))

  stop.ons.names <- sapply(stop.ons, function(stop.on) stop.on$code)
  if (!isSubset(restart.triggers, stop.ons.names)) {
    stopf("Only codes / short names of active stopping conditions allowed as restart trigger, but '%s' are no stopping conditions.", collapse(setdiff(restart.triggers, stop.ons.names), sep = ", "))
  }
  restart.multiplier <- getCMAESParameter(control, "restart.multiplier", 2)
  assertNumber(restart.multiplier, lower = 1, na.ok = FALSE, finite = TRUE)
  max.restarts <- getCMAESParameter(control, "max.restarts", 0L)
  assertInt(max.restarts)

  # FIXME: default value should be derived from bounds
  sigma <- getCMAESParameter(control, "sigma", 0.5)
  assertNumber(sigma, lower = 0L, finite = TRUE)

  # Precompute E||N(0,I)||
  chi.n <- sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n^2))


  # bookkeep best individual
  best.param <- rep(NA, n)
  best.fitness <- Inf

  # set initial distribution mean
  m <- start.point

  # logs
  log.population <- getCMAESParameter(control, "log.population", FALSE)
  assertFlag(log.population, na.ok = FALSE)
  population.trace <- list()

  # init some termination criteria stuff
  iter <- 0L
  n.evals <- 0L
  start.time <- Sys.time()


  # somehow dirty trick to "really quit" if stopping condition is met and
  # now more restart should be triggered.
  do.terminate <- FALSE

  if (DATA_VER == "zestaw1RestFunNOneTimesCalcMOD") {
    lambda <- getCMAESParameter(control, "lambda", 4L + floor(3 * log(n)))
    MID_REST_MOD <- 10 + 30 * n / lambda
  }

  for (run in 0:max.restarts) {
    # population and offspring size
    if (run == 0) {
      lambda <- getCMAESParameter(control, "lambda", 4L + floor(3 * log(n)))
      assertInt(lambda, lower = 4)
      mu <- getCMAESParameter(control, "mu", floor(lambda / 2))
      assertInt(mu)
    } else {
      lambda <- getCMAESParameter(control, "lambda", 4L + floor(3 * log(n)))
      lambda <- ceiling(restart.multiplier^run * lambda)
      if (DATA_VER == "incrSigmaLikePop") {
        sigma <- getCMAESParameter(control, "sigma", 0.5)
        sigma <- ceiling(restart.multiplier^run * sigma)
      }

      if (DATA_VER != "withoutMuIncr") {
        # withoutMuIncr:
        mu <- floor(lambda / 2)
      }

      # startRestBest10
      if (DATA_VER == "startRestBest10" || DATA_VER == "startRestBest100") {
        n.evals <- n.evals + NUMBER_OF_TRIALS_FOR_START_SOL
        m <- selectStartSol()
      } else {
        if ((run == 1 || run == 2) && DATA_VER == "rest1i2zaburzMShiftBy5") { # przy pierwszym restarcie
          tmp <- rnorm(n, m, sigma)
          katy <- runif(n, 0, 2 * pi)
          tmp <- tmp + 5 * cos(katy)
          m <- ifelse(tmp < lb, lb, ifelse(tmp > ub, ub, tmp))
        } else {
          m <- runif(n, lb, ub) # default
        }
      }
    }

    # path for covariance matrix C and stepsize sigma
    pc <- rep(0, n)
    ps <- rep(0, n)

    # initialize recombination weights

    if (DATA_VER == "zestaw1RestFunNOneTimesCalcMOD" || DATA_VER == "zestaw1NoMidRestFunN" || DATA_VER == "zestaw1RestFunN" || DATA_VER == "zestaw1NoMiddle" || DATA_VER == "weights1" || DATA_VER == "zestaw1") {
      weights <- getCMAESParameter(control, "weights", log(mu + 1) - log(1:mu))
    } else {
      weights <- getCMAESParameter(control, "weights", log(mu + 0.5) - log(1:mu))
    }

    #  weights     <- controlParam("weights", log(mu+1) - log(1:mu)) #TU różnica!!!!!!!!!!

    if (any(weights < 0)) {
      stopf("All weights need to be positive, but there are %i negative ones.", sum(which(weights < 0)))
    }
    if (length(weights) != mu) {
      stopf("You need to pass %i 'weights', but you passed %i.", mu, length(weights))
    }

    # normalize weights
    weights <- weights / sum(weights)

    # variance-effectiveness / variance effective selection mass of sum w_i x_i
    mu.eff <- sum(weights)^2 / sum(weights^2) # chosen such that mu.eff ~ lambda/4


    # step-size control
    cs <- (mu.eff + 2) / (n + mu.eff + 5)


    ds <- 1 + 2 * max(0, sqrt((mu.eff - 1) / (n + 1)) - 1) + cs # damping factor

    # covariance matrix Adaptation parameters
    cc <- (4 + mu.eff / n) / (n + 4 + 2 * mu.eff / n)
    c1 <- 2 / ((n + 1.3)^2 + mu.eff)
    alpha.mu <- 2L
    cmu <- min(1 - c1, alpha.mu * (mu.eff - 2 + 1 / mu.eff) / ((n + 2)^2 + mu.eff))

    # covariance matrix
    if (DATA_VER != "incrSigmaLikePop") {
      sigma <- getCMAESParameter(control, "sigma", 0.5)
    }
    B <- diag(n)
    D <- diag(n)
    BD <- B %*% D
    C <- BD %*% t(BD) # C = B D^2 B^T = B B^T, since D equals I_n
    Cinvsqrt <- B %*% diag(1 / sqrt(diag(D))) %*% t(B)

    # no restart trigger fired until now
    restarting <- FALSE
    naprSigma <- 1
    brakZmFit <- 1
    oldMFit <- Inf
    mFit <- Inf
    # break inner loop if terminating stopping condition active or
    # restart triggered
    while (!restarting) {
      iter <- iter + 1L
      if (n.evals %% mod2putProgr == 0) {
        cat("#")
      }

      # create new population of search points
      arz <- matrix(rnorm(n * lambda), ncol = lambda) # ~ N(0, I)
      ary <- BD %*% arz # ~ N(0, C)
      # oldArx=arx
      arx <- m + sigma * ary # ~ N(m, sigma^2 C)

      arx.repaired <- ifelse(arx < lb, lb, ifelse(arx > ub, ub, arx))

      # Prepare penalization based on distance to repaired points (see Eq. 51)
      penalty.alpha <- 1L
      penalties <- penalty.alpha * colSums((arx - arx.repaired)^2)
      penalties[is.infinite(penalties)] <- .Machine$double.max / 2
      arx <- arx.repaired

      fitn <- objective.fun(t(arx)) + penalties

      # update evaluation
      n.evals <- n.evals + lambda

      # order fitness values
      fitn.ordered.idx <- order(fitn, decreasing = FALSE)
      fitn.ordered <- fitn[fitn.ordered.idx]


      # update mean value / center of mass
      new.pop.idx <- fitn.ordered.idx[1:mu]
      x.best <- arx[, new.pop.idx, drop = FALSE]
      m.old <- m
      m <- drop(x.best %*% weights)

      y.best <- ary[, new.pop.idx, drop = FALSE]
      y.w <- drop(y.best %*% weights)
      z.best <- arz[, new.pop.idx, drop = FALSE]
      z.w <- drop(z.best %*% weights)

      if (DATA_VER == "zestaw1RestFunNOneTimesCalcMOD") {
        if (iter %% MID_REST_MOD == 0) {
          oldMfit <- mFit
          mFit <- objective.fun(m)
          n.evals <- n.evals + 1
          if (abs(oldMfit - mFit) < 1e-8) {
            print(paste("restart brakZmianFit", n.evals))
            restarting <- TRUE
          }
        }
      }

      if (DATA_VER == "zestaw1RestFunN") {
        MID_REST_MOD <- 10 + 30 * n / lambda
        if (iter %% MID_REST_MOD == 0) {
          oldMfit <- mFit
          mFit <- objective.fun(m)
          n.evals <- n.evals + 1
          if (abs(oldMfit - mFit) < 1e-8) {
            print(paste("restart brakZmianFit", n.evals))
            restarting <- TRUE
          }
        }
      }
      if (DATA_VER == "minEvalues15_noFitChanAft15" || DATA_VER == "zestaw1") {
        # nowosc:minEvalues15_noFitChanAft10
        if (iter %% 15 == 0) {
          oldMfit <- mFit
          mFit <- objective.fun(m)
          n.evals <- n.evals + 1
          if (abs(oldMfit - mFit) < 1e-8) {
            print(paste("restart brakZmianFit", n.evals))
            restarting <- TRUE
          }
        }
      }

      # log population
      if (log.population) {
        population.trace[[iter]] <- x.best
      }
      if (length(stateDumpAtEv) > 0 && n.evals >= stateDumpAtEv[cecDumpIndx]) {
        # z dokladnoscia co do ewaluacji!
        if (stateDumpAtEv[cecDumpIndx] < n.evals) {
          ileEvZaDuzo <- n.evals - stateDumpAtEv[cecDumpIndx]
          if (DATA_VER == "startRestBest10" || DATA_VER == "startRestBest100") {
            if (ileEvZaDuzo >= lambda) {
              ileEvZaDuzo <- lambda - 1
            }
          }

          tmpFitn <- fitn[1:(lambda - ileEvZaDuzo)]
          minIndx <- which.min(tmpFitn)
          best.fitness <- tmpFitn[minIndx]
          best.param <- arx[, minIndx]
        } else {
          best.fitness <- fitn.ordered[1]
          best.param <- arx[, fitn.ordered.idx[1]]
        }
        if (best.fitness < bestSoFarFit) {
          bestSoFarFit <- best.fitness
          bestSoFarSol <- best.param
          whenFoundBestSol <- n.evals
        }

        bestSoFarSolHist[, cecDumpIndx] <- bestSoFarSol
        bestSoFarFitHist[cecDumpIndx] <- bestSoFarFit
        whenFoundBestSolHist[cecDumpIndx] <- whenFoundBestSol

        if (!is.null(genAddSolsFun)) {
          addSols <- matrix(0, n, 4)
          n.evals <- n.evals + 1 # ToDo!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          addSols[, 4] <- m
          addSolsFit <- objFun(t(addSols))
          betterAdditSolsFitIndxs <- addSolsFit < bestSoFarAdditSolsFit
          bestSoFarAdditSolsFit[betterAdditSolsFitIndxs] <- addSolsFit[betterAdditSolsFitIndxs]
          bestSoFarAdditSols[, betterAdditSolsFitIndxs] <- addSols[, betterAdditSolsFitIndxs]
          whenFoundAdditBestSols[betterAdditSolsFitIndxs] <- n.evals

          bestSoFarAdditSolsHist[, , cecDumpIndx] <- bestSoFarAdditSols
          bestSoFarAdditSolsFitHist[, cecDumpIndx] <- bestSoFarAdditSolsFit
          whenFoundAdditBestSolsHist[, cecDumpIndx] <- whenFoundAdditBestSols
        }

        cecDumpIndx <- cecDumpIndx + 1
      } else { # if state Dump
        best.fitness <- fitn.ordered[1]
        best.param <- arx[, fitn.ordered.idx[1]]
        if (best.fitness < bestSoFarFit) {
          bestSoFarFit <- best.fitness
          bestSoFarSol <- best.param
          whenFoundBestSol <- n.evals
        }
      }







      # Update evolution path with cumulative step-size Adaptation (CSA) / path length control
      # For an explanation of the last factor see appendix A in https://www.lri.fr/~hansen/cmatutorial.pdf
      ps <- (1 - cs) * ps + sqrt(cs * (2 - cs) * mu.eff) * (Cinvsqrt %*% y.w)
      h.sigma <- as.integer(norm2(ps) / sqrt(1 - (1 - cs)^(2 * (iter + 1))) < chi.n * (1.4 + 2 / (n + 1)))

      # Update covariance matrix
      pc <- (1 - cc) * pc + h.sigma * sqrt(cc * (2 - cc) * mu.eff) * y.w
      y <- BD %*% z.best
      delta.h.sigma <- as.numeric((1 - h.sigma) * cc * (2 - cc) <= 1)
      C <- (1 - c1 - cmu) * C + c1 * (pc %*% t(pc) + delta.h.sigma * C) + cmu * y %*% diag(weights) %*% t(y)

      # Update step-size sigma
      sigma <- sigma * exp(cs / ds * ((norm2(ps) / chi.n) - 1))

      if (DATA_VER == "sig120really") {
        maxdx <- 130
      } else {
        if (DATA_VER == "zestaw1RestFunNOneTimesCalcMOD" || DATA_VER == "zestaw1NoMidRestFunN" || DATA_VER == "zestaw1RestFunN" || DATA_VER == "zestaw1NoMiddle" || DATA_VER == "maxdxNa4" || DATA_VER == "zestaw1") {
          maxdx <- 160 / 4
        } else {
          maxdx <- 160 / 2 # default
        }
      }

      if (any(sigma * sqrt(diag(C)) > maxdx)) {
        print(paste(naprSigma, "Naprawa sigma", sigma, min(maxdx / sqrt(diag(C)))))
        if (DATA_VER == "zestaw1RestFunNOneTimesCalcMOD" || DATA_VER == "zestaw1NoMidRestFunN" || DATA_VER == "zestaw1RestFunN" || DATA_VER == "zestaw1NoMiddle" || DATA_VER == "naprSigmaDiv4" || DATA_VER == "zestaw1") {
          SIGMA_DIVIDER <- 4
        } else {
          SIGMA_DIVIDER <- 2
        }
        sigma <- min(maxdx / sqrt(diag(C))) / SIGMA_DIVIDER
        naprSigma <- naprSigma + 1
        if (DATA_VER == "naprSigma2") {
          PROG <- 2
        } else {
          PROG <- 5
        }
        if (naprSigma > PROG) {
          print(paste("Restart przez napr sigma", n.evals))
          restarting <- TRUE
        }
      }

      # Finally do decomposition C = B D^2 B^T
      e <- eigen(C, symmetric = TRUE)


      # escape flat fitness values
      if (fitn.ordered[1L] == fitn.ordered[ceiling(0.7 * lambda)]) {
        print(paste("Flat fittness:", sigma, sigma * exp(0.2 + cs / ds), n.evals))

        if (DATA_VER == "zestaw1RestFunNOneTimesCalcMOD" || DATA_VER == "zestaw1NoMidRestFunN" || DATA_VER == "zestaw1RestFunN" || DATA_VER == "zestaw1NoMiddle" || DATA_VER == "flatFit2xSig" || DATA_VER == "zestaw1") {
          sigma <- sigma * 2 * exp(0.2 + cs / ds)
        } else {
          sigma <- sigma * exp(0.2 + cs / ds)
        }
      }

      # CHECK STOPPING CONDITIONS
      # =========================
      stop.obj <- checkStoppingConditions(stop.ons)
      B <- e$vectors
      D <- diag(sqrt(e$values))
      BD <- B %*% D
      Cinvsqrt <- B %*% diag(1 / diag(D)) %*% t(B) # update C^-1/2



      n.stop.codes <- length(stop.obj$codes)
      if (n.stop.codes > 0 && any(stop.obj$codes %in% list("maxEvals", "optValTol"))) {
        print(paste("Koniec powod:", stop.obj$codes))
        do.terminate <- TRUE
        break
      }
      if (max.restarts > 0L && any(stop.obj$codes %in% restart.triggers)) {
        print(paste("Restart:", n.evals, " powod:", stop.obj$codes))
        n.stop.codes <- sum(!(stop.obj$codes %in% restart.triggers))
        restarting <- TRUE
      }

      # check if CMA-ES should really quit, i.e., is there a stopping condition,
      # that is active and does not trigger a restart?
      if (!restarting && (n.stop.codes > 0L)) {
        print(paste("Koniec powod:", stop.obj$codes))
        do.terminate <- TRUE
        break
      }
    }

    # really quit without more restarts
    if (do.terminate) {
      break
    }
  }

  # callMonitor(monitor, "after")
  print(paste("koniec:", n.evals))



  cat("\n")
  toRet <- list(bestSoFarFit = bestSoFarFit, bestSoFarSol = bestSoFarSol, whenFoundBestSol = whenFoundBestSol)
  if (length(stateDumpAtEv) > 0) {
    toRet$bestSoFarSolHist <- bestSoFarSolHist
    toRet$bestSoFarFitHist <- bestSoFarFitHist
    toRet$whenFoundBestSolHist <- whenFoundBestSolHist

    if (!is.null(genAddSolsFun)) {
      toRet$bestSoFarAdditSolsHist <- bestSoFarAdditSolsHist
      toRet$bestSoFarAdditSolsFitHist <- bestSoFarAdditSolsFitHist
      toRet$whenFoundAdditBestSolsHist <- whenFoundAdditBestSolsHist
    }
  }
  return(toRet)
}

############################
############################
############################
############################

getCMAESParameter <- function(control, what, default) {
  return(coalesce(control[[what]], default))
}

norm2 <- function(x) {
  return(drop(sqrt(crossprod(x))))
}

# not used
stopOnNoChangeInSrFit <- function() {
  return(makeStoppingCondition(
    name = "noChangeInSrFit",
    message = "No chagne in fitness.",
    stop.fun = function(envir = parent.frame()) {
      e.values <- envir$e$values
      return(any(is.na(e.values)) || min(e.values) <= 1e-9 || max(e.values) > 1e14 * min(e.values))
    }
  ))
}

stopOnIndefCovMat <- function() {
  return(makeStoppingCondition(
    name = "indefCovMat",
    message = "Covariance matrix is not numerically positive definite.",
    stop.fun = function(envir = parent.frame()) {
      minEvalues <- 1e-15
      if (any(is.na(envir$e$values)) || min(envir$e$values) <= minEvalues || max(envir$e$values) > 1e14 * min(envir$e$values)) {
        if (any(envir$e$values <= minEvalues)) {
          envir$e$values[envir$e$values <= minEvalues] <- minEvalues #* 10 zestaw1mult10

          print("korekcja e$values na 1e-8")
        }
      }
      return(any(is.na(envir$e$values)) || max(envir$e$values) > 1e14 * min(envir$e$values))
    }
  ))
}

stopOnOptValue <- function(opt.value, tol = 1e-8) {
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

stopOnMaxEvals <- function(max.evals) {
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

# FIXME: default value is 1e-12 * sigma. Here we have no access to the sigma value.
stopOnTolX <- function(tol = 1e-9) {
  assertInt(tol, na.ok = FALSE)
  force(tol)
  return(makeStoppingCondition(
    name = "tolX",
    message = sprintf("Standard deviation below tolerance in all coordinates."),
    stop.fun = function(envir = parent.frame()) {
      return(all(diag(sqrt(envir$e$values)) < tol) && all((envir$sigma * envir$p.c) < tol))
    }
  ))
}


stopOnNoEffectAxis <- function() {
  return(makeStoppingCondition(
    name = "noEffectAxis",
    message = "Addition of 0.1 times sigma does not change mean value.",
    stop.fun = function(envir = parent.frame()) {
      ii <- (envir$iter %% envir$n) + 1L
      ui <- envir$e$vectors[, ii]
      lambdai <- sqrt(envir$e$values[ii])
      m <- envir$m
      return(sum((m - (m + 0.3 * envir$sigma * lambdai * ui))^2) < .Machine$double.eps)
    }
  ))
}

stopOnNoEffectCoord <- function() {
  return(makeStoppingCondition(
    name = "noEffectCoord",
    message = "Addition of 0.2 times sigma in any coordinate does not change mean value.",
    stop.fun = function(envir = parent.frame()) {
      m <- envir$m
      return(sum((m - (m + 0.2 * envir$sigma))^2) < .Machine$double.eps)
    }
  ))
}

stopOnCondCov <- function(tol = 1e14) {
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

makeStoppingCondition <- function(name, message, stop.fun, code = name, control = list()) {
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

shouldStop <- function(x, envir) {
  UseMethod("shouldStop")
}

shouldStop.cma_stopping_condition <- function(x, envir) {
  return(x$stop.fun(envir))
}

checkStoppingConditions <- function(stop.ons, envir = parent.frame()) {
  assertList(stop.ons, min.len = 1L, types = "cma_stopping_condition")
  stop.msgs <- character(0L)
  codes <- character(0L)
  for (stop.on in stop.ons) {
    if (shouldStop(stop.on, envir = envir)) {
      stop.msgs <- c(stop.msgs, stop.on$message)
      codes <- c(codes, stop.on$code)
      # since some stopping conditions need a "correct" covariance matrix
      # we stop here if the first condition is met (infefCovMat is first)
      break
    }
  }
  return(list(stop.msgs = stop.msgs, codes = codes))
}

#' @title Return list of default stopping conditions.
#'
#' @description Default stopping conditions which are active in the reference
#' implementation by Nico Hansen in Python.
#'
#' @return [\code{list}]
#' @export
getDefaultStoppingConditions <- function() {
  return(
    list(
      stopOnTolX(),
      stopOnNoEffectAxis(),
      stopOnNoEffectCoord(),
      stopOnCondCov()
    )
  )
}