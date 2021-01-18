noBounds <- function(mutants, lower, upper) {
  return(mutants)
}

bounceBack <- function(newPop, lower, upper) {
  # print(sum(newPop>upper))
  for (i in 1:(ncol(newPop))) {
    indxs <- which(newPop[, i] < lower)
    newPop[indxs, i] <- lower[indxs] + abs(lower[indxs] - newPop[indxs, i]) %% (upper[indxs] - lower[indxs])
    indxs <- which(newPop[, i] > upper)
    newPop[indxs, i] <- upper[indxs] - abs(upper[indxs] - newPop[indxs, i]) %% (upper[indxs] - lower[indxs])
  }
  return(newPop)
}

bounceBackAndPen <- function(newPop, lower, upper) {
  k <- 1
  pen <- numeric(ncol(newPop))
  for (i in 1:(ncol(newPop))) {
    indxs <- which(newPop[, i] < lower)
    pen[i] <- k * sum(lower[indxs] - newPop[indxs, i])
    newPop[indxs, i] <- lower[indxs] + abs(lower[indxs] - newPop[indxs, i]) %% (upper[indxs] - lower[indxs])
    indxs <- which(newPop[, i] > upper)
    pen[i] <- pen[i] + k * sum(newPop[indxs, i] - upper[indxs])
    newPop[indxs, i] <- upper[indxs] - abs(upper[indxs] - newPop[indxs, i]) %% (upper[indxs] - lower[indxs])
  }
  return(list(pop = newPop, pen = pen))
}
