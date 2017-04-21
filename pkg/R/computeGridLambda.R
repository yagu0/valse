#' computeGridLambda 
#'
#' Construct the data-driven grid for the regularization parameters used for the Lasso estimator
#'
#' @param phiInit value for phi
#' @param rhoInit\tvalue for rho
#' @param piInit\tvalue for pi
#' @param gamInit value for gamma
#' @param X matrix of covariates (of size n*p)
#' @param Y matrix of responses (of size n*m)
#' @param gamma power of weights in the penalty
#' @param mini minimum number of iterations in EM algorithm
#' @param maxi maximum number of iterations in EM algorithm
#' @param tau threshold to stop EM algorithm
#'
#' @return the grid of regularization parameters
#'
#' @export
computeGridLambda <- function(phiInit, rhoInit, piInit, gamInit, X, Y, gamma, mini, 
  maxi, tau, fast)
{
  n <- nrow(X)
  p <- ncol(X)
  m <- ncol(Y)
  k <- length(piInit)

  list_EMG <- EMGLLF(phiInit, rhoInit, piInit, gamInit, mini, maxi, gamma, lambda = 0, 
    X, Y, tau, fast)
  grid <- array(0, dim = c(p, m, k))
  for (i in 1:p)
  {
    for (j in 1:m)
      grid[i, j, ] <- abs(list_EMG$S[i, j, ])/(n * list_EMG$pi^gamma)
  }
  sort(unique(grid))
}
