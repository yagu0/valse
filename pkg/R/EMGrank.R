#' EMGrank
#'
#' Description de EMGrank
#'
#' @param Pi Parametre de proportion
#' @param Rho Parametre initial de variance renormalisé
#' @param mini Nombre minimal d'itérations dans l'algorithme EM
#' @param maxi Nombre maximal d'itérations dans l'algorithme EM
#' @param X Régresseurs
#' @param Y Réponse
#' @param eps Seuil pour accepter la convergence
#' @param rank Vecteur des rangs possibles
#'
#' @return A list ...
#'   phi : parametre de moyenne renormalisé, calculé par l'EM
#'   LLF : log vraisemblance associé à cet échantillon, pour les valeurs estimées des paramètres
#'
#' @export
EMGrank <- function(Pi, Rho, mini, maxi, X, Y, eps, rank, fast = TRUE)
{
  if (!fast)
  {
    # Function in R
    return(.EMGrank_R(Pi, Rho, mini, maxi, X, Y, eps, rank))
  }

  # Function in C
  n <- nrow(X)  #nombre d'echantillons
  p <- ncol(X)  #nombre de covariables
  m <- ncol(Y)  #taille de Y (multivarié)
  k <- length(Pi)  #nombre de composantes dans le mélange
  .Call("EMGrank", Pi, Rho, mini, maxi, X, Y, eps, as.integer(rank), phi = double(p * m * k), 
    LLF = double(1), n, p, m, k, PACKAGE = "valse")
}

# helper to always have matrices as arg (TODO: put this elsewhere? improve?)  -->
# Yes, we should use by-columns storage everywhere... [later!]
matricize <- function(X)
{
  if (!is.matrix(X)) 
    return(t(as.matrix(X)))
  return(X)
}

# R version - slow but easy to read
.EMGrank_R <- function(Pi, Rho, mini, maxi, X, Y, eps, rank)
{
  # matrix dimensions
  n <- nrow(X)
  p <- ncol(X)
  m <- ncol(Y)
  k <- length(Pi)

  # init outputs
  phi <- array(0, dim = c(p, m, k))
  Z <- rep(1, n)
  LLF <- 0

  # local variables
  Phi <- array(0, dim = c(p, m, k))
  deltaPhi <- c()
  sumDeltaPhi <- 0
  deltaPhiBufferSize <- 20
  
  # main loop
  ite <- 1
  while (ite <= mini || (ite <= maxi && sumDeltaPhi > eps))
  {
    # M step: update for Beta ( and then phi)
    for (r in 1:k)
    {
      Z_indice <- seq_len(n)[Z == r] #indices where Z == r
      if (length(Z_indice) == 0) 
        next
      # U,S,V = SVD of (t(Xr)Xr)^{-1} * t(Xr) * Yr
      s <- svd(MASS::ginv(crossprod(matricize(X[Z_indice, ]))) %*% 
                 crossprod(matricize(X[Z_indice, ]), matricize(Y[Z_indice, ])))
      S <- s$d
      # Set m-rank(r) singular values to zero, and recompose best rank(r) approximation
      # of the initial product
      if (rank[r] < length(S)) 
        S[(rank[r] + 1):length(S)] <- 0
      phi[, , r] <- s$u %*% diag(S) %*% t(s$v) %*% Rho[, , r]
    }

    # Step E and computation of the loglikelihood
    sumLogLLF2 <- 0
    for (i in seq_len(n))
    {
      sumLLF1 <- 0
      maxLogGamIR <- -Inf
      for (r in seq_len(k))
      {
        dotProduct <- tcrossprod(Y[i, ] %*% Rho[, , r] - X[i, ] %*% phi[, , r])
        logGamIR <- log(Pi[r]) + log(gdet(Rho[, , r])) - 0.5 * dotProduct
        # Z[i] = index of max (gam[i,])
        if (logGamIR > maxLogGamIR)
        {
          Z[i] <- r
          maxLogGamIR <- logGamIR
        }
        sumLLF1 <- sumLLF1 + exp(logGamIR)/(2 * pi)^(m/2)
      }
      sumLogLLF2 <- sumLogLLF2 + log(sumLLF1)
    }

    LLF <- -1/n * sumLogLLF2

    # update distance parameter to check algorithm convergence (delta(phi, Phi))
    deltaPhi <- c(deltaPhi, max((abs(phi - Phi))/(1 + abs(phi)))) #TODO: explain?
    if (length(deltaPhi) > deltaPhiBufferSize) 
      deltaPhi <- deltaPhi[2:length(deltaPhi)]
    sumDeltaPhi <- sum(abs(deltaPhi))

    # update other local variables
    Phi <- phi
    ite <- ite + 1
  }
  return(list(phi = phi, LLF = LLF))
}
