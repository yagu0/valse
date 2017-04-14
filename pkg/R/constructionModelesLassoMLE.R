#' constructionModelesLassoMLE 
#'
#' Construct a collection of models with the Lasso-MLE procedure.
#' 
#' @param phiInit an initialization for phi, get by initSmallEM.R
#' @param rhoInit an initialization for rho, get by initSmallEM.R
#' @param piInit an initialization for pi, get by initSmallEM.R
#' @param gamInit an initialization for gam, get by initSmallEM.R
#' @param mini integer, minimum number of iterations in the EM algorithm, by default = 10
#' @param maxi integer, maximum number of iterations in the EM algorithm, by default = 100
#' @param gamma integer for the power in the penaly, by default = 1
#' @param X matrix of covariates (of size n*p)
#' @param Y matrix of responses (of size n*m)
#' @param eps real, threshold to say the EM algorithm converges, by default = 1e-4
#' @param S output of selectVariables.R
#' @param ncores Number of cores, by default = 3
#' @param fast TRUE to use compiled C code, FALSE for R code only
#' @param verbose TRUE to show some execution traces
#' 
#' @return a list with several models, defined by phi, rho, pi, llh
#'
#' @export
constructionModelesLassoMLE <- function(phiInit, rhoInit, piInit, gamInit, mini, 
  maxi, gamma, X, Y, eps, S, ncores = 3, fast = TRUE, verbose = FALSE)
  {
  if (ncores > 1)
  {
    cl <- parallel::makeCluster(ncores, outfile = "")
    parallel::clusterExport(cl, envir = environment(), varlist = c("phiInit", 
      "rhoInit", "gamInit", "mini", "maxi", "gamma", "X", "Y", "eps", "S", 
      "ncores", "fast", "verbose"))
  }
  
  # Individual model computation
  computeAtLambda <- function(lambda)
  {
    if (ncores > 1) 
      require("valse")  #nodes start with an empty environment
    
    if (verbose) 
      print(paste("Computations for lambda=", lambda))
    
    n <- dim(X)[1]
    p <- dim(phiInit)[1]
    m <- dim(phiInit)[2]
    k <- dim(phiInit)[3]
    sel.lambda <- S[[lambda]]$selected
    # col.sel = which(colSums(sel.lambda)!=0) #if boolean matrix
    col.sel <- which(sapply(sel.lambda, length) > 0)  #if list of selected vars
    if (length(col.sel) == 0) 
      return(NULL)
    
    # lambda == 0 because we compute the EMV: no penalization here
    res <- EMGLLF(phiInit[col.sel, , ], rhoInit, piInit, gamInit, mini, maxi, 
      gamma, 0, X[, col.sel], Y, eps, fast)
    
    # Eval dimension from the result + selected
    phiLambda2 <- res$phi
    rhoLambda <- res$rho
    piLambda <- res$pi
    phiLambda <- array(0, dim = c(p, m, k))
    for (j in seq_along(col.sel)) phiLambda[col.sel[j], sel.lambda[[j]], ] <- phiLambda2[j, 
      sel.lambda[[j]], ]
    dimension <- length(unlist(sel.lambda))
    
    # Computation of the loglikelihood
    densite <- vector("double", n)
    for (r in 1:k)
    {
      if (length(col.sel) == 1)
      {
        delta <- (Y %*% rhoLambda[, , r] - (X[, col.sel] %*% t(phiLambda[col.sel, 
          , r])))
      } else delta <- (Y %*% rhoLambda[, , r] - (X[, col.sel] %*% phiLambda[col.sel, 
        , r]))
      densite <- densite + piLambda[r] * det(rhoLambda[, , r])/(sqrt(2 * base::pi))^m * 
        exp(-diag(tcrossprod(delta))/2)
    }
    llhLambda <- c(sum(log(densite)), (dimension + m + 1) * k - 1)
    list(phi = phiLambda, rho = rhoLambda, pi = piLambda, llh = llhLambda)
  }
  
  # For each lambda, computation of the parameters
  out <- if (ncores > 1) 
    parLapply(cl, 1:length(S), computeAtLambda) else lapply(1:length(S), computeAtLambda)
  
  if (ncores > 1) 
    parallel::stopCluster(cl)
  
  out
}
