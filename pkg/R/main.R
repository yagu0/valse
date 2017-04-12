#' valse
#'
#' Main function
#'
#' @param X matrix of covariates (of size n*p)
#' @param Y matrix of responses (of size n*m)
#' @param procedure among 'LassoMLE' or 'LassoRank'
#' @param selecMod method to select a model among 'DDSE', 'DJump', 'BIC' or 'AIC'
#' @param gamma integer for the power in the penaly, by default = 1
#' @param mini integer, minimum number of iterations in the EM algorithm, by default = 10
#' @param maxi integer, maximum number of iterations in the EM algorithm, by default = 100
#' @param eps real, threshold to say the EM algorithm converges, by default = 1e-4
#' @param kmin integer, minimum number of clusters, by default = 2
#' @param kmax integer, maximum number of clusters, by default = 10
#' @param rang.min integer, minimum rank in the low rank procedure, by default = 1
#' @param rang.max integer, maximum rank in the
#' @param ncores_outer Number of cores for the outer loop on k
#' @param ncores_inner Number of cores for the inner loop on lambda
#' @param size_coll_mod (Maximum) size of a collection of models
#' @param fast TRUE to use compiled C code, FALSE for R code only
#' @param verbose TRUE to show some execution traces
#'
#' @return a list with estimators of parameters
#'
#' @examples
#' #TODO: a few examples
#' @export
valse = function(X, Y, procedure='LassoMLE', selecMod='DDSE', gamma=1, mini=10, maxi=50,
                 eps=1e-4, kmin=2, kmax=2, rank.min=1, rank.max=10, ncores_outer=1, ncores_inner=1,
                 size_coll_mod=10, fast=TRUE, verbose=FALSE, plot = TRUE)
{
  p = dim(X)[2]
  m = dim(Y)[2]
  n = dim(X)[1]
  
  if (verbose)
    print("main loop: over all k and all lambda")
  
  if (ncores_outer > 1)
  {
    cl = parallel::makeCluster(ncores_outer, outfile='')
    parallel::clusterExport( cl=cl, envir=environment(), varlist=c("X","Y","procedure",
                                                                   "selecMod","gamma","mini","maxi","eps","kmin","kmax","rang.min","rang.max",
                                                                   "ncores_outer","ncores_inner","verbose","p","m") )
  }
  
  # Compute models with k components
  computeModels <- function(k)
  {
    if (ncores_outer > 1)
      require("valse") #nodes start with an empty environment
    
    if (verbose)
      print(paste("Parameters initialization for k =",k))
    #smallEM initializes parameters by k-means and regression model in each component,
    #doing this 20 times, and keeping the values maximizing the likelihood after 10
    #iterations of the EM algorithm.
    P = initSmallEM(k, X, Y)
    grid_lambda <- computeGridLambda(P$phiInit, P$rhoInit, P$piInit, P$gamInit, X, Y,
                                     gamma, mini, maxi, eps, fast)
    if (length(grid_lambda)>size_coll_mod)
      grid_lambda = grid_lambda[seq(1, length(grid_lambda), length.out = size_coll_mod)]
    
    if (verbose)
      print("Compute relevant parameters")
    #select variables according to each regularization parameter
    #from the grid: S$selected corresponding to selected variables
    S = selectVariables(P$phiInit, P$rhoInit, P$piInit, P$gamInit, mini, maxi, gamma,
                        grid_lambda, X, Y, 1e-8, eps, ncores_inner, fast) #TODO: 1e-8 as arg?! eps?
    
    if (procedure == 'LassoMLE')
    {
      if (verbose)
        print('run the procedure Lasso-MLE')
      #compute parameter estimations, with the Maximum Likelihood
      #Estimator, restricted on selected variables.
      models <- constructionModelesLassoMLE(P$phiInit, P$rhoInit, P$piInit, P$gamInit,
                                            mini, maxi, gamma, X, Y, thresh, eps, S, ncores_inner, fast, verbose)
    }
    else
    {
      if (verbose)
        print('run the procedure Lasso-Rank')
      #compute parameter estimations, with the Low Rank
      #Estimator, restricted on selected variables.
      models <- constructionModelesLassoRank(S$Pi, S$Rho, mini, maxi, X, Y, eps, S,
                                             rank.min, rank.max, ncores_inner, fast, verbose)
    }
    #warning! Some models are NULL after running selectVariables
    models = models[sapply(models, function(cell) !is.null(cell))]
    models
  }
  
  # List (index k) of lists (index lambda) of models
  models_list <-
    if (ncores_outer > 1)
      parLapply(cl, kmin:kmax, computeModels)
  else
    lapply(kmin:kmax, computeModels)
  if (ncores_outer > 1)
    parallel::stopCluster(cl)
  
  if (! requireNamespace("capushe", quietly=TRUE))
  {
    warning("'capushe' not available: returning all models")
    return (models_list)
  }
  
  # Get summary "tableauRecap" from models
  tableauRecap = do.call( rbind, lapply( seq_along(models_list), function(i) {
    models <- models_list[[i]]
    #For a collection of models (same k, several lambda):
    LLH <- sapply( models, function(model) model$llh[1] )
    k = length(models[[1]]$pi)
    sumPen = sapply(models, function(model)
      k*(dim(model$rho)[1]+sum(model$phi[,,1]!=0)+1)-1)
    data.frame(model=paste(i,".",seq_along(models),sep=""),
               pen=sumPen/n, complexity=sumPen, contrast=-LLH)
  } ) )
  
  print(tableauRecap)
  tableauRecap = tableauRecap[which(tableauRecap[,4]!= Inf),]
  modSel = capushe::capushe(tableauRecap, n)
  indModSel <-
    if (selecMod == 'DDSE')
      as.numeric(modSel@DDSE@model)
  else if (selecMod == 'Djump')
    as.numeric(modSel@Djump@model)
  else if (selecMod == 'BIC')
    modSel@BIC_capushe$model
  else if (selecMod == 'AIC')
    modSel@AIC_capushe$model
  
  mod = as.character(tableauRecap[indModSel,1])
  listMod = as.integer(unlist(strsplit(mod, "[.]")))
  modelSel = models_list[[listMod[1]]][[listMod[2]]]
  
  ##Affectations
  Gam = matrix(0, ncol = length(modelSel$pi), nrow = n)
  for (i in 1:n){
    for (r in 1:length(modelSel$pi)){
      sqNorm2 = sum( (Y[i,]%*%modelSel$rho[,,r]-X[i,]%*%modelSel$phi[,,r])^2 )
      Gam[i,r] = modelSel$pi[r] * exp(-0.5*sqNorm2)* det(modelSel$rho[,,r])
    }
  }
  Gam = Gam/rowSums(Gam)
  modelSel$affec = apply(Gam, 1,which.max)
  modelSel$proba = Gam
  
  if (plot){
    print(plot_valse(modelSel,n))
  }
  
  return(modelSel)
}
