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
#'
#' @return a list with estimators of parameters
#'
#' @examples
#' #TODO: a few examples
#' @export
valse = function(X, Y, procedure='LassoMLE', selecMod='DDSE', gamma=1, mini=10, maxi=50,
	eps=1e-4, kmin=2, kmax=2, rang.min=1, rang.max=10, ncores_outer=1, ncores_inner=3,
	verbose=FALSE)
{
  p = dim(X)[2]
  m = dim(Y)[2]
  n = dim(X)[1]

  if (verbose)
		print("main loop: over all k and all lambda")

	if (ncores_outer > 1)
	{
		cl = parallel::makeCluster(ncores_outer)
		parallel::clusterExport( cl=cl, envir=environment(), varlist=c("X","Y","procedure",
			"selecMod","gamma","mini","maxi","eps","kmin","kmax","rang.min","rang.max",
			"ncores_outer","ncores_inner","verbose","p","m","k","tableauRecap") )
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
			gamma, mini, maxi, eps)
		# TODO: 100 = magic number
    if (length(grid_lambda)>100)
      grid_lambda = grid_lambda[seq(1, length(grid_lambda), length.out = 100)]

		if (verbose)
			print("Compute relevant parameters")
    #select variables according to each regularization parameter
    #from the grid: S$selected corresponding to selected variables
    S = selectVariables(P$phiInit, P$rhoInit, P$piInit, P$gamInit, mini, maxi, gamma,
			grid_lambda, X, Y, 1e-8, eps, ncores_inner) #TODO: 1e-8 as arg?! eps?

    if (procedure == 'LassoMLE')
		{
      if (verbose)
				print('run the procedure Lasso-MLE')
      #compute parameter estimations, with the Maximum Likelihood
      #Estimator, restricted on selected variables.
      models <- constructionModelesLassoMLE(phiInit, rhoInit, piInit, gamInit, mini,
				maxi, gamma, X, Y, thresh, eps, S$selected, ncores_inner, verbose)
    }
		else
		{
      if (verbose)
				print('run the procedure Lasso-Rank')
      #compute parameter estimations, with the Low Rank
      #Estimator, restricted on selected variables.
      models <- constructionModelesLassoRank(S$Pi, S$Rho, mini, maxi, X, Y, eps, A1,
				rank.min, rank.max, ncores_inner, verbose)
    }
    models
  }

	# List (index k) of lists (index lambda) of models
	models_list <-
		if (ncores_k > 1)
			parLapply(cl, kmin:kmax, computeModels)
		else
			lapply(kmin:kmax, computeModels)
	if (ncores_k > 1)
		parallel::stopCluster(cl)

	if (! requireNamespace("capushe", quietly=TRUE))
	{
		warning("'capushe' not available: returning all models")
		return (models_list)
	}

	# Get summary "tableauRecap" from models ; TODO: jusqu'à ligne 114 à mon avis là c'est faux :/
	tableauRecap = t( sapply( models_list, function(models) {
		llh = do.call(rbind, lapply(models, function(model) model$llh)
    LLH = llh[-1,1]
    D = llh[-1,2]
		c(LLH, D, rep(k, length(model)), 1:length(model))
	) } ) )
	if (verbose)
		print('Model selection')
  tableauRecap = tableauRecap[rowSums(tableauRecap[, 2:4])!=0,]
  tableauRecap = tableauRecap[!is.infinite(tableauRecap[,1]),]
  data = cbind(1:dim(tableauRecap)[1], tableauRecap[,2], tableauRecap[,2], tableauRecap[,1])

  modSel = capushe::capushe(data, n)
  indModSel <-
		if (selecMod == 'DDSE')
			as.numeric(modSel@DDSE@model)
		else if (selecMod == 'Djump')
			as.numeric(modSel@Djump@model)
		else if (selecMod == 'BIC')
			modSel@BIC_capushe$model
		else if (selecMod == 'AIC')
			modSel@AIC_capushe$model
  model[[tableauRecap[indModSel,3]]][[tableauRecap[indModSel,4]]]
}
