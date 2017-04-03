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
valse = function(X,Y,procedure = 'LassoMLE',selecMod = 'DDSE',gamma = 1,mini = 10,
                 maxi = 50,eps = 1e-4,kmin = 2,kmax = 2,
                 rang.min = 1,rang.max = 10, ncores_k=1, ncores_lambda=3, verbose=FALSE)
{
  p = dim(X)[2]
  m = dim(Y)[2]
  n = dim(X)[1]

  tableauRecap = list()
  if (verbose)
		print("main loop: over all k and all lambda")

	if (ncores_k > 1)
	{
		cl = parallel::makeCluster(ncores_k)
		parallel::clusterExport( cl=cl, envir=environment(), varlist=c("X","Y","procedure",
			"selecMod","gamma","mini","maxi","eps","kmin","kmax","rang.min","rang.max",
			"ncores_k","ncores_lambda","verbose","p","m","k","tableauRecap") )
	}

	# Compute model with k components
	computeModel <- function(k)
	{
		if (ncores_k > 1)
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
    #from the grid: A1 corresponding to selected variables, and
    #A2 corresponding to unselected variables.
    S = selectVariables(P$phiInit,P$rhoInit,P$piInit,P$gamInit,mini,maxi,gamma,
			grid_lambda,X,Y,1e-8,eps,ncores_lambda)

    if (procedure == 'LassoMLE')
		{
      if (verbose)
				print('run the procedure Lasso-MLE')
      #compute parameter estimations, with the Maximum Likelihood
      #Estimator, restricted on selected variables.
      model = constructionModelesLassoMLE(phiInit, rhoInit, piInit, gamInit, mini,
				maxi, gamma, X, Y, thresh, eps, S$selected)
      llh = matrix(ncol = 2)
      for (l in seq_along(model[[k]]))
        llh = rbind(llh, model[[k]][[l]]$llh)
      LLH = llh[-1,1]
      D = llh[-1,2]
    }
		else
		{
      if (verbose)
				print('run the procedure Lasso-Rank')
      #compute parameter estimations, with the Low Rank
      #Estimator, restricted on selected variables.
      model = constructionModelesLassoRank(S$Pi, S$Rho, mini, maxi, X, Y, eps, A1,
				rank.min, rank.max)

      ################################################
      ### Regarder la SUITE  
      phi = runProcedure2()$phi
      Phi2 = Phi
      if (dim(Phi2)[1] == 0)
        Phi[, , 1:k,] <- phi
      else
      {
        Phi <- array(0, dim = c(p, m, kmax, dim(Phi2)[4] + dim(phi)[4]))
        Phi[, , 1:(dim(Phi2)[3]), 1:(dim(Phi2)[4])] <<- Phi2
        Phi[, , 1:k,-(1:(dim(Phi2)[4]))] <<- phi
      }
    }
    tableauRecap[[k]] = matrix(c(LLH, D, rep(k, length(model[[k]])), 1:length(model[[k]])), ncol = 4))
  }

	model <-
		if (ncores_k > 1)
			parLapply(cl, kmin:kmax, computeModel)
		else
			lapply(kmin:kmax, computeModel)
	if (ncores_k > 1)
		parallel::stopCluster(cl)

	if (verbose)
		print('Model selection')
	tableauRecap = do.call( rbind, tableaurecap ) #stack list cells into a matrix
  tableauRecap = tableauRecap[rowSums(tableauRecap[, 2:4])!=0,]
  tableauRecap = tableauRecap[(tableauRecap[,1])!=Inf,]
  data = cbind(1:dim(tableauRecap)[1], tableauRecap[,2], tableauRecap[,2], tableauRecap[,1])

	require(capushe)
  modSel = capushe(data, n)
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
