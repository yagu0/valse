#' constructionModelesLassoMLE
#'
#' TODO: description
#'
#' @param ...
#'
#' @return ...
#'
#' export
constructionModelesLassoMLE = function(phiInit, rhoInit, piInit, gamInit, mini, maxi,
	gamma, X, Y, seuil, tau, selected, ncores=3, verbose=FALSE)
{
  if (ncores > 1)
	{
    cl = parallel::makeCluster(ncores)
    parallel::clusterExport( cl, envir=environment(),
			varlist=c("phiInit","rhoInit","gamInit","mini","maxi","gamma","X","Y","seuil",
			"tau","selected","ncores","verbose") )
	}

	# Individual model computation
	computeAtLambda <- function(lambda)
	{
		if (ncores > 1)
			require("valse") #// nodes start with an ampty environment

    if (verbose)
			print(paste("Computations for lambda=",lambda))

		n = dim(X)[1]
		p = dim(phiInit)[1]
		m = dim(phiInit)[2]
		k = dim(phiInit)[3]

		sel.lambda = selected[[lambda]]
#		col.sel = which(colSums(sel.lambda)!=0) #if boolean matrix
		col.sel <- which( sapply(sel.lambda,length) > 0 ) #if list of selected vars

		if (length(col.sel) == 0)
			return (NULL)

		# lambda == 0 because we compute the EMV: no penalization here
		res = EMGLLF(phiInit[col.sel,,],rhoInit,piInit,gamInit,mini,maxi,gamma,0,
			X[,col.sel],Y,tau)
		
		# Eval dimension from the result + selected
		phiLambda2 = res_EM$phi
		rhoLambda = res_EM$rho
		piLambda = res_EM$pi
    phiLambda = array(0, dim = c(p,m,k))
		for (j in seq_along(col.sel))
			phiLambda[col.sel[j],,] = phiLambda2[j,,]

		dimension = 0
		for (j in 1:p)
		{
			b = setdiff(1:m, sel.lambda[,j])
			if (length(b) > 0)
				phiLambda[j,b,] = 0.0
			dimension = dimension + sum(sel.lambda[,j]!=0)
		}

		# on veut calculer la vraisemblance avec toutes nos estimations
		densite = vector("double",n)
		for (r in 1:k)
		{
			delta = Y%*%rhoLambda[,,r] - (X[, col.sel]%*%phiLambda[col.sel,,r])
			densite = densite + piLambda[r] *
				det(rhoLambda[,,r])/(sqrt(2*base::pi))^m * exp(-tcrossprod(delta)/2.0)
		}
		llhLambda = c( sum(log(densite)), (dimension+m+1)*k-1 )
		list("phi"= phiLambda, "rho"= rhoLambda, "pi"= piLambda, "llh" = llhLambda)
	}

	#Pour chaque lambda de la grille, on calcule les coefficients
  out =
		if (ncores > 1)
			parLapply(cl, glambda, computeAtLambda)
		else
			lapply(glambda, computeAtLambda)

	if (ncores > 1)
    parallel::stopCluster(cl)

	out
}
