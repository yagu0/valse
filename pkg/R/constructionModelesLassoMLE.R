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
	gamma, X, Y, thresh, tau, S, ncores=3, artefact = 1e3, verbose=FALSE)
{
	if (ncores > 1)
	{
		cl = parallel::makeCluster(ncores, outfile='')
		parallel::clusterExport( cl, envir=environment(),
			varlist=c("phiInit","rhoInit","gamInit","mini","maxi","gamma","X","Y","thresh",
			"tau","S","ncores","verbose") )
	}

	# Individual model computation
	computeAtLambda <- function(lambda)
	{
		if (ncores > 1)
			require("valse") #nodes start with an empty environment

		if (verbose)
			print(paste("Computations for lambda=",lambda))

		n = dim(X)[1]
		p = dim(phiInit)[1]
		m = dim(phiInit)[2]
		k = dim(phiInit)[3]

		sel.lambda = S[[lambda]]$selected
#		col.sel = which(colSums(sel.lambda)!=0) #if boolean matrix
		col.sel <- which( sapply(sel.lambda,length) > 0 ) #if list of selected vars

		if (length(col.sel) == 0)
			return (NULL)

		# lambda == 0 because we compute the EMV: no penalization here
		res = EMGLLF(phiInit[col.sel,,],rhoInit,piInit,gamInit,mini,maxi,gamma,0,
			X[,col.sel],Y,tau)
		
		# Eval dimension from the result + selected
		phiLambda2 = res$phi
		rhoLambda = res$rho
		piLambda = res$pi
		phiLambda = array(0, dim = c(p,m,k))
		for (j in seq_along(col.sel))
			phiLambda[col.sel[j],,] = phiLambda2[j,,]
		dimension = length(unlist(sel.lambda))

		# Computation of the loglikelihood
		densite = vector("double",n)
		for (r in 1:k)
		{
			delta = (Y%*%rhoLambda[,,r] - (X[, col.sel]%*%phiLambda[col.sel,,r]))/artefact
			print(max(delta))
			densite = densite + piLambda[r] *
				det(rhoLambda[,,r])/(sqrt(2*base::pi))^m * exp(-tcrossprod(delta)/2.0)
		}
		llhLambda = c( sum(artefact^2 * log(densite)), (dimension+m+1)*k-1 )
		list("phi"= phiLambda, "rho"= rhoLambda, "pi"= piLambda, "llh" = llhLambda)
	}

	# For each lambda, computation of the parameters
	out =
		if (ncores > 1)
			parLapply(cl, 1:length(S), computeAtLambda)
		else
			lapply(1:length(S), computeAtLambda)

	if (ncores > 1)
		parallel::stopCluster(cl)

	out
}
