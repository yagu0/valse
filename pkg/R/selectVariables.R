#' selectVariables
#'
#' It is a function which construct, for a given lambda, the sets of relevant variables.
#'
#' @param phiInit an initial estimator for phi (size: p*m*k)
#' @param rhoInit an initial estimator for rho (size: m*m*k)
#' @param piInit	an initial estimator for pi (size : k)
#' @param gamInit an initial estimator for gamma
#' @param mini		minimum number of iterations in EM algorithm
#' @param maxi		maximum number of iterations in EM algorithm
#' @param gamma	 power in the penalty
#' @param glambda grid of regularization parameters
#' @param X			 matrix of regressors
#' @param Y			 matrix of responses
#' @param thres	 threshold to consider a coefficient to be equal to 0
#' @param tau		 threshold to say that EM algorithm has converged
#'
#' @return a list of outputs, for each lambda in grid: selected,Rho,Pi
#'
#' @examples TODO
#'
#' @export
#'
selectVariables = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,
	X,Y,thresh,tau, ncores=1) #ncores==1 ==> no //
{
	if (ncores > 1)
	{
		cl = parallel::makeCluster(ncores)
		parallel::clusterExport(cl=cl,
			varlist=c("phiInit","rhoInit","gamInit","mini","maxi","glambda","X","Y","thresh","tau"),
			envir=environment())
	}

	# Calcul pour un lambda
	computeCoefs <-function(lambda)
	{
		params = EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,lambda,X,Y,tau)

		p = dim(phiInit)[1]
		m = dim(phiInit)[2]

		#selectedVariables: list where element j contains vector of selected variables in [1,m]
		selectedVariables = lapply(1:p, function(j) {
			#from boolean matrix mxk of selected variables obtain the corresponding boolean m-vector,
			#and finally return the corresponding indices
			seq_len(m)[ apply( abs(params$phi[j,,]) > thresh, 1, any ) ]
		})

		list("selected"=selectedVariables,"Rho"=params$Rho,"Pi"=params$Pi)
	}

	# Pour chaque lambda de la grille, on calcule les coefficients
	out <-
		if (ncores > 1)
			parLapply(cl, seq_along(glambda, computeCoefs)
		else
			lapply(seq_along(glambda), computeCoefs)
	if (ncores > 1)
		parallel::stopCluster(cl)
	out
}
