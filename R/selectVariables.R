#' selectVaribles
#' It is a function which construct, for a given lambda, the sets of
#' relevant variables and irrelevant variables.
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
#' @return
#' @export
#'
#' @examples
selectVariables <- function(phiInit,rhoInit,piInit,gamInit,
	mini,maxi,gamma,glambda,X,Y,thres,tau)
{
	dimphi <- dim(phiInit)
	p <- dimPhi[1]
	m <- dimPhi[2]
	k <- dimPhi[3]
	L <- length(glambda);
	A1 <- array(0, dim <- c(p,m+1,L))
	A2 <- array(0, dim <- c(p,m+1,L))
	Rho <- array(0, dim <- c(m,m,k,L))
	Pi <- array(0, dim <- c(k,L));

	# For every lambda in gridLambda, comutation of the coefficients
	for (lambdaIndex in c(1:L))
	{
		Res <- EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,
			gamma,glambda[lambdaIndex],X,Y,tau);
		phi <- Res$phi
		rho <- Res$rho
		pi <- Res$pi

		# If a coefficient is larger than the threshold, we keep it
		selectedVariables <- array(0, dim = c(p,m))
		discardedVariables <- array(0, dim = c(p,m))
		atLeastOneSelectedVariable <- false
		for (j in c(1:p))
		{
			cpt <- 1
			cpt2 <-1
			for (mm in c(1:m))
			{
				if (max(abs(phi[j,mm,])) > thres)
				{
					selectedVariables[j,cpt] <- mm
					cpt <- cpt+1
					atLeastOneSelectedVariable <- true
				} else
				{
					discardedVariables[j,cpt2] <- mm
					cpt2 <- cpt2+1
				}
			}
		}

		# If no coefficients have been selected, we provide the zero matrix
		# We delete zero coefficients: vec = indices of zero values
		if (atLeastOneSelectedVariable)
		{
			vec <- c()
			for (j in c(1:p))
			{
				if (selectedVariables(j,1) != 0)
					vec <- c(vec,j)
				# Else ( NOTE: [auder] else ?! TODO: explain? )
				# we provide the indices of relevant coefficients
				A1[,1,lambdaIndex] <- c(vec,rep(0,p-length(vec)))
				A1[1:length(vec),2:(m+1),lambdaIndex] <- selectedVariables[vec,]
				A2[,1,lambdaIndex] <- 1:p
				A2[,2:(m+1),lambdaIndex] <- discardedVariables
				Rho[,,,lambdaIndex] <- rho
				Pi[,lambdaIndex] <- pi
			}
		}
	}

	return(res = list(A1 = A1, A2 = A2 , Rho = Rho, Pi = Pi))
}
