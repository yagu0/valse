#' EMGLLF
#'
#' Description de EMGLLF
#'
#' @param phiInit Parametre initial de moyenne renormalisé
#' @param rhoInit Parametre initial de variance renormalisé
#' @param piInit Parametre initial des proportions
#' @param gamInit Paramètre initial des probabilités a posteriori de chaque échantillon
#' @param mini Nombre minimal d'itérations dans l'algorithme EM
#' @param maxi Nombre maximal d'itérations dans l'algorithme EM
#' @param gamma Puissance des proportions dans la pénalisation pour un Lasso adaptatif
#' @param lambda Valeur du paramètre de régularisation du Lasso
#' @param X Régresseurs
#' @param Y Réponse
#' @param tau Seuil pour accepter la convergence
#'
#' @return A list ... phi,rho,pi,LLF,S,affec:
#'   phi : parametre de moyenne renormalisé, calculé par l'EM
#'   rho : parametre de variance renormalisé, calculé par l'EM
#'   pi : parametre des proportions renormalisé, calculé par l'EM
#'   LLF : log vraisemblance associée à cet échantillon, pour les valeurs estimées des paramètres
#'   S : ... affec : ...
#'
#' @export
EMGLLF <- function(phiInit, rhoInit, piInit, gamInit,
	mini, maxi, gamma, lambda, X, Y, tau, fast=TRUE)
{
	if (!fast)
	{
		# Function in R
		return (.EMGLLF_R(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,lambda,X,Y,tau))
	}

	# Function in C
	n = nrow(X) #nombre d'echantillons
	p = ncol(X) #nombre de covariables
	m = ncol(Y) #taille de Y (multivarié)
	k = length(piInit) #nombre de composantes dans le mélange
	.Call("EMGLLF",
		phiInit, rhoInit, piInit, gamInit, mini, maxi, gamma, lambda, X, Y, tau,
		phi=double(p*m*k), rho=double(m*m*k), pi=double(k), LLF=double(maxi),
			S=double(p*m*k), affec=integer(n),
		n, p, m, k,
		PACKAGE="valse")
}

# R version - slow but easy to read
.EMGLLF_R = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,lambda,X,Y,tau)
{
	# Matrix dimensions
	n = dim(X)[1]
	p = dim(phiInit)[1]
	m = dim(phiInit)[2]
	k = dim(phiInit)[3]

	# Outputs
	phi = phiInit
	rho = rhoInit
	pi = piInit
	llh = -Inf
	S = array(0, dim=c(p,m,k))

	# Algorithm variables
	gam = gamInit
	Gram2 = array(0, dim=c(p,p,k))
	ps2 = array(0, dim=c(p,m,k))
	X2 = array(0, dim=c(n,p,k))
	Y2 = array(0, dim=c(n,m,k))
	EPS = 1e-15

	for (ite in 1:maxi)
	{
		# Remember last pi,rho,phi values for exit condition in the end of loop
		Phi = phi
		Rho = rho
		Pi = pi

		# Calcul associé à Y et X
		for (r in 1:k)
		{
			for (mm in 1:m)
				Y2[,mm,r] = sqrt(gam[,r]) * Y[,mm]
			for (i in 1:n)
				X2[i,,r] = sqrt(gam[i,r]) * X[i,]
			for (mm in 1:m)
				ps2[,mm,r] = crossprod(X2[,,r],Y2[,mm,r])
			for (j in 1:p)
			{
				for (s in 1:p)
					Gram2[j,s,r] = crossprod(X2[,j,r], X2[,s,r])
			}
		}

		##########
		#Etape M #
		##########

		# Pour pi
		b = sapply( 1:k, function(r) sum(abs(phi[,,r])) )
		gam2 = colSums(gam)
		a = sum(gam %*% log(pi))

		# Tant que les props sont negatives
		kk = 0
		pi2AllPositive = FALSE
		while (!pi2AllPositive)
		{
			pi2 = pi + 0.1^kk * ((1/n)*gam2 - pi)
			pi2AllPositive = all(pi2 >= 0)
			kk = kk+1
		}

		# t(m) la plus grande valeur dans la grille O.1^k tel que ce soit décroissante ou constante
		while( kk < 1000 && -a/n + lambda * sum(pi^gamma * b) <
			-sum(gam2 * log(pi2))/n + lambda * sum(pi2^gamma * b) )
		{
			pi2 = pi + 0.1^kk * (1/n*gam2 - pi)
			kk = kk + 1
		}
		t = 0.1^kk
		pi = (pi + t*(pi2-pi)) / sum(pi + t*(pi2-pi))

		#Pour phi et rho
		for (r in 1:k)
		{
			for (mm in 1:m)
			{
				ps = 0
				for (i in 1:n)
					ps = ps + Y2[i,mm,r] * sum(X2[i,,r] * phi[,mm,r])
				nY2 = sum(Y2[,mm,r]^2)
				rho[mm,mm,r] = (ps+sqrt(ps^2+4*nY2*gam2[r])) / (2*nY2)
			}
		}

		for (r in 1:k)
		{
			for (j in 1:p)
			{
				for (mm in 1:m)
				{
					S[j,mm,r] = -rho[mm,mm,r]*ps2[j,mm,r] + sum(phi[-j,mm,r] * Gram2[j,-j,r])
					if (abs(S[j,mm,r]) <= n*lambda*(pi[r]^gamma))
						phi[j,mm,r]=0
					else if(S[j,mm,r] > n*lambda*(pi[r]^gamma))
						phi[j,mm,r] = (n*lambda*(pi[r]^gamma)-S[j,mm,r]) / Gram2[j,j,r]
					else
						phi[j,mm,r] = -(n*lambda*(pi[r]^gamma)+S[j,mm,r]) / Gram2[j,j,r]
				}
			}
		}

		##########
		#Etape E #
		##########

		# Precompute det(rho[,,r]) for r in 1...k
		detRho = sapply(1:k, function(r) det(rho[,,r]))

		sumLogLLH = 0
		for (i in 1:n)
		{
			# Update gam[,]
			sumGamI = 0
			for (r in 1:k)
			{
				gam[i,r] = pi[r]*exp(-0.5*sum((Y[i,]%*%rho[,,r]-X[i,]%*%phi[,,r])^2))*detRho[r]
				sumGamI = sumGamI + gam[i,r]
			}
			sumLogLLH = sumLogLLH + log(sumGamI) - log((2*base::pi)^(m/2))
			if (sumGamI > EPS) #else: gam[i,] is already ~=0
				gam[i,] = gam[i,] / sumGamI
		}

		sumPen = sum(pi^gamma * b)
		last_llh = llh
		llh = -sumLogLLH/n + lambda*sumPen
		dist = ifelse( ite == 1, llh, (llh-last_llh) / (1+abs(llh)) )
		Dist1 = max( (abs(phi-Phi)) / (1+abs(phi)) )
		Dist2 = max( (abs(rho-Rho)) / (1+abs(rho)) )
		Dist3 = max( (abs(pi-Pi)) / (1+abs(Pi)) )
		dist2 = max(Dist1,Dist2,Dist3)

		if (ite >= mini && (dist >= tau || dist2 >= sqrt(tau)))
			break
	}

	affec = apply(gam, 1, which.max)
	list( "phi"=phi, "rho"=rho, "pi"=pi, "llh"=llh, "S"=S, "affec"=affec )
}
