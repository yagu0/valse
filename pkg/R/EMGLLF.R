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
	mini, maxi, gamma, lambda, X, Y, tau)
{
	#TEMPORARY: use R version
	return (EMGLLF_R(phiInit, rhoInit, piInit, gamInit,mini, maxi, gamma, lambda, X, Y, tau))

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
