#' EMGrank
#'
#' Description de EMGrank
#'
#' @param phiInit ...
#' @param Pi Parametre de proportion
#' @param Rho Parametre initial de variance renormalisé
#' @param mini Nombre minimal d'itérations dans l'algorithme EM
#' @param maxi Nombre maximal d'itérations dans l'algorithme EM
#' @param X Régresseurs
#' @param Y Réponse
#' @param tau Seuil pour accepter la convergence
#' @param rank Vecteur des rangs possibles
#'
#' @return A list ...
#'   phi : parametre de moyenne renormalisé, calculé par l'EM
#'   LLF : log vraisemblance associé à cet échantillon, pour les valeurs estimées des paramètres
#'
#' @export
EMGrank <- function(Pi, Rho, mini, maxi, X, Y, tau, rank)
{
	#TEMPORARY: use R version
	return (EMGrank_R(Pi, Rho, mini, maxi, X, Y, tau, rank))

	n = nrow(X) #nombre d'echantillons
	p = ncol(X) #nombre de covariables
	m = ncol(Y) #taille de Y (multivarié)
	k = length(Pi) #nombre de composantes dans le mélange
	.Call("EMGrank",
		Pi, Rho, mini, maxi, X, Y, tau, rank,
		phi=double(p*m*k), LLF=double(1),
		n, p, m, k,
		PACKAGE="valse")
}
