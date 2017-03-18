#' EMGLLF
#'
#' Description de EMGLLF
#'
#' @param phiInit ...
#'
#' @return ...
#'
#' @examples
#' ...
#' ...
#' @export
EMGLLF <- function(phiInit, rhoInit, piInit, gamInit,
	mini, maxi, gamma, lambda, X, Y, tau)
{
	.Call("EMGLLF", phiInit, rhoInit, piInit, gamInit,
		mini, maxi, gamma, lambda, X, Y, tau, PACKAGE="valse")
}
