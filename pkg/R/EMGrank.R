#' EMGrank
#'
#' Description de EMGrank
#'
#' @param Pi ...
#'
#' @return ...
#'
#' @examples
#' ...
#' ...
#' @export
EMGrank <- function(Pi, Rho, mini, maxi, X, Y, tau, rank)
{
	.Call("EMGrank", Pi, Rho, mini, maxi, X, Y, tau, rank, PACKAGE="valse")
}
