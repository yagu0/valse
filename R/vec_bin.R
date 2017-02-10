#' A function needed in initSmallEM
#'
#' @param X vector with integer values
#' @param r integer
#'
#' @return a list with Z (a binary vector of size the size of X) and indices where Z is equal to 1
vec_bin = function(X,r)
{
	Z = rep(0,length(X))
	indice = c()
	j = 1
	for (i in 1:length(X))
	{
		if(X[i] == r)
		{
			Z[i] = 1
			indice[j] = i
			j=j+1
		} else
			Z[i] = 0
	}
	return (list(Z=Z,indice=indice))
}
