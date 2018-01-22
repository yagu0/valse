# ...
gdet <- function(M)
{
	if (is.matrix(M))
		return (det(M))
	return (M[1]) #numeric, double
}
