covariance = function(p,a)
{
	A = matrix(a, p,p)
	for (i in 1:p)
		A[i,] = A[i,]^abs(i-(1:p))

	return (A)
}
