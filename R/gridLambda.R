gridLambda = function(phiInit, rhoInit, piInit, gamInit, X, Y, gamma, mini, maxi, tau)
{
	n = nrow(X)
	p = dim(phiInit)[1]
	m = dim(phiInit)[2]
	k = dim(phiInit)[3]

	list_EMG = .Call("EMGLLF",phiInit,rhoInit,piInit,gamInit,mini,maxi,1,0,X,Y,tau)

	grid = array(0, dim=c(p,m,k))
	for (i in 1:p)
	{
		for (j in 1:m)
			grid[i,j,] = abs(list_EMG$S[i,j,]) / (n*list_EMG$pi^gamma)
	}
	grid = unique(grid)
	grid = grid[grid <=1]

	return(grid)
}
