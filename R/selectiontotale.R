#Return a list of outputs, for each lambda in grid: selected,Rho,Pi
selectiontotale = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau)
{
	cl = parallel::makeCluster( parallel::detectCores() / 4 )
	parallel::clusterExport(cl=cl,
		varlist=c("phiInit","rhoInit","gamInit","mini","maxi","glambda","X","Y","seuil","tau"),
		envir=environment())
	#Pour chaque lambda de la grille, on calcule les coefficients
	out = parLapply( 1:L, function(lambdaindex)
	{
		params = EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda[lambdaIndex],X,Y,tau)

		p = dim(phiInit)[1]
		m = dim(phiInit)[2]
		#selectedVariables: list where element j contains vector of selected variables in [1,m]
		selectedVariables = lapply(1:p, function(j) {
			#from boolean matrix mxk of selected variables obtain the corresponding boolean m-vector,
			#and finally return the corresponding indices
			seq_len(m)[ apply( abs(params$phi[j,,]) > seuil, 1, any ) ]
		})

		list("selected"=selectedVariables,"Rho"=params$Rho,"Pi"=params$Pi)
	})
	parallel::stopCluster(cl)
}
