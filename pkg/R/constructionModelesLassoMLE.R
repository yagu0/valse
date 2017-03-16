constructionModelesLassoMLE = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,
	X,Y,seuil,tau,selected)
{
	#TODO: parameter ncores (chaque tâche peut aussi demander du parallélisme...)
	cl = parallel::makeCluster( parallel::detectCores() / 4 )
	parallel::clusterExport(cl=cl,
		varlist=c("phiInit","rhoInit","gamInit","mini","maxi","glambda","X","Y","seuil","tau"),
		envir=environment())
	#Pour chaque lambda de la grille, on calcule les coefficients
	out = parLapply( seq_along(glambda), function(lambdaindex)
	{
		n = dim(X)[1]
		p = dim(phiInit)[1]
		m = dim(phiInit)[2]
		k = dim(phiInit)[3]

		#TODO: phiInit[selected] et X[selected] sont bien sûr faux; par quoi remplacer ?
		#lambda == 0 c'est normal ? -> ED : oui, ici on calcule le maximum de vraisembance, donc on ne pénalise plus
    res = EMGLLF(phiInit[selected],rhoInit,piInit,gamInit,mini,maxi,gamma,0.,X[selected],Y,tau)

		#comment évaluer la dimension à partir du résultat et de [not]selected ?
    #dimension = ...

    #on veut calculer la vraisemblance avec toutes nos estimations
		densite = vector("double",n)
		for (r in 1:k)
		{
			delta = Y%*%rho[,,r] - (X[selected]%*%res$phi[selected,,r])
			densite = densite + pi[r] *
				det(rho[,,r])/(sqrt(2*base::pi))^m * exp(-tcrossprod(delta)/2.0)
		}
		llh = c( sum(log(densite[,lambdaIndex])), (dimension+m+1)*k-1 )
		list("phi"=res$phi, "rho"=res$rho, "pi"=res$pi, "llh" = llh)
	})
	parallel::stopCluster(cl)
	out
}
