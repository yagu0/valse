constructionModelesLassoMLE = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,
                                       X,Y,seuil,tau,selected, parallel = FALSE)
{
  if (parallel) {
    #TODO: parameter ncores (chaque tâche peut aussi demander du parallélisme...)
    cl = parallel::makeCluster( parallel::detectCores() / 4 )
    parallel::clusterExport(cl=cl,
                            varlist=c("phiInit","rhoInit","gamInit","mini","maxi","X","Y","seuil","tau"),
                            envir=environment())
    #Pour chaque lambda de la grille, on calcule les coefficients
    out = parLapply( seq_along(glambda), function(lambda)
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
      llh = c( sum(log(densite[,lambda])), (dimension+m+1)*k-1 )
      list("phi"=res$phi, "rho"=res$rho, "pi"=res$pi, "llh" = llh)
    })
    parallel::stopCluster(cl)
    out
  }
  else {
    #Pour chaque lambda de la grille, on calcule les coefficients
    n = dim(X)[1]
    p = dim(phiInit)[1]
    m = dim(phiInit)[2]
    k = dim(phiInit)[3]
    L = length(selected)
    phi = list()
    phiLambda = array(0, dim = c(p,m,k))
    rho = list()
    pi = list()
    llh = list()
    
    out = lapply( seq_along(selected), function(lambda)
    {
      print(lambda)
      sel.lambda = selected[[lambda]]
      col.sel = which(colSums(sel.lambda)!=0)
      if (length(col.sel)>0){
        res_EM = EMGLLF(phiInit[col.sel,,],rhoInit,piInit,gamInit,mini,maxi,gamma,0.,X[,col.sel],Y,tau)
        phiLambda2 = res_EM$phi
        rhoLambda = res_EM$rho
        piLambda = res_EM$pi
        for (j in 1:length(col.sel)){
          phiLambda[col.sel[j],,] = phiLambda2[j,,]
        }
        
        dimension = 0
        for (j in 1:p){
          b = setdiff(1:m, sel.lambda[,j])
          if (length(b) > 0){
            phiLambda[j,b,] = 0.0
          }
          dimension = dimension + sum(sel.lambda[,j]!=0)
        }
        
        #on veut calculer la vraisemblance avec toutes nos estimations
        densite = vector("double",n)
        for (r in 1:k)
        {
          delta = Y%*%rhoLambda[,,r] - (X[, col.sel]%*%phiLambda[col.sel,,r])
          densite = densite + piLambda[r] *
            det(rhoLambda[,,r])/(sqrt(2*base::pi))^m * exp(-tcrossprod(delta)/2.0)
        }
        llhLambda = c( sum(log(densite)), (dimension+m+1)*k-1 )
        list("phi"= phiLambda, "rho"= rhoLambda, "pi"= piLambda, "llh" = llhLambda)
      }
    }
    )
    return(out)
  }
}
