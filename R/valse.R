#' Main function
#'
#' @param X matrix of covariates (of size n*p)
#' @param Y matrix of responses (of size n*m)
#' @param procedure among 'LassoMLE' or 'LassoRank'
#' @param selecMod method to select a model among 'SlopeHeuristic', 'BIC', 'AIC'
#' @param gamma integer for the power in the penaly, by default = 1
#' @param mini integer, minimum number of iterations in the EM algorithm, by default = 10
#' @param maxi integer, maximum number of iterations in the EM algorithm, by default = 100
#' @param eps real, threshold to say the EM algorithm converges, by default = 1e-4
#' @param kmin integer, minimum number of clusters, by default = 2
#' @param kmax integer, maximum number of clusters, by default = 10
#' @param rang.min integer, minimum rank in the low rank procedure, by default = 1
#' @param rang.max integer, maximum rank in the
#' @return a list with estimators of parameters
#' @export
#-----------------------------------------------------------------------
valse = function(X,Y,procedure,selecMod,gamma = 1,mini = 10,
                 maxi = 100,eps = 1e-4,kmin = 2,kmax = 10,
                 rang.min = 1,rang.max = 10) {
  ##################################
  #core workflow: compute all models
  ##################################
  
  p = dim(phiInit)[1]
  m = dim(phiInit)[2]
  
  print("main loop: over all k and all lambda")
  for (k in kmin:kmax)
  {
    print(k)
    
    print("Parameters initialization")
    #smallEM initializes parameters by k-means and regression model in each component,
    #doing this 20 times, and keeping the values maximizing the likelihood after 10
    #iterations of the EM algorithm.
    init = initSmallEM(k, X, Y)
    phiInit <<- init$phiInit
    rhoInit <<- init$rhoInit
    piInit	<<- init$piInit
    gamInit <<- init$gamInit
    
    gridLambda <<- gridLambda(phiInit, rhoInit, piInit, tauInit, X, Y, gamma, mini, maxi, eps)
    
    print("Compute relevant parameters")
    #select variables according to each regularization parameter
    #from the grid: A1 corresponding to selected variables, and
    #A2 corresponding to unselected variables.
    params = selectiontotale(phiInit,rhoInit,piInit,tauInit,
                             mini,maxi,gamma,gridLambda,
                             X,Y,thresh,eps)
    A1 <<- params$A1
    A2 <<- params$A2
    Rho <<- params$Rho
    Pi <<- params$Pi
    
    if (procedure == 'LassoMLE') {
      print('run the procedure Lasso-MLE')
      #compute parameter estimations, with the Maximum Likelihood
      #Estimator, restricted on selected variables.
      model = constructionModelesLassoMLE(
        phiInit, rhoInit,piInit,tauInit,mini,maxi,
        gamma,gridLambda,X,Y,thresh,eps,A1,A2)
      ################################################
      ### Regarder la SUITE
      r1 = runProcedure1()
      Phi2 = Phi
      Rho2 = Rho
      Pi2 = Pi
      
      if (is.null(dim(Phi2)))
        #test was: size(Phi2) == 0
      {
        Phi[, , 1:k] <<- r1$phi
        Rho[, , 1:k] <<- r1$rho
        Pi[1:k,] <<- r1$pi
      } else
      {
        Phi <<-
          array(0., dim = c(p, m, kmax, dim(Phi2)[4] + dim(r1$phi)[4]))
        Phi[, , 1:(dim(Phi2)[3]), 1:(dim(Phi2)[4])] <<- Phi2
        Phi[, , 1:k, dim(Phi2)[4] + 1] <<- r1$phi
        Rho <<-
          array(0., dim = c(m, m, kmax, dim(Rho2)[4] + dim(r1$rho)[4]))
        Rho[, , 1:(dim(Rho2)[3]), 1:(dim(Rho2)[4])] <<- Rho2
        Rho[, , 1:k, dim(Rho2)[4] + 1] <<- r1$rho
        Pi <<- array(0., dim = c(kmax, dim(Pi2)[2] + dim(r1$pi)[2]))
        Pi[1:nrow(Pi2), 1:ncol(Pi2)] <<- Pi2
        Pi[1:k, ncol(Pi2) + 1] <<- r1$pi
      }
    } else {
      print('run the procedure Lasso-Rank')
      #compute parameter estimations, with the Low Rank
      #Estimator, restricted on selected variables.
      model = constructionModelesLassoRank(Pi, Rho, mini, maxi, X, Y, eps,
                                           A1, rank.min, rank.max)
      
      ################################################
      ### Regarder la SUITE  
      phi = runProcedure2()$phi
      Phi2 = Phi
      if (dim(Phi2)[1] == 0)
      {
        Phi[, , 1:k,] <<- phi
      } else
      {
        Phi <<- array(0, dim = c(p, m, kmax, dim(Phi2)[4] + dim(phi)[4]))
        Phi[, , 1:(dim(Phi2)[3]), 1:(dim(Phi2)[4])] <<- Phi2
        Phi[, , 1:k,-(1:(dim(Phi2)[4]))] <<- phi
      }
    }
  }
  print('Model selection')
  if (selecMod == 'SlopeHeuristic') {
    
  } else if (selecMod == 'BIC') {
    
  } else if (selecMod == 'AIC') {
    
  }
}
