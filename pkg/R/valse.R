#' Main function
#'
#' @param X matrix of covariates (of size n*p)
#' @param Y matrix of responses (of size n*m)
#' @param procedure among 'LassoMLE' or 'LassoRank'
#' @param selecMod method to select a model among 'DDSE', 'DJump', 'BIC' or 'AIC'
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
valse = function(X,Y,procedure = 'LassoMLE',selecMod = 'DDSE',gamma = 1,mini = 10,
                 maxi = 50,eps = 1e-4,kmin = 2,kmax = 3,
                 rang.min = 1,rang.max = 10) {
  ##################################
  #core workflow: compute all models
  ##################################
  
  p = dim(X)[2]
  m = dim(Y)[2]
  n = dim(X)[1]
  
  model = list()
  tableauRecap = array(0, dim=c(1000,4))
  cpt = 0
  print("main loop: over all k and all lambda")
  
  for (k in kmin:kmax){
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
    source('~/valse/pkg/R/gridLambda.R')
    grid_lambda <<- gridLambda(phiInit, rhoInit, piInit, gamInit, X, Y, gamma, mini, maxi, eps)
    
    # if (length(grid_lambda)>50){
    #   grid_lambda = grid_lambda[seq(1, length(grid_lambda), length.out = 50)]
    # }
    print("Compute relevant parameters")
    #select variables according to each regularization parameter
    #from the grid: A1 corresponding to selected variables, and
    #A2 corresponding to unselected variables.
    
    params = selectiontotale(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,grid_lambda,X,Y,1e-8,eps)
    #params2 = selectVariables(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,grid_lambda[seq(1,length(grid_lambda), by=3)],X,Y,1e-8,eps)
    ## etrange : params et params 2 sont différents ...
    
    selected <<- params$selected
    Rho <<- params$Rho
    Pi <<- params$Pi
    
    if (procedure == 'LassoMLE') {
      print('run the procedure Lasso-MLE')
      #compute parameter estimations, with the Maximum Likelihood
      #Estimator, restricted on selected variables.
      model[[k]] = constructionModelesLassoMLE(phiInit, rhoInit,piInit,gamInit,mini,maxi,gamma,X,Y,thresh,eps,selected)
      llh = matrix(ncol = 2)
      for (l in seq_along(model[[k]])){
        llh = rbind(llh, model[[k]][[l]]$llh)
      }
      LLH = llh[-1,1]
      D = llh[-1,2]
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
    tableauRecap[(cpt+1):(cpt+length(model[[k]])), ] = matrix(c(LLH, D, rep(k, length(model[[k]])), 1:length(model[[k]])), ncol = 4)
    cpt = cpt+length(model[[k]])
  }
  print('Model selection')
  tableauRecap = tableauRecap[rowSums(tableauRecap[, 2:4])!=0,]
  tableauRecap = tableauRecap[(tableauRecap[,1])!=Inf,]
  data = cbind(1:dim(tableauRecap)[1], tableauRecap[,2], tableauRecap[,2], tableauRecap[,1])
  require(capushe)
  modSel = capushe(data, n)
  if (selecMod == 'DDSE') {
    indModSel = as.numeric(modSel@DDSE@model)
  } else if (selecMod == 'Djump') {
    indModSel = as.numeric(modSel@Djump@model)
  } else if (selecMod == 'BIC') {
    indModSel = modSel@BIC_capushe$model
  } else if (selecMod == 'AIC') {
    indModSel = modSel@AIC_capushe$model
  }
  return(model[[tableauRecap[indModSel,3]]][[tableauRecap[indModSel,4]]])
}
