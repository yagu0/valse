generateRunSaveTest_constructionModelesLassoMLE = function(n, p, m, k, mini, maxi, gamma, glambda, varargin){
  #set defaults for optional inputs
  optargs = c(200 15 10 3 5 10 1.0 list(0.0,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.85,0.99))
  #replace defaults by user parameters
  optargs[1:length(varargin)] = varargin
  n = optargs[1]
  p = optargs[2]
  m = optargs[3]
  k = optargs[4]
  mini = optargs[5]
  maxi = optargs[6]
  gamma = optargs[7]
  glambda = optargs[8]
  tau = 1e-6
  seuil = 1e-15
  L = length(glambda)
  
  #Generate phiInit,piInit,...
  basicInit = basicInitParameters(n, p, m, k)
  phiInit = basicInit$phiInit
  rhoInit = basicInit$rhoInit
  piInit = basicInit$piInit
  gamInit = basicInit$gamInit
  
  #Generate X and Y
  generateIOdef = generateIOdefault(n, p, m, k)
  X = generateIOdef$X
  Y = generateIOdef$Y

  A2 = array(0, dim=c(p, m+1, L))
  A1 = array(0, dim=c(p, m+1, L))
  for(i in 1:L){
    for(j in 1:p){
      A2[j, 1, i] = j
      A1[j, 1, i] = j
    }
    for(k in 1:5){
      A2[k,2,i] = k
      A1[k,2,i] = k
    }
  }

  testFolder = 'data/'
  dir.create(testFolder)
  delimiter = ' '
  
  
  #save inputs
  write(strcat(testFolder,'phiInit'), reshape(phiInit,1), delimiter)
  write(strcat(testFolder,'rhoInit'), reshape(rhoInit,1), delimiter)
  write(strcat(testFolder,'piInit'), piInit, delimiter)
  write(strcat(testFolder,'gamInit'), reshape(gamInit,1), delimiter)
  write(strcat(testFolder,'mini'), mini, delimiter)
  write(strcat(testFolder,'maxi'), maxi, delimiter)
  write(strcat(testFolder,'gamma'), gamma, delimiter)
  write(strcat(testFolder,'glambda'), glambda, delimiter)
  write(strcat(testFolder,'X'), reshape(X,1), delimiter)
  write(strcat(testFolder,'Y'), reshape(Y,1), delimiter)
  mwrite(strcat(testFolder,'seuil'), seuil, delimiter)
  write(strcat(testFolder,'tau'), tau, delimiter)
  write(strcat(testFolder,'A1'), reshape(A1,1), delimiter)
  write(strcat(testFolder,'A2'), reshape(A2,1), delimiter)
  write(strcat(testFolder,'dimensions'), [n,p,m,k,L], delimiter)
  
  construct_LME = constructionModelesLassoMLE(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,A1,A2)
  phi = construct_LME$phi
  rho = construct_LME$rho
  pi = construct_LME$pi
  lvraisemblance = construct_LME$lvraisemblance
  
  #save output
  write(strcat(testFolder,'phi'), reshape(phi,1), delimiter);
  write(strcat(testFolder,'rho'), reshape(rho,1), delimiter);
  write(strcat(testFolder,'pi'), reshape(pi,1), delimiter);
  write(strcat(testFolder,'lvraisemblance'), reshape(lvraisemblance,1), delimiter);
}
