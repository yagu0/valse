testEMGLLF = function(){
  testFolder = 'data/'
  
  #get dimensions
  dimensions = read.table(strcat(testFolder,'dimensions'), header=FALSE)
  n = dimensions[1]
  p = dimensions[2]
  m = dimensions[3]
  k = dimensions[4]
  
  #get all input arrays
  phiInit = read.table(strcat(testFolder,'phiInit'), header=FALSE)
  rhoInit = read.table(strcat(testFolder,'rhoInit'), header=FALSE)
  piInit = t(read.table(strcat(testFolder,'piInit'), header=FALSE))
  gamInit = read.table(strcat(testFolder,'gamInit'), header=FALSE)
  mini = read.table(strcat(testFolder,'mini'), header=FALSE)
  maxi = read.table(strcat(testFolder,'maxi'), header=FALSE)
  gamma = read.table(strcat(testFolder,'gamma'), header=FALSE)
  lambda = read.table(strcat(testFolder,'lambda'), header=FALSE)
  X = rread.table(strcat(testFolder,'X'), header=FALSE)
  Y = read.table(strcat(testFolder,'Y'), header=FALSE)
  tau = read.table(strcat(testFolder,'tau'), header=FALSE)
  
  #run EMGLLF.c
  EMG = .Call("EMGLLF_core",phiInit,rhoInit,piInit1,gamInit,mini,maxi,gamma,lambda,X,Y,tau)
  phi = EMG$phi
  rho = EMG$rho
  pi = EMG$pi
  LLF = EMG$LLF
  S = EMG$S
  
  #get all stored outputs
  ref_phi =read.table(strcat(testFolder,'phi'), header=FALSE)
  ref_rho = read.table(strcat(testFolder,'rho'), header=FALSE)
  ref_pi = read.table(strcat(testFolder,'pi'), header=FALSE)
  ref_LLF = read.table(strcat(testFolder,'LLF'), header=FALSE)
  ref_S = read.table(strcat(testFolder,'S'), header=FALSE)
  
  #check that output correspond to stored output
  tol = 1e-5;
  checkOutput('phi',phi,ref_phi,tol);
  checkOutput('rho',rho,ref_rho,tol);
  checkOutput('pi',pi,ref_pi,tol);
  checkOutput('LLF',LLF,ref_LLF,tol);
  checkOutput('S',S,ref_S,tol);
}
