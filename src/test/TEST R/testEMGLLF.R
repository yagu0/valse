testEMGLLF = function(n,p,m,k, phiInit, rhoInit, gamInit, mini, maxi, gamma, lambda, X, Y, tau, ref_phi, ref_rho, ref_pi, ref_LLf, ref_S){
  testFolder = 'data/'
  delimiter = '\n'
  
  EMG = .Call("EMGLLF_core",phiInit,rhoInit,piInit1,gamInit,mini,maxi,gamma,lambda,X,Y,tau)
  phi = EMG$phi
  rho = EMG$rho
  pi = EMG$pi
  LLF = EMG$LLF
  S = EMG$S
  
  
  tol = 1e-5;
  checkOutput('phi',phi,ref_phi,tol);
  checkOutput('rho',rho,ref_rho,tol);
  checkOutput('pi',pi,ref_pi,tol);
  checkOutput('LLF',LLF,ref_LLF,tol);
  checkOutput('S',S,ref_S,tol);
}
