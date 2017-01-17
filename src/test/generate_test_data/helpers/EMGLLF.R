EMGLLF = function(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,lambda,X,Y,tau){
  #matrix dimensions
  n = dim(X)[1]
  p = dim[phiInit][1]
  m = dim[phiInit][2]
  k = dim[phiInit][3]
  
  #init outputs
  phi = phiInit
  rho = rhoInit
  Pi = piInit
  LLF = rep(0, maxi)
  S = array(0, dim=c(p,m,k))
  
  
  gam = gamInit
  Gram2 = array(0, dim=c(p,p,k))
  ps2 = array(0, dim=c(p,m,k))
  b = rep(0, k)
  pen = matrix(0, maxi, k)
  X2 = array(0, dim=c(n,p,k))
  Y2 = array(0, dim=c(p,m,k))
  dist = 0
  dist2 = 0
  ite = 1
  Pi2 = rep(0, k)
  ps = matrix(0, m,k)
  nY2 = matrix(0, m,k)
  ps1 = array(0, dim=c(n,m,k))
  nY21 = array(0, dim=c(n,m,k))
  Gam = matrix(0, n,k)
  EPS = 1E-15
  
  while(ite <= mini || (ite<= maxi && (dist>= tau || dist2 >= sqrt(tau)))){
    Phi = phi
    Rho = rho
    PI = Pi
    #calcul associé à Y et X
    for(r in 1:k){
      for(mm in 1:m){
        Y2[,mm,r] = sqrt(gam[,r]) .* Y[,mm]
      }
      for(i in 1:n){
        X2[i,,r] = X[i,] .* sqrt(gam[i,r])
      }
      for(mm in 1:m){
        ps2[,mm,r] = crossprod(X2[,,r],Y2[,mm,r])
      }
      for(j in 1:p){
        for(s in 1:p){
          Gram2[j,s,r] = tcrossprod(X2[,j,r], X2[,s,r])
        }
      }
    }
    
    ##########
    #Etape M #
    ##########
    
    #pour pi
    for(r in 1:k){
      b[r] = sum(sum(abs(phi[,,r])))
    }
    gam2 = sum(gam[1,])  #BIG DOUTE
    a = sum(gam*t(log(Pi)))
    
    #tant que les props sont negatives
    kk = 0
    pi2AllPositive = FALSE
    while(pi2AllPositive == FALSE){
      pi2 = pi + 0.1^kk * ((1/n)*gam2 - pi)
      pi2AllPositive = TRUE
      for(r in 1:k){
        if(pi2[r] < 0){
          pi2AllPositive = false;
          break
        }
      }
      kk = kk+1
    }
    
    #t[m]la plus grande valeur dans la grille O.1^k tel que ce soit
    #décroissante ou constante
    while((-1/n*a+lambda*((pi.^gamma)*b))<(-1/n*gam2*t(log(pi2))+lambda.*(pi2.^gamma)*b) && kk<1000){
      pi2 = pi+0.1^kk*(1/n*gam2-pi)
      kk = kk+1
    }
    t = 0.1^(kk)
    pi = (pi+t*(pi2-pi)) / sum(pi+t*(pi2-pi))
    
    #Pour phi et rho
    for(r in 1:k){
      for(mm in 1:m){
        for(i in 1:n){
          ps1[i,mm,r] = Y2[i,mm,r] * dot(X2(i,:,r), phi(:,mm,r))
          nY21[i,mm,r] = (Y2[i,mm,r])^2
        }
        ps[mm,r] = sum(ps1(:,mm,r));
        nY2[mm,r] = sum(nY21(:,mm,r));
        rho[mm,mm,r] = ((ps[mm,r]+sqrt(ps[mm,r]^2+4*nY2[mm,r]*(gam2[r])))/(2*nY2[mm,r]))
      }
    }
    for(r in 1:k){
      for(j in 1:p){
        for(mm in 1:m){
          S[j,mm,r] = -rho[mm,mm,r]*ps2[j,mm,r] + dot(phi[1:j-1,mm,r],Gram2[j,1:j-1,r])  + dot(phi[j+1:p,mm,r],Gram2[j,j+1:p,r])
          if(abs(S(j,mm,r)) <= n*lambda*(pi(r)^gamma))
            phi[j,mm,r]=0
          else{
            if(S[j,mm,r]> n*lambda*(Pi[r]^gamma))
              phi[j,mm,r] = (n*lambda*(Pi[r]^gamma)-S[j,mm,r])/Gram2[j,j,r]
          else
            phi[j,mm,r] = -(n*lambda*(Pi[r]^gamma)+S[j,mm,r])/Gram2[j,j,r]
          }
        }
      }
    }
    
    ##########
    #Etape E #
    ##########
    sumLogLLF2 = 0
    for(i in 1:n){
      #precompute dot products to numerically adjust their values
      dotProducts = rep(0,k)
      for(r in 1:k){
        dotProducts[r] = tcrossprod(Y[i,]%*%rho[,,r]-X[i,]%*%phi[,,r])
      }
      shift = 0.5*min(dotProducts)
    
      #compute Gam(:,:) using shift determined above
      sumLLF1 = 0.0;
      for(r in 1:k){
        Gam[i,r] = Pi[r]*det(rho[,,r])*exp(-0.5*dotProducts[r] + shift)
        sumLLF1 = sumLLF1 + Gam[i,r]/(2*pi)^(m/2)
      }
      sumLogLLF2 = sumLogLLF2 + log(sumLLF1)
      sumGamI = sum(Gam[i,])
      if(sumGamI > EPS)
        gam[i,] = Gam[i,] / sumGamI
      else
        gam[i,] = rep(0,k) 
    }
    
    
    sumPen = 0
    for(r in 1:k){
      sumPen = sumPen + Pi[r].^gamma^b[r]
    }
    LLF[ite] = -(1/n)*sumLogLLF2 + lambda*sumPen
    
    if(ite == 1)
      dist = LLF[ite]
    else
      dist = (LLF[ite]-LLF[ite-1])/(1+abs(LLF[ite]))
    
    Dist1=max(max(max((abs(phi-Phi))./(1+abs(phi)))))
    Dist2=max(max(max((abs(rho-Rho))./(1+abs(rho)))))
    Dist3=max(max((abs(Pi-PI))./(1+abs(PI))))
    dist2=max([Dist1,Dist2,Dist3])
    
    ite=ite+1
  }
    
  Pi = transpose(Pi)
  return(list(phi=phi, rho=rho, Pi=Pi, LLF=LLF, S=S))
}