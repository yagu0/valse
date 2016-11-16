function[phi,rho,pi,lvraisemblance] = constructionModelesLassoMLE(...
	phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,A1,A2)

	PI = 4.0 * atan(1.0);
	
	%get matrix sizes
	n = size(X, 1);
	[p,m,k] = size(phiInit);
	L = length(glambda);
	
	%output parameters
	phi = zeros(p,m,k,L);
	rho = zeros(m,m,k,L);
	pi = zeros(k,L);
	lvraisemblance = zeros(L,2);
	
	for lambdaIndex=1:L
		% Procedure Lasso-MLE  
		a = A1(:,1,lambdaIndex);
		a(a==0) = [];
		if length(a) == 0
			continue;
		end
		[phiLambda,rhoLambda,piLambda,~,~] = EMGLLF(...
			phiInit(a,:,:),rhoInit,piInit,gamInit,mini,maxi,gamma,0,X(:,a),Y,tau);
		
		for j=1:length(a)
			phi(a(j),:,:,lambdaIndex) = phiLambda(j,:,:);
		end
		rho(:,:,:,lambdaIndex) = rhoLambda;
		pi(:,lambdaIndex) = piLambda;
		
		dimension = 0;
		for j=1:p
			b = A2(j,2:end,lambdaIndex);
			b(b==0) = [];
			if length(b) > 0
				phi(A2(j,1,lambdaIndex),b,:,lambdaIndex) = 0.0;
			end
			c = A1(j,2:end,lambdaIndex);
			c(c==0) = [];
			dimension = dimension + length(c);
		end
		
		%on veut calculer l'EMV avec toutes nos estimations
		densite = zeros(n,L);
		for i=1:n
			for r=1:k
				delta = Y(i,:)*rho(:,:,r,lambdaIndex) - (X(i,a)*(phi(a,:,r,lambdaIndex)));
				densite(i,lambdaIndex) = densite(i,lambdaIndex) +...
					pi(r,lambdaIndex)*det(rho(:,:,r,lambdaIndex))/(sqrt(2*PI))^m*exp(-dot(delta,delta)/2.0);
			end
		end
		lvraisemblance(lambdaIndex,1) = sum(log(densite(:,lambdaIndex)));
		lvraisemblance(lambdaIndex,2) = (dimension+m+1)*k-1;
	end

end
