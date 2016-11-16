function[gridLambda] = grillelambda(phiInit,rhoInit,piInit,gamInit,X,Y,gamma,mini,maxi,tau)

	n = size(X, 1);
	[p,m,k] = size(phiInit);
	[phi,rho,pi,~,S] = EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,1,0,X,Y,tau);

	gridLambda = zeros(p,m,k);
	for j=1:p
		for mm=1:m
			gridLambda(j,mm,:) = abs(reshape(S(j,mm,:),[k,1])) ./ (n*pi.^gamma);
		end
	end
	
	gridLambda = unique(gridLambda);
	gridLambda(gridLambda()>1) = [];

end
