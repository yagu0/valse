function[phiInit,rhoInit,piInit,gamInit] = basicInitParameters(n,p,m,k)

	phiInit = zeros(p,m,k);
	
	piInit = (1.0/k) * ones(1,k);
	
	rhoInit = zeros(m,m,k);
	for r=1:k
		rhoInit(:,:,r) = eye(m,m);
	end
	
	gamInit = 0.1 * ones(n,k);
	R = random('unid',k,n,1);
	for i=1:n
		gamInit(i,R(i)) = 0.9;
	end
	gamInit = gamInit / (sum(gamInit(1,:)));

end
