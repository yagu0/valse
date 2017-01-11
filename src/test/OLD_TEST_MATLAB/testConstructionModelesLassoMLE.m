function[] = testConstructionModelesLassoMLE()

	testFolder = 'data/';
	delimiter = '\n';

	%get dimensions
	dimensions = dlmread(strcat(testFolder,'dimensions'), delimiter);
	n = dimensions(1);
	p = dimensions(2);
	m = dimensions(3);
	k = dimensions(4);
	L = dimensions(5);

	%get all input arrays
	phiInit = reshape(dlmread(strcat(testFolder,'phiInit'), delimiter), p, m, k);
	rhoInit = reshape(dlmread(strcat(testFolder,'rhoInit'), delimiter), m, m, k);
	piInit = transpose(dlmread(strcat(testFolder,'piInit'), delimiter));
	gamInit = reshape(dlmread(strcat(testFolder,'gamInit'), delimiter), n, k);
	mini = int64(dlmread(strcat(testFolder,'mini'), delimiter));
	maxi = int64(dlmread(strcat(testFolder,'maxi'), delimiter));
	gamma = dlmread(strcat(testFolder,'gamma'), delimiter);
	glambda = dlmread(strcat(testFolder,'glambda'), delimiter);
	X = reshape(dlmread(strcat(testFolder,'X'), delimiter), n, p);
	Y = reshape(dlmread(strcat(testFolder,'Y'), delimiter), n, m);
	seuil = dlmread(strcat(testFolder,'seuil'), delimiter);
	tau = dlmread(strcat(testFolder,'tau'), delimiter);
	A1 = int64(reshape(dlmread(strcat(testFolder,'A1'), delimiter), p, m+1, L));
	A2 = int64(reshape(dlmread(strcat(testFolder,'A2'), delimiter), p, m+1, L));

	%run constructionModelesLassoMLE.m
	[phi,rho,pi,lvraisemblance] = constructionModelesLassoMLE(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,A1,A2);

	%get all stored outputs
	ref_phi = reshape(dlmread(strcat(testFolder,'phi'), delimiter), p, m, k, L);
	ref_rho = reshape(dlmread(strcat(testFolder,'rho'), delimiter), m, m, k, L);
	ref_pi = reshape(dlmread(strcat(testFolder,'pi'), delimiter), k, L);
	ref_lvraisemblance = reshape(dlmread(strcat(testFolder,'lvraisemblance'), delimiter), L, 2);

	%check that output correspond to stored output
	tol = 1e-5;
	checkOutput('phi',phi,ref_phi,tol);
	checkOutput('rho',rho,ref_rho,tol);
	checkOutput('pi',pi,ref_pi,tol);
	checkOutput('lvraisemblance',lvraisemblance,ref_lvraisemblance,tol);

end
