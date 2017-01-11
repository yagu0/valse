function[] = testEMGLLF()

	testFolder = 'data/';
	delimiter = '\n';

	%get dimensions
	dimensions = dlmread(strcat(testFolder,'dimensions'), delimiter);
	n = dimensions(1);
	p = dimensions(2);
	m = dimensions(3);
	k = dimensions(4);

	%get all input arrays
	phiInit = reshape(dlmread(strcat(testFolder,'phiInit'), delimiter), p, m, k);
	rhoInit = reshape(dlmread(strcat(testFolder,'rhoInit'), delimiter), m, m, k);
	piInit = transpose(dlmread(strcat(testFolder,'piInit'), delimiter));
	gamInit = reshape(dlmread(strcat(testFolder,'gamInit'), delimiter), n, k);
	mini = int64(dlmread(strcat(testFolder,'mini'), delimiter));
	maxi = int64(dlmread(strcat(testFolder,'maxi'), delimiter));
	gamma = dlmread(strcat(testFolder,'gamma'), delimiter);
	lambda = dlmread(strcat(testFolder,'lambda'), delimiter);
	X = reshape(dlmread(strcat(testFolder,'X'), delimiter), n, p);
	Y = reshape(dlmread(strcat(testFolder,'Y'), delimiter), n, m);
	tau = dlmread(strcat(testFolder,'tau'), delimiter);

	%run EMGLLF.m
	[phi,rho,pi,LLF,S] = EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,lambda,X,Y,tau);

	%get all stored outputs
	ref_phi = reshape(dlmread(strcat(testFolder,'phi'), delimiter), p, m, k);
	ref_rho = reshape(dlmread(strcat(testFolder,'rho'), delimiter), m, m, k);
	ref_pi = dlmread(strcat(testFolder,'pi'), delimiter);
	ref_LLF = dlmread(strcat(testFolder,'LLF'), delimiter);
	ref_S = reshape(dlmread(strcat(testFolder,'S'), delimiter), p, m, k);

	%check that output correspond to stored output
	tol = 1e-5;
	checkOutput('phi',phi,ref_phi,tol);
	checkOutput('rho',rho,ref_rho,tol);
	checkOutput('pi',pi,ref_pi,tol);
	checkOutput('LLF',LLF,ref_LLF,tol);
	checkOutput('S',S,ref_S,tol);

end
