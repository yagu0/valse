function[] = testSelectiontotale()

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

	%run constructionModelesLassoMLE.m
	[A1,A2,Rho,Pi] = selectiontotale(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau);

	%get all stored outputs
	ref_A1 = int64(reshape(dlmread(strcat(testFolder,'A1'), delimiter), p, m+1, L));
	ref_A2 = int64(reshape(dlmread(strcat(testFolder,'A2'), delimiter), p, m+1, L));
	ref_Rho = reshape(dlmread(strcat(testFolder,'Rho'), delimiter), m, m, k, L);
	ref_Pi = reshape(dlmread(strcat(testFolder,'Pi'), delimiter), k, L);

	%check that output correspond to stored output
	tol = 1e-5;
	checkOutput('A1',A1,ref_A1,tol);
	checkOutput('A2',A2,ref_A2,tol);
	checkOutput('Rho',Rho,ref_Rho,tol);
	checkOutput('Pi',Pi,ref_Pi,tol);

end
