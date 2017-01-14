function[] = generateRunSaveTest_EMGLLF(n, p, m, k, mini, maxi, gamma, lambda, varargin)

	%set defaults for optional inputs
	optargs = {200 15 10 3 5 10 1.0 0.5};
	%replace defaults by user parameters
	optargs(1:length(varargin)) = varargin;
	[n, p, m, k, mini, maxi, gamma, lambda] = optargs{:};
	tau = 1e-6;
	mini = int64(mini);
	maxi = int64(maxi);

	%Generate phiInit,piInit,...
	[phiInit,rhoInit,piInit,gamInit] = basicInitParameters(n, p, m, k);

	%Generate X and Y
	[X, Y, ~] = generateIOdefault(n, p, m, k);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	testFolder = 'data/';
	mkdir(testFolder);
	delimiter = ' ';

	%save inputs
	dlmwrite(strcat(testFolder,'phiInit'), reshape(phiInit,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'rhoInit'), reshape(rhoInit,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'piInit'), piInit, delimiter);
	dlmwrite(strcat(testFolder,'gamInit'), reshape(gamInit,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'mini'), mini, delimiter);
	dlmwrite(strcat(testFolder,'maxi'), maxi, delimiter);
	dlmwrite(strcat(testFolder,'gamma'), gamma, delimiter);
	dlmwrite(strcat(testFolder,'lambda'), lambda, delimiter);
	dlmwrite(strcat(testFolder,'X'), reshape(X,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'Y'), reshape(Y,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'tau'), tau, delimiter);
	dlmwrite(strcat(testFolder,'dimensions'), [n,p,m,k], delimiter);

	[phi,rho,pi,LLF,S] = EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,lambda,X,Y,tau);

	%save output
	dlmwrite(strcat(testFolder,'phi'), reshape(phi,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'rho'), reshape(rho,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'pi'), pi, delimiter);
	dlmwrite(strcat(testFolder,'LLF'), LLF, delimiter);
	dlmwrite(strcat(testFolder,'S'), reshape(S,1,[]), delimiter);

end
