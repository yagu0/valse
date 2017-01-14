function[] = generateRunSaveTest_constructionModelesLassoMLE(n, p, m, k, mini, maxi, gamma, glambda, varargin)

	%set defaults for optional inputs
	optargs = {200 15 10 3 5 10 1.0 [0.0,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.85,0.99]};
	%replace defaults by user parameters
	optargs(1:length(varargin)) = varargin;
	[n, p, m, k, mini, maxi, gamma, glambda] = optargs{:};
	tau = 1e-6;
	seuil = 1e-15;
	mini = int64(mini);
	maxi = int64(maxi);
	L = length(glambda);

	%Generate phiInit,piInit,...
	[phiInit,rhoInit,piInit,gamInit] = basicInitParameters(n, p, m, k);

	%Generate X and Y
	[X, Y, ~] = generateIOdefault(n, p, m, k);

	A2 = zeros(p,m+1,L,'int64');
    for i=1:L
        A2(:,1,i) = 1:p;
        A2(1:5,2,i) = 1:5;
    end
	A1 = zeros(p,m+1,L,'int64');
    for i=1:L
        A1(:,1,i) = 1:p;
        A1(1:5,2,i) = 1:5;
    end

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
	dlmwrite(strcat(testFolder,'glambda'), glambda, delimiter);
	dlmwrite(strcat(testFolder,'X'), reshape(X,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'Y'), reshape(Y,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'seuil'), seuil, delimiter);
	dlmwrite(strcat(testFolder,'tau'), tau, delimiter);
	dlmwrite(strcat(testFolder,'A1'), reshape(A1,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'A2'), reshape(A2,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'dimensions'), [n,p,m,k,L], delimiter);

	[phi,rho,pi,lvraisemblance] = constructionModelesLassoMLE(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,A1,A2);

	%save output
	dlmwrite(strcat(testFolder,'phi'), reshape(phi,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'rho'), reshape(rho,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'pi'), reshape(pi,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'lvraisemblance'), reshape(lvraisemblance,1,[]), delimiter);

end
