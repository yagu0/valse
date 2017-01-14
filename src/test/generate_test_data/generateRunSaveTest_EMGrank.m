function[] = generateRunSaveTest_EMGrank(n, p, m, k, mini, maxi, gamma, rank, varargin)

	%set defaults for optional inputs
	optargs = {200 15 10 3 5 10 1.0 1:3};
	%replace defaults by user parameters
	optargs(1:length(varargin)) = varargin;
	[n, p, m, k, mini, maxi, gamma, rank] = optargs{:};
	mini = int64(mini);
	maxi = int64(maxi);
	rank = int64(rank);
	tau = 1e-6;

	Pi = (1.0/k)*ones(1,k);
	Rho = zeros(m,m,k);
	for r=1:k
		Rho(:,:,r) = eye(m);
	end

	%Generate X and Y
	[X, Y, ~] = generateIOdefault(n, p, m, k);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	testFolder = 'data/';
	mkdir(testFolder);
	delimiter = ' ';

	%save inputs
	dlmwrite(strcat(testFolder,'Rho'), reshape(Rho,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'Pi'), Pi, delimiter);
	dlmwrite(strcat(testFolder,'mini'), mini, delimiter);
	dlmwrite(strcat(testFolder,'maxi'), maxi, delimiter);
	dlmwrite(strcat(testFolder,'X'), reshape(X,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'Y'), reshape(Y,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'tau'), tau, delimiter);
	dlmwrite(strcat(testFolder,'rank'), rank, delimiter);
	dlmwrite(strcat(testFolder,'dimensions'), [n,p,m,k], delimiter);;

	[phi,LLF] = EMGrank(Pi,Rho,mini,maxi,X,Y,tau,rank);

	%save output
	dlmwrite(strcat(testFolder,'phi'), reshape(phi,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'LLF'), LLF, delimiter);

end
