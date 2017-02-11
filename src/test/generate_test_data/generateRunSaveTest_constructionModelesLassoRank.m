function[] = generateRunSaveTest_constructionModelesLassoRank(n, p, m, k, L, mini, maxi, gamma, rangmin, rangmax, varargin)

	%set defaults for optional inputs
	optargs = {200 15 10 3 12 5 10 1.0 3 6};
	%replace defaults by user parameters
	optargs(1:length(varargin)) = varargin;
	[n, p, m, k, L, mini, maxi, gamma, rangmin, rangmax] = optargs{:};
	mini = int64(mini);
	maxi = int64(maxi);
	rangmin = int64(rangmin);
	rangmax = int64(rangmax);
	tau = 1e-6;

	Pi = zeros(k,L);
	for l=1:L
		Pi(:,l) = (1.0/k)*ones(1,k);
	end
	Rho = zeros(m,m,k,L);
	for l=1:L
		for r=1:k
			Rho(:,:,r,l) = eye(m);
		end
	end

	%Generate X and Y
	[X, Y, ~] = generateIOdefault(n, p, m, k);

	A1 = zeros(p,L,'int64');
    for i=1:L
        A1(:,i) = 1:p;
    end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	testFolder = 'data/';
	mkdir(testFolder);
	delimiter = ' ';

	%save inputs
	dlmwrite(strcat(testFolder,'Rho'), reshape(Rho,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'Pi'), reshape(Pi,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'mini'), mini, delimiter);
	dlmwrite(strcat(testFolder,'maxi'), maxi, delimiter);
	dlmwrite(strcat(testFolder,'X'), reshape(X,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'Y'), reshape(Y,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'tau'), tau, delimiter);
	dlmwrite(strcat(testFolder,'A1'), reshape(A1,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'rangmin'), rangmin, delimiter);
	dlmwrite(strcat(testFolder,'rangmax'), rangmax, delimiter);
	dlmwrite(strcat(testFolder,'dimensions'), [n,p,m,k,L], delimiter);

	[phi,llh] = constructionModelesLassoRank(Pi,Rho,mini,maxi,X,Y,tau,A1,rangmin,rangmax);

	%save output
	Size = (rangmax-rangmin+1)^k;
	dlmwrite(strcat(testFolder,'phi'), reshape(phi,1,[]), delimiter);
	dlmwrite(strcat(testFolder,'llh'), reshape(llh,1,[]), delimiter);

end
