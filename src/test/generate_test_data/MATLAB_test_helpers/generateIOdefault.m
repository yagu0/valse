%call generateIO with default parameters (random means, covariances = identity, equirepartition)
function[X,Y,Z] = generateIOdefault(n, p, m, k)

	rangeX = 100;
	meanX = rangeX * (1 - 2*rand(k, p));
	covX = zeros(p,p,k);
	covY = zeros(m,m,k);
	for r=1:k
		covX(:,:,r) = eye(p);
		covY(:,:,r) = eye(m);
	end
	pi = (1/k) * ones(1,k);
	
	%initialize beta to a random number of non-zero random value
	beta = zeros(p,m,k);
	for j=1:p
		nonZeroCount = ceil(m*rand(1));
		beta(j,1:nonZeroCount,:) = rand(nonZeroCount, k);
	end
	
	[X,Y,Z] = generateIO(meanX, covX, covY, pi, beta, n);
	
end
