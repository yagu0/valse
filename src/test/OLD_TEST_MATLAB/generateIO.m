%X is generated following a gaussian mixture \sum pi_r N(meanX_k, covX_r)
%Y is generated then, with Y_i ~ \sum pi_r N(Beta_r.X_i, covY_r)
function[X,Y,Z] = generateIO(meanX, covX, covY, pi, beta, n)

	[p, ~, k] = size(covX);
	[m, ~, ~] = size(covY);
	if exist('octave_config_info')
		%Octave statistics package	doesn't have gmdistribution()
		X = zeros(n, p);
		Z = zeros(n);
		cs = cumsum(pi);
		for i=1:n
			%TODO: vectorize ? http://stackoverflow.com/questions/2977497/weighted-random-numbers-in-matlab
			tmpRand01 = rand();
			[~,Z(i)] = min(cs - tmpRand01 >= 0);
			X(i,:) = mvnrnd(meanX(Z(i),:), covX(:,:,Z(i)), 1);
		end
	else
		gmDistX = gmdistribution(meanX, covX, pi);
		[X, Z] = random(gmDistX, n);
	end
	
	Y = zeros(n, m);
	BX = zeros(n,m,k);
	for i=1:n
		for r=1:k
			%compute beta_r . X_i
			BXir = zeros(1, m);
			for mm=1:m
				BXir(mm) = dot(X(i,:), beta(:,mm,r));
			end
			%add pi(r) * N(beta_r . X_i, covY) to Y_i
			Y(i,:) = Y(i,:) + pi(r) * mvnrnd(BXir, covY(:,:,r), 1);
		end
	end

end
