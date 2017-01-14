function[phi,LLF] = EMGrank(Pi,Rho,mini,maxi,X,Y,tau,rank)

	% get matrix sizes
	[~,m,k] = size(Rho);
	[n,p] = size(X);

	% allocate output matrices
	phi = zeros(p,m,k);
	Z = ones(n,1,'int64');
	LLF = 0.0;

	% local variables
	Phi = zeros(p,m,k);
	deltaPhi = 0.0;
	deltaPhi = [];
	sumDeltaPhi = 0.0;
	deltaPhiBufferSize = 20;

	%main loop (at least mini iterations)
	ite = int64(1);
	while ite<=mini || (ite<=maxi && sumDeltaPhi>tau)

		%M step: Mise à jour de Beta (et donc phi)
		for r=1:k
			if (sum(Z==r) == 0)
				continue;
			end
			%U,S,V = SVD of (t(Xr)Xr)^{-1} * t(Xr) * Yr
			[U,S,V] = svd(pinv(transpose(X(Z==r,:))*X(Z==r,:))*transpose(X(Z==r,:))*Y(Z==r,:));
			%Set m-rank(r) singular values to zero, and recompose
			%best rank(r) approximation of the initial product
			S(rank(r)+1:end,:) = 0;
			phi(:,:,r) = U * S * transpose(V) * Rho(:,:,r);
		end

		%Etape E et calcul de LLF
		sumLogLLF2 = 0.0;
		for i=1:n
			sumLLF1 = 0.0;
			maxLogGamIR = -Inf;
			for r=1:k
				dotProduct = (Y(i,:)*Rho(:,:,r)-X(i,:)*phi(:,:,r)) * transpose(Y(i,:)*Rho(:,:,r)-X(i,:)*phi(:,:,r));
				logGamIR = log(Pi(r)) + log(det(Rho(:,:,r))) - 0.5*dotProduct;
				%Z(i) = index of max (gam(i,:))
				if logGamIR > maxLogGamIR
					Z(i) = r;
					maxLogGamIR = logGamIR;
				end
				sumLLF1 = sumLLF1 + exp(logGamIR) / (2*pi)^(m/2);
			end
			sumLogLLF2 = sumLogLLF2 + log(sumLLF1);
		end

		LLF = -1/n * sumLogLLF2;

		% update distance parameter to check algorithm convergence (delta(phi, Phi))
		deltaPhi = [ deltaPhi, max(max(max((abs(phi-Phi))./(1+abs(phi))))) ];
		if length(deltaPhi) > deltaPhiBufferSize
			deltaPhi = deltaPhi(2:length(deltaPhi));
		end
		sumDeltaPhi = sum(abs(deltaPhi));

		% update other local variables
		Phi = phi;
		ite = ite+1;

	end

end
