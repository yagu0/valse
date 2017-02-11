function[phi,llh] = constructionModelesLassoRank(Pi,Rho,mini,maxi,X,Y,tau,A1,rangmin,rangmax)

	PI = 4.0 * atan(1.0);

	%get matrix sizes
	[n,p] = size(X);
	[~,m,k,~] = size(Rho);
	L = size(A1, 2); %A1 est p x m+1 x L ou p x L ?!

	%On cherche les rangs possiblement intéressants
	deltaRank = rangmax - rangmin + 1;
	Size = deltaRank^k;
	Rank = zeros(Size,k,'int64');
	for r=1:k
		%On veut le tableau de toutes les combinaisons de rangs possibles
		%Dans la première colonne : on répète (rangmax-rangmin)^(k-1) chaque chiffre : ca remplit la colonne
		%Dans la deuxieme : on répète (rangmax-rangmin)^(k-2) chaque chiffre, et on fait ca (rangmax-rangmin)^2 fois
		%...
		%Dans la dernière, on répète chaque chiffre une fois, et on fait ca (rangmin-rangmax)^(k-1) fois.
		Rank(:,r) = rangmin + reshape(repmat(0:(deltaRank-1), deltaRank^(k-r), deltaRank^(r-1)), Size, 1);
	end

	%output parameters
	phi = zeros(p,m,k,L*Size);
	llh = zeros(L*Size,2);
	for lambdaIndex=1:L
		%On ne garde que les colonnes actives
		%active sera l'ensemble des variables informatives
		active = A1(:,lambdaIndex);
		active(active==0) = [];
		if length(active) > 0
			for j=1:Size
				[phiLambda,LLF] = EMGrank(Pi(:,lambdaIndex),Rho(:,:,:,lambdaIndex),mini,maxi,X(:,active),Y,tau,Rank(j,:));
				llh((lambdaIndex-1)*Size+j,:) = [LLF, sum(Rank(j,:) .* (length(active)-Rank(j,:)+m))];
				phi(active,:,:,(lambdaIndex-1)*Size+j) = phiLambda;
			end
		end
	end

end
