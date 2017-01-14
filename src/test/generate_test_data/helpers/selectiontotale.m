function[A1,A2,Rho,Pi] = selectiontotale(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau)

	[p,m,k] = size(phiInit);
	L = length(glambda);
	A1 = zeros(p,m+1,L,'int64');
	A2 = zeros(p,m+1,L,'int64');
	Rho = zeros(m,m,k,L);
	Pi = zeros(k,L);

	%Pour chaque lambda de la grille, on calcule les coefficients
	for lambdaIndex=1:L
		[phi,rho,pi,~,~] = EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda(lambdaIndex),X,Y,tau);

		%Si un des coefficients est supérieur au seuil, on garde cette variable
		selectedVariables = zeros(p,m);
		discardedVariables = zeros(p,m);
		atLeastOneSelectedVariable = false;
		for j=1:p
			cpt=1;
			cpt2=1;
			for mm=1:m
				if max(abs(phi(j,mm,:))) > seuil
					selectedVariables(j,cpt) = mm;
					cpt = cpt+1;
					atLeastOneSelectedVariable = true;
				else
					discardedVariables(j,cpt2) = mm;
					cpt2 = cpt2+1;
				end
			end
		end

		%Si aucun des coefficients n'a été gardé on renvoit la matrice nulle
		%Et si on enlevait ces colonnes de zéro ??? Indices des colonnes vides
		if atLeastOneSelectedVariable
			vec = [];
			for j=1:p
				if selectedVariables(j,1) ~= 0
					vec = [vec;j];
				end
			end

			%Sinon on renvoit les numéros des coefficients utiles
			A1(:,1,lambdaIndex) = [vec;zeros(p-length(vec),1)];
			A1(1:length(vec),2:m+1,lambdaIndex) = selectedVariables(vec,:);
			A2(:,1,lambdaIndex) = 1:p;
			A2(:,2:m+1,lambdaIndex) = discardedVariables;
			Rho(:,:,:,lambdaIndex) = rho;
			Pi(:,lambdaIndex) = pi;
		end

	end

end
