function[phi,rho,pi,LLF,S] = EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,lambda,X,Y,tau)

	%Get matrices dimensions
	PI = 4.0 * atan(1.0);
	n = size(X, 1);
	[p,m,k] = size(phiInit);

	%Initialize outputs
	phi = phiInit;
	rho = rhoInit;
	pi = piInit;
	LLF = zeros(maxi,1);
	S = zeros(p,m,k);

	%Other local variables
	%NOTE: variables order is always n,p,m,k
	gam = gamInit;
	Gram2 = zeros(p,p,k);
	ps2 = zeros(p,m,k);
	b = zeros(k,1);
	pen = zeros(maxi,k);
	X2 = zeros(n,p,k);
	Y2 = zeros(n,m,k);
	dist = 0;
	dist2 = 0;
	ite = 1;
	pi2 = zeros(k,1);
	ps = zeros(m,k);
	nY2 = zeros(m,k);
	ps1 = zeros(n,m,k);
	nY21 = zeros(n,m,k);
	Gam = zeros(n,k);
	EPS = 1e-15;

	while ite<=mini || (ite<=maxi && (dist>=tau || dist2>=sqrt(tau)))

		Phi = phi;
		Rho = rho;
		Pi = pi;

		%Calculs associés à Y et X
		for r=1:k
			for mm=1:m
				Y2(:,mm,r) = sqrt(gam(:,r)) .* Y(:,mm);
			end
			for i=1:n
				X2(i,:,r) = X(i,:) .* sqrt(gam(i,r));
			end
			for mm=1:m
				ps2(:,mm,r) = transpose(X2(:,:,r)) * Y2(:,mm,r);
			end
			for j=1:p
				for s=1:p
					Gram2(j,s,r) = dot(X2(:,j,r), X2(:,s,r));
				end
			end
		end

		%%%%%%%%%%
		%Etape M %
		%%%%%%%%%%

		%Pour pi
		for r=1:k
			b(r) = sum(sum(abs(phi(:,:,r))));
		end
		gam2 = sum(gam,1);
		a = sum(gam*transpose(log(pi)));

		%tant que les proportions sont negatives
		kk = 0;
		pi2AllPositive = false;
		while ~pi2AllPositive
			pi2 = pi + 0.1^kk * ((1/n)*gam2 - pi);
			pi2AllPositive = true;
			for r=1:k
				if pi2(r) < 0
					pi2AllPositive = false;
					break;
				end
			end
			kk = kk+1;
		end

		%t(m) la plus grande valeur dans la grille O.1^k tel que ce soit
		%décroissante ou constante
		while (-1/n*a+lambda*((pi.^gamma)*b))<(-1/n*gam2*transpose(log(pi2))+lambda.*(pi2.^gamma)*b) && kk<1000
			pi2 = pi+0.1^kk*(1/n*gam2-pi);
			kk = kk+1;
		end
		t = 0.1^(kk);
		pi = (pi+t*(pi2-pi)) / sum(pi+t*(pi2-pi));

		%Pour phi et rho
		for r=1:k
			for mm=1:m
				for i=1:n
					ps1(i,mm,r) = Y2(i,mm,r) * dot(X2(i,:,r), phi(:,mm,r));
					nY21(i,mm,r) = (Y2(i,mm,r))^2;
				end
				ps(mm,r) = sum(ps1(:,mm,r));
				nY2(mm,r) = sum(nY21(:,mm,r));
				rho(mm,mm,r) = ((ps(mm,r)+sqrt(ps(mm,r)^2+4*nY2(mm,r)*(gam2(r))))/(2*nY2(mm,r)));
			end
		end
		for r=1:k
			for j=1:p
				for mm=1:m
					S(j,mm,r) = -rho(mm,mm,r)*ps2(j,mm,r) + dot(phi(1:j-1,mm,r),Gram2(j,1:j-1,r)')...
						+ dot(phi(j+1:p,mm,r),Gram2(j,j+1:p,r)');
					if abs(S(j,mm,r)) <= n*lambda*(pi(r)^gamma)
						phi(j,mm,r)=0;
					else
						if S(j,mm,r)> n*lambda*(pi(r)^gamma)
							phi(j,mm,r)=(n*lambda*(pi(r)^gamma)-S(j,mm,r))/Gram2(j,j,r);
						else
							phi(j,mm,r)=-(n*lambda*(pi(r)^gamma)+S(j,mm,r))/Gram2(j,j,r);
						end
					end
				end
			end
		end

		%%%%%%%%%%
		%Etape E %
		%%%%%%%%%%

		sumLogLLF2 = 0.0;
		for i=1:n
			%precompute dot products to numerically adjust their values
			dotProducts = zeros(k,1);
			for r=1:k
				dotProducts(r)= (Y(i,:)*rho(:,:,r)-X(i,:)*phi(:,:,r)) * transpose(Y(i,:)*rho(:,:,r)-X(i,:)*phi(:,:,r));
			end
			shift = 0.5*min(dotProducts);

			%compute Gam(:,:) using shift determined above
			sumLLF1 = 0.0;
			for r=1:k
				Gam(i,r) = pi(r)*det(rho(:,:,r))*exp(-0.5*dotProducts(r) + shift);
				sumLLF1 = sumLLF1 + Gam(i,r)/(2*PI)^(m/2);
			end
			sumLogLLF2 = sumLogLLF2 + log(sumLLF1);
			sumGamI = sum(Gam(i,:));
			if sumGamI > EPS
				gam(i,:) = Gam(i,:) / sumGamI;
			else
				gam(i,:) = zeros(k,1);
			end
		end

		sumPen = 0.0;
		for r=1:k
			sumPen = sumPen + pi(r).^gamma .* b(r);
		end
		LLF(ite) = -(1/n)*sumLogLLF2 + lambda*sumPen;

		if ite == 1
			dist = LLF(ite);
		else
			dist = (LLF(ite)-LLF(ite-1))/(1+abs(LLF(ite)));
		end

		Dist1=max(max(max((abs(phi-Phi))./(1+abs(phi)))));
		Dist2=max(max(max((abs(rho-Rho))./(1+abs(rho)))));
		Dist3=max(max((abs(pi-Pi))./(1+abs(Pi))));
		dist2=max([Dist1,Dist2,Dist3]);

		ite=ite+1;
	end

	pi = transpose(pi);

end
