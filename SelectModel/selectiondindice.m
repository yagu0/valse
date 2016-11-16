%utile dans selectiontotale.m, remarque quels coefficients sont nuls et lesquels ne le sont pas
function[A,B]=selectiondindice(phi,seuil)

	[pp,m,~]=size(phi);
	A=zeros(pp,m);
	B=zeros(pp,m);
	for j=1:pp
		cpt=0;cpt2=0;
		for mm=1:m
			if (max(phi(j,mm,:)) > seuil)
				cpt=cpt+1;
				A(j,cpt)=mm;
			else cpt2=cpt2+1;
				 B(j,cpt2)=mm;
			end
		end
	end

end 
