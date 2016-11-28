function[B1,ind,rho,pi]=suppressionmodelesegaux2(B1,rho,pi)

	ind=[];
	nombreLambda=size(B1,2);
	for l=1:nombreLambda
		for ll=1:l-1
			if B1(:,l)==B1(:,ll)
				ind=[ind l];
			end
		end
	end
	ind=unique(ind);
	B1(:,ind)=[];
	rho(:,:,:,ind)=[];
	pi(:,ind)=[];

end
