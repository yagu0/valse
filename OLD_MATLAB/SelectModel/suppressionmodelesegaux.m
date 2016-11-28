function[B1,B2,glambda,ind,rho,pi]=suppressionmodelesegaux(B1,B2,glambda,rho,pi)

	ind=[];
	for l=1:length(glambda) 
		for ll=1:l-1
			if B1(:,:,l)==B1(:,:,ll)
				ind=[ind l];
			end
		end
	end
	ind=unique(ind);
	B1(:,:,ind)=[];
	glambda(ind)=[];
	B2(:,:,ind)=[];
	rho(:,:,:,ind)=[];
	pi(:,ind)=[];

end
