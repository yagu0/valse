function [indice,D1]=selectionmodele(vraisemblance)

	D=vraisemblance(:,2);
	[D1]=unique(D);
	indice=ones(1,length(D1));
	%On ne sélectionne que celui qui maximise : l'EMV
	if length(D1)>2
		for i=1:length(D1)
			a=[];
			for j=1:length(D)
				if D(j)==D1(i)
					a=[a,vraisemblance(j,1)];
				end
			end
			b=max(a);
			indice(i)=find(vraisemblance(:,1)==b,1);
		end
	end

end
