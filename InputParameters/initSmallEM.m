function[phiInit,rhoInit,piInit,gamInit] = initSmallEM(k,x,y,tau)
[n,m]=size(y);
gamInit1=zeros(n,k,20);
for repet=1:20
    Zinit1(:,repet)=clusterdata(y,k);
    for r=1:k
        betaInit1(:,:,r,repet)=pinv(transpose(x(Zinit1(:,repet)==r,:))*x(Zinit1(:,repet)==r,:))*transpose(x(Zinit1(:,repet)==r,:))*y(Zinit1(:,repet)==r,:);
        sigmaInit1(:,:,r,repet)=eye(m);
        phiInit1(:,:,r,repet)=betaInit1(:,:,r,repet)/sigmaInit1(:,:,r,repet);
        rhoInit1(:,:,r,repet)=inv(sigmaInit1(:,:,r,repet));
        piInit1(repet,r)=sum(Zinit1(:,repet)==r)/n;
    end
    for i=1:n
        for r=1:k
            dotProduct = (y(i,:)*rhoInit1(:,:,r,repet)-x(i,:)*phiInit1(:,:,r,repet)) * transpose(y(i,:)*rhoInit1(:,:,r,repet)-x(i,:)*phiInit1(:,:,r,repet));
            Gam(i,r) = piInit1(repet,r)*det(rhoInit1(:,:,r,repet))*exp(-0.5*dotProduct);
        end
        sumGamI = sum(Gam(i,:));
        gamInit1(i,:,repet) = Gam(i,:) / sumGamI;
    end
    miniInit=int64(10);
    maxiInit=int64(11);
    
    [~,~,~,LLFEssai,~] = EMGLLF(phiInit1(:,:,:,repet),rhoInit1(:,:,:,repet),piInit1(repet,:),gamInit1(:,:,repet),miniInit,maxiInit,1,0,x,y,tau);
    LLFinit1(repet)=LLFEssai(end);
end
[~,b]=max(LLFinit1);

phiInit=phiInit1(:,:,:,b);
rhoInit=rhoInit1(:,:,:,b);
piInit=piInit1(b,:);
gamInit=gamInit1(:,:,b);
end