ychap=zeros(n,m,k,94);
%Ychap=zeros(n,m,l,94);

for LL=1:10
    for i=1:n
        for r=1:2
            ychap(i,:,r,LL)=x(i,:)*(phitrue(:,:,r,LL)*inv(rhotrue(:,:,r,LL)))';
            Ychap(i,:,r,LL)=waverec(ychap(i,:,r,LL)',L,'sym4');
        end
    end
end
for i=1:n
    for r=1:2
        for LL=1:10
            RMSE(i,LL,r)=sqrt(sum((donneesC(78:96,i)-Ychap(i,:,r,LL)').^2));
        end
    end
end
LL=88;
plot(Ychap(89,:,2,LL),'r')
hold on
plot(donnees(49:96,89)/100,'b')
hold off


for LL=1:10
    Z(LL,:)=principeMAP(Y,X,phiLassoMLE(:,:,:,LL),rhoLassoMLE(:,:,:,LL),piLassoMLE(:,LL),3.14);
    sum((Z(1:100)==1)+(Z(101:200)==2))
end