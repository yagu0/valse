beta = betaLassoMLE(:,:,1:2,ind);
for i=1:338
    for r=1:2
        YChap(:,i,r) = X(i,:) * beta(:,:,r);
        yChap(:,i,r) = waverec([zeros(1,3),YChap(:,i,r)',zeros(1,36)],L,'haar');
    end
end
