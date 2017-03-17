ind=indiceLassoMLE([2,3,4,160,165])
cpt=0;
for L2=ind
    cpt = cpt+1
%    figure(L2)
%    hold on
%    for r=1:max(b(:,L2))
%        subplot(2,ceil(max(b(:,L2))/2),r)
%        plot(X2Final(:,b(:,L2)==r),c(r))
%     end
%     hold off
    figure(1)
    subplot(2,3,cpt)
    hold on
    for r=1:max(b(:,L2))
        plot([mean(X2Final(:,b(:,L2)==r),2);mean(Y2Final(:,(b(:,L2)==r)),2)],c(r))
        title('signal moyen de chaque classe, pour chaque modele')
    end
    hold off
    
end