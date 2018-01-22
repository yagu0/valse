%Saut de dimension
figure(1)
vec=[0,10:100:1210];

for r=vec
    tableauBis(:,:,r+1) = tableauLassoMLE;
    [aaaa,bbbb]= min(tableauLassoMLE(:,4)+r*tableauLassoMLE(:,2));
    tableauBis(:,4,r+1) = tableauLassoMLE(:,4)+r*tableauLassoMLE(:,2);
    reponse(r+1)=tableauLassoMLE(bbbb,3);
end

[x1,y1] = stairs(vec,reponse(vec+1))
plot(x1,y1,'LineWidth',2)
xlabel('$\kappa$','Interpreter','LaTex','FontSize', 45)
ylabel('Model dimension','FontSize', 45)
set(gca, 'FontSize', 40, 'fontName','Times');
[Cx1,IA,IC] = unique(y1,'stable')
hold on
plot([x1(IA(2))-0.001,x1(IA(2))],[0,Cx1(2)],'--','LineWidth',2)
plot([2*x1(IA(2))-0.001,2*x1(IA(2))],[0,Cx1(3)],'--','LineWidth',2)
%plot([0,2*x1(IA(2))],[Cx1(3),Cx1(3)],'--','LineWidth',2)
set(gca,'Box','off')
set(gca,'YTick',unique([0:800:700]))
set(gca,'XTick',unique([0:1300:1200]))
text(x1(IA(2))-20,-20,'$\hat{\kappa}$','Interpreter','LaTex','FontSize',45)
text(2*x1(IA(2))-20,-20,'$2 \hat{\kappa}$','Interpreter','LaTex','FontSize',45)
%text(-80,Cx1(3),'$\hat{m}$','Interpreter','LaTex','FontSize',45)
axis([0 1200 0 650])
hold off

%LLF pénalisé

% figure(2)
% plot(tableauLassoMLE(:,3),tableauLassoMLE(:,4)+0*tableauLassoMLE(:,2),'.','markersize', 20)
% xlabel('Model dimension','FontSize', 30)
% ylabel('Log-likelihood','Interpreter','LaTex','FontSize', 30)
% set(gca, 'FontSize', 20, 'fontName','Times');
figure(3)
ind = find(tableauLassoMLE(:,4) >4000);
tableauLassoMLE(ind,:)=[];
plot(tableauLassoMLE(2:end,3),tableauLassoMLE(2:end,4)+1100*tableauLassoMLE(2:end,2),'.','markersize', 30)
indiceInt = find(tableauLassoMLE(:,4)+1100*tableauLassoMLE(:,2)<-4401);
indiceBis = find(tableauLassoMLE(:,3)==53);
hold on
plot(tableauLassoMLE(indiceInt,3),tableauLassoMLE(indiceInt,4)+1100*tableauLassoMLE(indiceInt,2),'sr','markersize', 30,'linewidth',5)
plot(tableauLassoMLE(indiceBis,3),tableauLassoMLE(indiceBis,4)+1100*tableauLassoMLE(indiceBis,2),'dg','markersize', 30,'linewidth',5)
xlabel('Model dimension','Interpreter','LaTex','FontSize', 45)
ylabel('Penalized log-likelihood for $2 \hat{\kappa}$','Interpreter','LaTex','FontSize', 45)
set(gca, 'FontSize', 30, 'fontName','Times');
set(gcf,'Units','normal')
set(gca,'Position',[.1 .1 .88 .85])
hold off
% figure(4)
% plot(tableauLassoMLE(:,3),tableauLassoMLE(:,4)+3000*tableauLassoMLE(:,2),'.','markersize', 20)
% xlabel('Model dimension','FontSize', 30)
% ylabel('Penalized log-likelihood for $\kappa$ = 3000','Interpreter','LaTex','FontSize', 30)
% set(gca, 'FontSize', 20, 'fontName','Times');
% 
% figure(5)
% n2 = find(piLassoMLE(3,:)~=0,1)-1;
% n3 = find(piLassoMLE(4,:)~=0,1)-1;
% plot(tableauLassoMLE(:,3),tableauLassoMLE(:,4)+0*tableauLassoMLE(:,2),'.','markersize', 20)
% xlabel('Model dimension','FontSize', 30)
% ylabel('Log-likelihood','Interpreter','LaTex','FontSize', 30)
% set(gca, 'FontSize', 20, 'fontName','Times');