%Dessins prétraitements
script

XC=X2(:,1)-mean(X2(:,1));

[C1,L1] = wavedec(X2(:,1),4,'haar');
xProj1=C(4:12);

[C2,L2]=wavedec(XC',4,'haar');
xProj2=C(1:12);

Xrecon1 = waverec([zeros(1,3),xProj1',zeros(1,36)],L1,'haar');
Xrecon2 = waverec([xProj2',zeros(1,36)],L2,'haar');

plot(X2(:,1),'LineWidth',2)
hold on
plot(Xrecon1,'r-o','LineWidth',2)
plot(Xrecon2,'g-s','LineWidth',2)
xlabel('Instant of the day','FontSize', 30)
ylabel(['Load consumption and its reconstructions' sprintf('\n') ' after projecting and preprocessing'],'FontSize', 30)
set(gca, 'FontSize', 20, 'fontName','Times');
legend('Original signal','Reconstruction after projection and preprocessing 1','Reconstruction after projection and preprocessing 2','Location','NorthWest')