ondelettes
niveauColores
figure(2)
subplot(6,1,1)
plot(1:7,courbe(1:7),'b','LineWidth',2)
hold on
plot(8:48,courbe(8:48),'b','LineWidth',2)
set(gca,'ytick',[])
ylabel('z','FontSize', 30)
set(gca, 'FontSize', 20)
axis([0 48 -1.6 1.6])
set(B, 'Position', [.91 .11 .03 .8150])
subplot(6,1,2:6)
colormap([cmapVert([1:1:end],:);cmapBleu(end:-1:1,:)]);
        %colormap(gr(end:-1:1,:));
indice = [3,3,6,12,24];
indice2=[0,cumsum(indice)];
beta = zeros(5,48);
for j=1:indice(1)
        for l=1:2^4
               beta(1,(j-1)*2^4+l) = C(indice2(1)+j);
        end
end
for r=2:5
    for j=1:indice(r)
        for l=1:2^(5-(r-1))
               beta(r,(j-1)*2^(5-(r-1))+l) = C(indice2(r)+j);
        end
    end
end


imagesc(beta)
ylabel('d_1               d_2               d_3               d_4               a_4','FontSize', 30)
set(gca, 'FontSize', 20)
%gr = gray(64);



set(gca,'xtick',[],'ytick',[])
set(gca, 'FontSize', 20)

% axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
% c=colorbar ('FontSize',18);
B=colorbar;
set(B, 'Position', [.91 .11 .03 .8150])
%set(gca,'Position',[.1 0.1 0.9 0.9])
set(gca, 'FontSize', 15)
% figure(3)    
% subplot(6,1,1)
% plot(courbe,'LineWidth',2)
% axis([0 48 -0.7 0.7])
% subplot(6,1,2)
% plot([8,24,40],C(1:3),'+','LineWidth',2)
% axis([0 48 -1.6 1.6])
% subplot(6,1,3)
% plot([8:16:48],C(4:6),'+','LineWidth',2)
% axis([0 48 -0.6 0.6])
% subplot(6,1,4)
% plot([4:8:48],C(7:12),'+','LineWidth',2)
% axis([0 48 -0.4 0.4])
% subplot(6,1,5)
% plot([2:4:48],C(13:24),'+','LineWidth',2)
% axis([0 48 -0.2 0.2])
% subplot(6,1,6)
% plot([1:2:48],C(25:48),'+','LineWidth',2)
% axis([0 48 -0.2 0.2])