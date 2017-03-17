niveauColores

for L=ind
    for r=1:max(b(:,L))
        betaLassoMLE(:,:,r)=inv(rhoLassoMLE(:,:,r,L)) *  phiLassoMLE(:,:,r,L);
    end
end
dessinRho = zeros(max(K),size(rhoLassoMLE,1),max(ind));
cpt=1;
gr = gray(64);
  
 for L=ind
    k=size(find(piLassoMLE(:,L)~=0),1)
    figure(1)
    %for r=1:k 
    r=1
        subplot(1,k+1,r)
        cpt=cpt+1;
        colormap(gr(end:-1:1,:));
        imagesc(abs(betaLassoMLE(:,:,r)))
        set(gca,'xtick',[],'ytick',[])
        set(gca, 'FontSize', 30)
        title(['$\hat{\underline{\beta}_1}$'],'FontSize', 45,'Interpreter','latex')
        h=colorbar;
        set(gca, 'FontSize', 30)
        r=2
        subplot(1,k+1,r)
        cpt=cpt+1;
        colormap(gr(end:-1:1,:));
        imagesc(abs(betaLassoMLE(:,:,r)))
        set(gca,'xtick',[],'ytick',[])
        set(gca, 'FontSize', 30)
        title(['$\hat{\underline{\beta}_2}$'],'FontSize', 45,'Interpreter','latex')
        h=colorbar;
        set(gca, 'FontSize', 30)
        
    %end
    subplot(1,k+1,k+1)
    gr = gray(64);
    colormap(gr(end:-1:1,:));
    imagesc(abs(betaLassoMLE(:,:,1)-betaLassoMLE(:,:,2)))
        set(gca,'xtick',[],'ytick',[])
        title('$|\hat{\underline{\beta}_1}-\hat{\underline{\beta}}_2|$','FontSize',45,'Interpreter','latex')
        h=colorbar;
        set(gca, 'FontSize', 30)
    figure(L+k+1)
    for r=1:k
        for z=1:size(rhoLassoMLE,1)
            dessinRho(r,z,L) = inv(rhoLassoMLE(z,z,r,L)^2);
        end
    end
    
    %for r=1:k
    r=1
         subplot(k,1,r)
%         if r==1
%             colormap(cmapRouge(end:-1:1,:));
%         else if r==2 
%             colormap(cmapBleu(end:-1:1,:));
%             else colormap(cmapVert(end:-1:1,:));
%             end
%         end
        imagesc(dessinRho(r,:,L))     
        set(gca,'xtick',[],'ytick',[])
        xlabh = get(gca,'XLabel');
        set(xlabh,'Position',get(xlabh,'Position') + [0 .15 0])
        title('$\hat{\Sigma}_1$','FontSize', 55,'Interpreter','latex')
        xlabel('d_4[1]          d_4[2]          d_4[3]          d_3[1]          d_3[2]          d_3[3]          d_3[4]          d_3[5]          d_3[6]','FontSize', 30)
         set(gca, 'FontSize', 45)
        set(gca,'Position',[0.01 0.1+0.5*(2-r) 0.9 0.2])  
        
        r=2
        subplot(k,1,r)
%         if r==1
%             colormap(cmapRouge(end:-1:1,:));
%         else if r==2 
%             colormap(cmapBleu(end:-1:1,:));
%             else colormap(cmapVert(end:-1:1,:));
%             end
%         end
        imagesc(dessinRho(r,:,L))     
        set(gca,'xtick',[],'ytick',[])
        xlabh = get(gca,'XLabel');
        set(xlabh,'Position',get(xlabh,'Position') + [0 .15 0])
        title('$\hat{\Sigma}_2$','FontSize', 55,'Interpreter','latex')
        xlabel('d_4[1]          d_4[2]          d_4[3]          d_3[1]          d_3[2]          d_3[3]          d_3[4]          d_3[5]          d_3[6]','FontSize', 30)
         set(gca, 'FontSize', 45)
        set(gca,'Position',[0.01 0.1+0.5*(2-r) 0.9 0.2])
    %end
    %h=colorbar;
    gr = gray(64);
    colormap(gr(end:-1:1,:));
    axes('Position', [0.14 0.1 0.8 0.8], 'Visible', 'off');
%     c=colorbar;
%     set(gca, 'FontSize', 40)
%     caxis([0 max(max(dessinRho(:,:,L)))])
end