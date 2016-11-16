%Selection de modèle dans la procedure Lasso Rank
vraisLassoRank=vraisLassoRank';
% Selection de modele : voici les parametres selectionnes pour chaque
% dimension, et l'indice associe
[indiceLassoRank,D1LassoRank]=selectionmodele(vraisLassoRank);
tableauLassoRank=enregistrerdonnees(D1LassoRank,vraisLassoRank(indiceLassoRank,:),n);
%On veut tracer la courbe des pentes
%figure(2)
%plot(D1LassoRank/n,vraisLassoRank(indiceLassoRank),'.')
save tableauLassoRank.mat tableauLassoRank
%On conclut avec Capushe !
