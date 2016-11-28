%Selection de modèle dans la procédure Lasso-MLE
% Selection de modele : voici les parametres selectionnes pour chaque
% dimension, et l'indice associe
[indiceLassoMLE,D1LassoMLE] = selectionmodele(vraisLassoMLE);

%on veut tracer la courbe
%figure(2)
%plot(D1LassoMLE/n,vraisLassoMLE(indiceLassoMLE),'.')
%on veut appliquer  Capushe : il nous faut un tableau de données
tableauLassoMLE = enregistrerdonnees(D1LassoMLE,vraisLassoMLE(indiceLassoMLE,:),n);
save tableauLassoMLE.mat tableauLassoMLE
%On conclut avec Capushe !
