## TODO: turn this code into R

classdef selmix < handle
    
    properties (SetAccess = private)
        %always user defined
        X; % regression data (size n*p, where n is the number of observations, and p is the number of regressors)
        Y; % response data (size n*m, where n is the number of observations, and m is the number of responses)
        
        %optionally user defined some default values
        gamma; % power in the penalty
        mini; % minimum number of iterations for EM algorithm
        maxi; % maximum number of iterations for EM algorithm
        eps; % threshold for stopping EM algorithm
        kmin; % minimum number of components in the mixture
        kmax; % maximum number of components in the mixture
        rangmin;
        rangmax;
        
        %computed through the workflow
        phiInit; % initialisation for the reparametrized conditional mean parameter
        rhoInit; % initialisation for the reparametrized variance parameter
        piInit; % initialisation for the proportions
        tauInit; % initialisation for the allocations probabilities in each component
        gridLambda = []; % values for the regularization parameter grid
        A1; %je ne crois pas vraiment qu'il faille les mettre en sortie, d'autant plus qu'on construit une matrice A1 et A2 pour chaque k, et elles sont grandes, donc ca coute un peu cher ...
        A2;
        Phi; % collection of estimations for the reparametrized conditional mean parameters
        Rho; % collection of estimations for the reparametrized variance parameters
        Pi; % collection of estimations for the proportions parameters
    end
    
    properties (Constant)
        %immutable
        seuil = 1e-15;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%
        %initialize main object
        %%%%%%%%%%%%%%%%%%%%%%%
        
        function sx = selmix(X,Y,varargin)
            %set defaults for optional inputs
            optargs = {1.0 5 10 1e-6 2 3 2 3};
            %replace defaults by user parameters
            optargs(1:length(varargin)) = varargin;
            sx.X = X;
            sx.Y = Y;
            [sx.gamma,sx.mini,sx.maxi,sx.eps,sx.kmin,sx.kmax,sx.rangmin,sx.rangmax] = optargs{:};
            %			z = somme(sx.X,sx.Y);
            %			sx.Z=z;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %core workflow: compute all models
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function initParameters(sx,k)
            [phi0,rho0,pi0,tau0] = initSmallEM(k,sx.X,sx.Y,sx.eps); %smallEM initializes parameters by k-means and regression model in each component, doing this 20 times, and keeping the values maximizing the likelihood after 10 iterations of the EM algorithm.
            sx.phiInit = phi0;
            sx.rhoInit = rho0;
            sx.piInit  = pi0;
            sx.tauInit = tau0;
        end
        
        function computeGridLambda(sx)
            [sx.gridLambda] = grillelambda(sx.phiInit,sx.rhoInit,sx.piInit,sx.tauInit,sx.X,sx.Y,sx.gamma,sx.mini,sx.maxi,sx.eps);
            % computation of the regularization grid, according to explicit
            % formulae given by EM algorithm.
        end
        
        function computeRelevantParameters(sx)
            [sx.A1,sx.A2,sx.Rho,sx.Pi] = selectiontotale(sx.phiInit,sx.rhoInit,sx.piInit,sx.tauInit,sx.mini,sx.maxi,sx.gamma,sx.gridLambda,sx.X,sx.Y,sx.seuil,sx.eps);
            %select variables according to each regularization parameter
            %from the grid: sx.A1 corresponding to selected variables, and
            %sx.A2 corresponding to unselected variables.
        end
        
        function [sx,phi,rho,pi]=runProcedure1(sx)
            [phi,rho,pi,~] = constructionModelesLassoMLE(sx.phiInit,sx.rhoInit,sx.piInit,sx.tauInit,sx.mini,sx.maxi,sx.gamma,sx.gridLambda,sx.X,sx.Y,sx.seuil,sx.eps,sx.A1,sx.A2);
            %compute parameter estimations, with the Maximum Likelihood
            %Estimator, restricted on selected variables.
        end
        
        function [phi] =runProcedure2(sx)
            [phi,~] = constructionModelesLassoRank(sx.Pi,sx.Rho,sx.mini,sx.maxi,sx.X,sx.Y,sx.eps,sx.A1,sx.rangmin,sx.rangmax);
            %compute parameter estimations, with the Low Rank
            %Estimator, restricted on selected variables.
        end
        
        % main loop: over all k and all lambda
        function run(sx,procedure) % Run the all procedure, 1 with the
            %maximum likelihood refitting, and 2 with the Low Rank refitting.
            [p,m,~]=size(sx.phiInit);
            for k=sx.kmin:sx.kmax
                k
                initParameters(sx,k);
                computeGridLambda(sx);
                computeRelevantParameters(sx);
                if (procedure == 1)
                    [~,phi,rho,pi] = runProcedure1(sx);
                    Phi2 = sx.Phi;
                    Rho2 = sx.Rho;
                    Pi2 = sx.Pi;
                    p = size(sx.X,2);
                    m = size(sx.Y,2);
                    if size(Phi2) == 0
                        sx.Phi(:,:,1:k,:) = phi;
                        sx.Rho(:,:,1:k,:) = rho;
                        sx.Pi(1:k,:) = pi;
                    else
                        sx.Phi = zeros(p,m,sx.kmax,size(Phi2,4)+size(phi,4));
                        sx.Phi(:,:,1:size(Phi2,3),1:size(Phi2,4)) = Phi2;
                        sx.Phi(:,:,1:k,size(Phi2,4)+1:end) = phi;
                        sx.Rho = zeros(m,m,sx.kmax,size(Rho2,4)+size(rho,4));
                        sx.Rho(:,:,1:size(Rho2,3),1:size(Rho2,4)) = Rho2;
                        sx.Rho(:,:,1:k,size(Rho2,4)+1:end) = rho;
                        sx.Pi = zeros(sx.kmax,size(Pi2,2)+size(pi,2));
                        sx.Pi(1:size(Pi2,1),1:size(Pi2,2)) = Pi2;
                        sx.Pi(1:k,size(Pi2,2)+1:end) = pi;
                    end
                else
                    [phi] = runProcedure2(sx);
                    phi
                    Phi2 = sx.Phi;
                    if size(Phi2,1) == 0
                        sx.Phi(:,:,1:k,:) = phi;
                    else
                        size(Phi2)
                        sx.Phi = zeros(p,m,sx.kmax,size(Phi2,4)+size(phi,4));
                        size(sx.Phi)
                        sx.Phi(:,:,1:size(Phi2,3),1:size(Phi2,4)) = Phi2;
                        sx.Phi(:,:,1:k,size(Phi2,4)+1:end) = phi;
                    end
                    
                end
                
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %pruning: select only one (or a few best ?!) model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 		function[model] selectModel(sx)
        % 			%TODO
        % 			%model = sxModel(...);
        % 		end
        
    end
    
end


%%%%%%%%%%%%%
%OLD VERSION:

%~ for k=K
%~ % On initialise
%~ initialisation
%~ disp('Initialisé')
%~ % On construit la grille des lambdas : variables informatives
%~ [glambda,phiEmv,rhoEmv,piEmv]=grillelambda(phiInit,rhoInit,piInit,tauInit,x,y,gamma,mini,maxi,tau);
%~ glambda=glambda(1:3);
%~ disp('glambda construite')
%~ % On trouve les variables informatives pour chaque lambda : S est la
%~ % matrice des coefficients des variables informatives
%~ % on parallelise à l interieur du selectiontotale()
%~ [B1,B2,rhoLasso,piLasso]=selectiontotale(phiInit,rhoInit,piInit,tauInit,mini,maxi,gamma,glambda,x,y,10^-15,tau);
%~ %S1 les variables informatives, S2 celles qui ne le sont pas
%~ [B1bis,B2bis,glambda2bis,ind,rhoLasso,piLasso]=suppressionmodelesegaux(B1,B2,glambda,rhoLasso,piLasso);
%~ dessinVariablesSelectionnees;
%~ disp('Le Lasso est fait')
%~ % Pour chaque lambda ainsi construit, on veut calculer l'EMV pour la procédure Lasso-MLE
%~ %On obtient une collection de modèles pour Lasso-MLE
%~ % ICI AUSSI ON PEUT PARALLELISER a l interieur de constructionModelesLassoMLE
%~ nombreLambda=size(B1bis,3);
%~ %[phiLassoMLE,rhoLassoMLE,piLassoMLE,vraisLassoMLE]=constructionModelesLassoMLE(phiInit,rhoInit,piInit, tauInit,mini,maxi,gamma,glambda2bis,X,Y,seuil,tau,B1bis,B2bis)
%~ %Pour Lasso-Rank
%~ %on ne garde que les colonnes actives
%~ B1ter=B1bis(:,1,:);
%~ %		[B1ter,ind,rhoLasso,piLasso]=suppressionmodelesegaux2(B1ter,rhoLasso,piLasso)
%~ %Pour chaque lambda, on veut calculer l'estimateur low rank pour la procédure Lasso-Rank
%~ %On obtient une collection de modèles pour Lasso-Rank
%~ %ICI AUSSI ON PEUT PARALLELISER a linterieur de constructionModelesLassoRank
%~ nombreLambda2=size(B1ter,2);
%~ [phi,lvraisemblance,Z] = constructionModelesLassoRank(...
%~ piEmv,rhoEmv,mini,maxi,X,Y,tau,B1ter,2,4);
%~ disp('On a construit les modèles pour les deux procédures')
%~ end
%~ %selectionModelesLassoMLE;
%~ %selectionModelesLassoRank;
%~ disp('On a sélectionné les modèles dans les deux procédures')
