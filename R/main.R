SelMix = setRefClass(
	Class = "selmix",

	fields = c(
		# User defined

		# regression data (size n*p, where n is the number of observations,
		# and p is the number of regressors)
		X = "numeric",
		# response data (size n*m, where n is the number of observations, 
		# and m is the number of responses)
		Y = "numeric",

		# Optionally user defined (some default values)

		# power in the penalty
		gamma = "double",
		# minimum number of iterations for EM algorithm
		mini = "integer",
		# maximum number of iterations for EM algorithm
		maxi = "integer",
		# threshold for stopping EM algorithm
		eps = "double",
		# minimum number of components in the mixture
		kmin = "integer",
		# maximum number of components in the mixture
		kmax = "integer",
		rangmin = "integer",
		rangmax = "integer",
		
		# Computed through the workflow

		# initialisation for the reparametrized conditional mean parameter
		phiInit,
		# initialisation for the reparametrized variance parameter
		rhoInit,
		# initialisation for the proportions
		piInit,
		# initialisation for the allocations probabilities in each component
		tauInit,
		# values for the regularization parameter grid
		gridLambda = [];
		# je ne crois pas vraiment qu'il faille les mettre en sortie, d'autant plus qu'on construit
		# une matrice A1 et A2 pour chaque k, et elles sont grandes, donc ca coute un peu cher ...
		A1,
		A2,
		# collection of estimations for the reparametrized conditional mean parameters
		Phi,
		# collection of estimations for the reparametrized variance parameters
		Rho,
		# collection of estimations for the proportions parameters
		Pi,

		#immutable
		seuil = 1e-15;
	),

	methods = list(
		#######################
		#initialize main object
		#######################
		initialize = function(X,Y,...)
		{
			"Initialize SelMix object"

			callSuper(...)

			X <<- X;
			Y <<- Y;
			gamma <<- ifelse (hasArg("gamma"), gamma, 1.)
			mini <<- ifelse (hasArg("mini"), mini, as.integer(5))
			maxi <<- ifelse (hasArg("maxi"), maxi, as.integer(10))
			eps <<- ifelse (hasArg("eps"), eps, 1e-6)
			kmin <<- ifelse (hasArg("kmin"), kmin, as.integer(2))
			kmax <<- ifelse (hasArg("kmax"), kmax, as.integer(3))
			rangmin <<- ifelse (hasArg("rangmin"), rangmin, as.integer(2))
			rangmax <<- ifelse (hasArg("rangmax"), rangmax, as.integer(3))
		},

		##################################
		#core workflow: compute all models
		##################################

		initParameters = function(k)
		{
			"Parameters initialization"

			#smallEM initializes parameters by k-means and regression model in each component,
			#doing this 20 times, and keeping the values maximizing the likelihood after 10
			#iterations of the EM algorithm.
			init = initSmallEM(k,sx.X,sx.Y,sx.eps);
			phiInit <<- init$phi0;
			rhoInit <<- init$rho0;
			piInit	<<- init$pi0;
			tauInit <<- init$tau0;
		},

		computeGridLambda = function()
		{
			"computation of the regularization grid"
			#(according to explicit formula given by EM algorithm)

			gridLambda <<- grillelambda(sx.phiInit,sx.rhoInit,sx.piInit,sx.tauInit,sx.X,sx.Y,
				sx.gamma,sx.mini,sx.maxi,sx.eps);
		},

		computeRelevantParameters = function()
		{
			"Compute relevant parameters"

			#select variables according to each regularization parameter
			#from the grid: sx.A1 corresponding to selected variables, and
			#sx.A2 corresponding to unselected variables.
			params = selectiontotale(sx.phiInit,sx.rhoInit,sx.piInit,sx.tauInit,sx.mini,sx.maxi,
				sx.gamma,sx.gridLambda,sx.X,sx.Y,sx.seuil,sx.eps);
			A1 <<- params$A1
			A2 <<- params$A2
			Rho <<- params$Rho
			Pi <<- params$Pi
		},

		runProcedure1 = function()
		{
			"Run procedure 1 [EMGLLF]"

			#compute parameter estimations, with the Maximum Likelihood
			#Estimator, restricted on selected variables.
			res = constructionModelesLassoMLE(sx.phiInit,sx.rhoInit,sx.piInit,sx.tauInit,
				sx.mini,sx.maxi,sx.gamma,sx.gridLambda,sx.X,sx.Y,sx.seuil,sx.eps,sx.A1,sx.A2);
			return (list( phi=res$phi, rho=res$rho, pi=res$pi))
		},

		runProcedure2 = function()
		{
			"Run procedure 2 [EMGrank]"

			#compute parameter estimations, with the Low Rank
			#Estimator, restricted on selected variables.
			return (constructionModelesLassoRank(sx.Pi,sx.Rho,sx.mini,sx.maxi,sx.X,sx.Y,sx.eps,
				sx.A1,sx.rangmin,sx.rangmax)$phi)
		},

		run = function(procedure)
		{
			"main loop: over all k and all lambda"

			# Run the all procedure, 1 with the
			#maximum likelihood refitting, and 2 with the Low Rank refitting.
				p = dim(phiInit)[1]
				m = dim(phiInit)[2]
				for (k in kmin:kmax)
				{
					print(k)
					initParameters(k)
					computeGridLambda()
					computeRelevantParameters()
					if (procedure == 1)
					{
						r1 = runProcedure1(sx)
						Phi2 = Phi
						Rho2 = Rho
						Pi2 = Pi
						p = ncol(X)
						m = ncol(Y)
						if size(Phi2) == 0 #TODO: continue translation MATLAB --> R
						sx.Phi(:,:,1:k,:) = r1$phi;
						sx.Rho(:,:,1:k,:) = r1$rho;
						sx.Pi(1:k,:) = r1$pi;
						else
						sx.Phi = zeros(p,m,sx.kmax,size(Phi2,4)+size(r1$phi,4));
						sx.Phi(:,:,1:size(Phi2,3),1:size(Phi2,4)) = Phi2;
						sx.Phi(:,:,1:k,size(Phi2,4)+1:end) = r1$phi;
						sx.Rho = zeros(m,m,sx.kmax,size(Rho2,4)+size(r1$rho,4));
						sx.Rho(:,:,1:size(Rho2,3),1:size(Rho2,4)) = Rho2;
						sx.Rho(:,:,1:k,size(Rho2,4)+1:end) = r1$rho;
						sx.Pi = zeros(sx.kmax,size(Pi2,2)+size(r1$pi,2));
						sx.Pi(1:size(Pi2,1),1:size(Pi2,2)) = Pi2;
						sx.Pi(1:k,size(Pi2,2)+1:end) = r1$pi;
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
		
		##################################################
		#pruning: select only one (or a few best ?!) model
		##################################################
		#
		# 		function[model] selectModel(sx)
		# 			#TODO
		# 			#model = sxModel(...);
		# 		end
		
		)
)
