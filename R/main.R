Valse = setRefClass(
	Class = "Valse",

	fields = c(
		# User defined

		# regression data (size n*p, where n is the number of observations,
		# and p is the number of regressors)
		X = "matrix",
		# response data (size n*m, where n is the number of observations,
		# and m is the number of responses)
		Y = "matrix",

		# Optionally user defined (some default values)

		# power in the penalty
		gamma = "numeric",
		# minimum number of iterations for EM algorithm
		mini = "integer",
		# maximum number of iterations for EM algorithm
		maxi = "integer",
		# threshold for stopping EM algorithm
		eps = "numeric",
		# minimum number of components in the mixture
		kmin = "integer",
		# maximum number of components in the mixture
		kmax = "integer",
		rangmin = "integer",
		rangmax = "integer",
		
		# Computed through the workflow

		# initialisation for the reparametrized conditional mean parameter
		phiInit = "numeric",
		# initialisation for the reparametrized variance parameter
		rhoInit = "numeric",
		# initialisation for the proportions
		piInit = "numeric",
		# initialisation for the allocations probabilities in each component
		tauInit = "numeric",
		# values for the regularization parameter grid
		gridLambda = "numeric",
		# je ne crois pas vraiment qu'il faille les mettre en sortie, d'autant plus qu'on construit
		# une matrice A1 et A2 pour chaque k, et elles sont grandes, donc ca coute un peu cher ...
		A1 = "integer",
		A2 = "integer",
		# collection of estimations for the reparametrized conditional mean parameters
		Phi = "numeric",
		# collection of estimations for the reparametrized variance parameters
		Rho = "numeric",
		# collection of estimations for the proportions parameters
		Pi = "numeric",

		#immutable (TODO:?)
		seuil = "numeric"
	),

	methods = list(
		#######################
		#initialize main object
		#######################
		initialize = function(X,Y,...)
		{
			"Initialize Valse object"

			callSuper(...)

			X <<- X
			Y <<- Y
			gamma <<- ifelse (hasArg("gamma"), gamma, 1.)
			mini <<- ifelse (hasArg("mini"), mini, as.integer(5))
			maxi <<- ifelse (hasArg("maxi"), maxi, as.integer(10))
			eps <<- ifelse (hasArg("eps"), eps, 1e-6)
			kmin <<- ifelse (hasArg("kmin"), kmin, as.integer(2))
			kmax <<- ifelse (hasArg("kmax"), kmax, as.integer(3))
			rangmin <<- ifelse (hasArg("rangmin"), rangmin, as.integer(2))
			rangmax <<- ifelse (hasArg("rangmax"), rangmax, as.integer(3))
			seuil <<- 1e-15 #immutable (TODO:?)
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
			init = initSmallEM(k,X,Y,eps)
			phiInit <<- init$phi0
			rhoInit <<- init$rho0
			piInit	<<- init$pi0
			tauInit <<- init$tau0
		},

		computeGridLambda = function()
		{
			"computation of the regularization grid"
			#(according to explicit formula given by EM algorithm)

			gridLambda <<- gridLambda(phiInit,rhoInit,piInit,tauInit,X,Y,gamma,mini,maxi,eps)
		},

		computeRelevantParameters = function()
		{
			"Compute relevant parameters"

			#select variables according to each regularization parameter
			#from the grid: A1 corresponding to selected variables, and
			#A2 corresponding to unselected variables.
			params = selectiontotale(
				phiInit,rhoInit,piInit,tauInit,mini,maxi,gamma,gridLambda,X,Y,seuil,eps)
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
			return ( constructionModelesLassoMLE(
				phiInit,rhoInit,piInit,tauInit,mini,maxi,gamma,gridLambda,X,Y,seuil,eps,A1,A2) )
		},

		runProcedure2 = function()
		{
			"Run procedure 2 [EMGrank]"

			#compute parameter estimations, with the Low Rank
			#Estimator, restricted on selected variables.
			return ( constructionModelesLassoRank(Pi,Rho,mini,maxi,X,Y,eps,
				A1,rangmin,rangmax) )
		},

		run = function()
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
					r1 = runProcedure1()
					Phi2 = Phi
					Rho2 = Rho
					Pi2 = Pi
					p = ncol(X)
					m = ncol(Y)
					if (is.null(dim(Phi2))) #test was: size(Phi2) == 0
					{
						Phi[,,1:k] <<- r1$phi
						Rho[,,1:k] <<- r1$rho
						Pi[1:k,] <<- r1$pi
					} else
					{
						Phi <<- array(0., dim=c(p,m,kmax,dim(Phi2)[4]+dim(r1$phi)[4]))
						Phi[,,1:(dim(Phi2)[3]),1:(dim(Phi2)[4])] <<- Phi2
						Phi[,,1:k,dim(Phi2)[4]+1] <<- r1$phi
						Rho <<- array(0., dim=c(m,m,kmax,dim(Rho2)[4]+dim(r1$rho)[4]))
						Rho[,,1:(dim(Rho2)[3]),1:(dim(Rho2)[4])] <<- Rho2
						Rho[,,1:k,dim(Rho2)[4]+1] <<- r1$rho
						Pi <<- array(0., dim=c(kmax,dim(Pi2)[2]+dim(r1$pi)[2]))
						Pi[1:nrow(Pi2),1:ncol(Pi2)] <<- Pi2
						Pi[1:k,ncol(Pi2)+1] <<- r1$pi
					}
				} else
				{
					phi = runProcedure2()$phi
					Phi2 = Phi
					if (dim(Phi2)[1] == 0)
					{
						Phi[,,1:k,] <<- phi
					} else
					{
						Phi <<- array(0., dim=c(p,m,kmax,dim(Phi2)[4]+dim(phi)[4]))
						Phi[,,1:(dim(Phi2)[3]),1:(dim(Phi2)[4])] <<- Phi2
						Phi[,,1:k,-(1:(dim(Phi2)[4]))] <<- phi
					}
				}
			}
		}

		##################################################
		#TODO: pruning: select only one (or a few best ?!) model
		##################################################
		#
		# 		function[model] selectModel(
		# 			#TODO
		# 			#model = odel(...)
		# 		end

		)
)
