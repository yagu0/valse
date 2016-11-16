# model SELECTion

This code is the applied part of the PhD thesis of [Emilie Devijver](http://www.math.u-psud.fr/~devijver/).

## Description

The function selmix delivers a multivariate Gaussian mixture in regression model collection. 
According to the parameter estimation, we can compute classical model selection criterion, as BIC or AIC, or slope heuristic, using the CAPUSHE package. 
The methodology used is described in 'Model-Based Clustering for High-Dimensional Data. Application to Functional Data.', 
available at [this location](https://hal.archives-ouvertes.fr/hal-01060063)

## Arguments

Regressors, denoted by X (of size n x p) and responses, denoted by Y (of size n x q) are must-have arguments. 

Optionally, we could add

* gamma: weight power in the Lasso penalty (according to Stadler et al., $\gamma \in \{0,1/2,1\}$;
* mini: the minimum number of iterations;
* maxi: the maximum number of iterations;
* tau: the threshold for stopping EM algorithm;
* kmin and kmax: the bounds of interesting number of components,
* rangmin and rangmax: the bounds of interesting rank values.

## Usage

	objet = selmix(X,Y)
	objet.run(index)

For index=1, it computes the Lasso-MLE procedure.
For index=2, it computes the Lasso-Rank procedure.

/!\ Be careful to the current path /!\

## Values

* phiInit, rhoInit, piInit, gamInit: the initialization of the matrices phi, rho, pi and gamma,
* gridLambda: grid of regularization parameters used to select relevant variables (if kmax-kmin=0, it is, if not, it is the last grid of regularization parameters)
* A1,A2: indices of variables selected or not selected (matrices of size (p+1) x q x size(gridLambda))
* Phi,Rho,Pi: estimations of each parameter thanks to the procedure LassoMLE if compute index=1, and thanks to the procedure LassoRank if computed index=2.


## Example

	n=10;
	p=10;
	q=5;
	X=randn(n,p);
	Y=randn(n,q);

	objet=selmix(X,Y);
	objet.run(1);
	objet.run(2);
