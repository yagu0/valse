#ifndef valse_constructionModelesLassoMLE_H
#define valse_constructionModelesLassoMLE_H

void constructionModelesLassoMLE_core(
	// IN parameters
	const double* phiInit,
	const double* rhoInit,
	const double* piInit,
	const double* gamInit,
	int mini,
	int maxi,
	double gamma,
	const double* glambda,
	const double* X,
	const double* Y,
	double seuil,
	double tau,
	const int* A1,
	const int* A2,
	// OUT parameters
	double* phi,
	double* rho,
	double* pi,
	double* lvraisemblance,
	// additional size parameters
	int n,
	int p,
	int m,
	int k,
	int L);

#endif
