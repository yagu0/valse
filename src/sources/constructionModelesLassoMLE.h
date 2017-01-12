#ifndef valse_constructionModelesLassoMLE_H
#define valse_constructionModelesLassoMLE_H

void constructionModelesLassoMLE_core(
	// IN parameters
	const float* phiInit,
	const float* rhoInit,
	const float* piInit,
	const float* gamInit,
	int mini,
	int maxi,
	float gamma,
	const float* glambda,
	const float* X,
	const float* Y,
	float seuil,
	float tau,
	const int* A1,
	const int* A2,
	// OUT parameters
	float* phi,
	float* rho,
	float* pi,
	float* lvraisemblance,
	// additional size parameters
	int n,
	int p,
	int m,
	int k,
	int L);

#endif
