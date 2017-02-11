#ifndef valse_constructionModelesLassoMLE_H
#define valse_constructionModelesLassoMLE_H

#include "utils.h"

void constructionModelesLassoMLE_core(
	// IN parameters
	const Real* phiInit,
	const Real* rhoInit,
	const Real* piInit,
	const Real* gamInit,
	int mini,
	int maxi,
	Real gamma,
	const Real* glambda,
	const Real* X,
	const Real* Y,
	Real seuil,
	Real tau,
	const int* A1,
	const int* A2,
	// OUT parameters
	Real* phi,
	Real* rho,
	Real* pi,
	Real* llh,
	// additional size parameters
	int n,
	int p,
	int m,
	int k,
	int L);

#endif
