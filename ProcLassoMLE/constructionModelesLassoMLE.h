#ifndef select_constructionModelesLassoMLE_H
#define select_constructionModelesLassoMLE_H

#include "ioutils.h"

void constructionModelesLassoMLE(
	// IN parameters 
	const Real* phiInit,
	const Real* rhoInit,
	const Real* piInit,
	const Real* gamInit,
	Int mini,
	Int maxi,
	Real gamma,
	const Real* glambda,
	const Real* X,
	const Real* Y,
	Real seuil,
	Real tau,
	const Int* A1,
	const Int* A2,
	// OUT parameters
	Real* phi,
    Real* rho,
	Real* pi,
    Real* lvraisemblance,
	// additional size parameters
	mwSize n,
	mwSize p,
	mwSize m,
	mwSize k,
	mwSize L);

#endif
