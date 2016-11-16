#ifndef select_EMGLLF_H
#define select_EMGLLF_H

#include "ioutils.h"

void EMGLLF(
	// IN parameters
	const Real* phiInit,
	const Real* rhoInit,
	const Real* piInit,
	const Real* gamInit,
	Int mini,
	Int maxi,
	Real gamma,
	Real lambda,
	const Real* X,
	const Real* Y,
	Real tau,
	// OUT parameters
	Real* phi,
	Real* rho,
	Real* pi,
	Real* LLF,
	Real* S,
	// additional size parameters
	mwSize n, 
	mwSize p, 
	mwSize m, 
	mwSize k);

#endif
