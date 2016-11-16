#ifndef select_selectiontotale_H
#define select_selectiontotale_H

#include "ioutils.h"

// Main job on raw inputs (after transformation from mxArray)
void selectiontotale(
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
	// OUT parameters
	Int* A1,
	Int* A2,
	Real* Rho,
	Real* Pi,
	// additional size parameters
	mwSize n,    
	mwSize p,
	mwSize m,
	mwSize k,
	mwSize L);
		
#endif
