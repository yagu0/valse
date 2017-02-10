#ifndef valse_selectiontotale_H
#define valse_selectiontotale_H

#include "utils.h"

// Main job on raw inputs (after transformation from mxArray)
void selectiontotale_core(
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
	// OUT parameters
	int* A1,
	int* A2,
	Real* Rho,
	Real* Pi,
	// additional size parameters
	int n,
	int p,
	int m,
	int k,
	int L);

#endif
