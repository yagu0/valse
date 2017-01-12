#ifndef valse_selectiontotale_H
#define valse_selectiontotale_H

// Main job on raw inputs (after transformation from mxArray)
void selectiontotale_core(
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
	// OUT parameters
	int* A1,
	int* A2,
	float* Rho,
	float* Pi,
	// additional size parameters
	int n,
	int p,
	int m,
	int k,
	int L);

#endif
