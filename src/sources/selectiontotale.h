#ifndef valse_selectiontotale_H
#define valse_selectiontotale_H

// Main job on raw inputs (after transformation from mxArray)
void selectiontotale(
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
	// OUT parameters
	int* A1,
	int* A2,
	double* Rho,
	double* Pi,
	// additional size parameters
	int n,
	int p,
	int m,
	int k,
	int L);

#endif
