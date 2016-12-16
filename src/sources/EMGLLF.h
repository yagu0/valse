#ifndef valse_EMGLLF_H
#define valse_EMGLLF_H

void EMGLLF(
	// IN parameters
	const double* phiInit,
	const double* rhoInit,
	const double* piInit,
	const double* gamInit,
	int mini,
	int maxi,
	double gamma,
	double lambda,
	const double* X,
	const double* Y,
	double tau,
	// OUT parameters
	double* phi,
	double* rho,
	double* pi,
	double* LLF,
	double* S,
	// additional size parameters
	int n,
	int p,
	int m,
	int k);

#endif
