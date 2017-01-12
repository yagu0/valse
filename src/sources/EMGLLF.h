#ifndef valse_EMGLLF_H
#define valse_EMGLLF_H

void EMGLLF_core(
	// IN parameters
	const float* phiInit,
	const float* rhoInit,
	const float* piInit,
	const float* gamInit,
	int mini,
	int maxi,
	float gamma,
	float lambda,
	const float* X,
	const float* Y,
	float tau,
	// OUT parameters
	float* phi,
	float* rho,
	float* pi,
	float* LLF,
	float* S,
	// additional size parameters
	int n,
	int p,
	int m,
	int k);

#endif
